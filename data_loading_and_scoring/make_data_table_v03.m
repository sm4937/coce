function [data,failed_data] = make_data_table_v03(raw_data)

%filename = ['\data\' filename];

raw_data = sortrows(raw_data,'trial_index'); %re-sort scrambled data
data = table;
colnames = raw_data.Properties.VariableNames;
default_length = 32;
subjnum = unique(raw_data.subjnum);
NFC = NaN; SAPS = NaN;

TOTidx = find(ismember(colnames,{'TOT'}));
overallidx = find(ismember(colnames,{'overall'}));
displayidx = find(ismember(colnames,{'task_displayed'}));
timesteps = find(ismember(colnames,{'time_elapsed'}));
rtidx = find(ismember(colnames,{'rt'}));
detectidx = find(ismember(colnames,{'detect'}));
nbackidx = find(ismember(colnames,{'nback'}));
perfidx = find(ismember(colnames,{'perf'}));
value = find(ismember(colnames,{'value'}));
offer = find(ismember(colnames,{'offer'}));
version_idx = ismember(colnames,{'exp_version'});
taskidx = find(ismember(colnames,{'task'}));
nidx = find(ismember(colnames,{'n'}));
% classic accuracy column doesn't work here since a withheld response is
% correct so let's do
raw_data.new_correct = categorical(raw_data.key_press) == categorical(raw_data.correct_key); %compare these hits/false alarms/misses
raw_data.new_correct(raw_data.correct_key==90) = true;
% so this is the whole thing, which is unfortunate because it includes
% feedback and instructions and BDM and stuff
% can be trimmed as you go
accidx = find(ismember(colnames,{'new_correct'})); 

%% extract basic performance information

%clean the detect and n-back columns up, it's messy because of the zero's in it become
%"undefined"s
detect = raw_data.detect;
valid = (detect==categorical({'1'})|detect==categorical({'NULL'}));
raw_data.detect(~valid) = categorical({'0'});
raw_data.detect = string(raw_data.detect);

nback = raw_data.nback;
valid = (nback==categorical({'1'})|nback==categorical({'NULL'}));
raw_data.nback(~valid) = categorical({'0'});
raw_data.nback = string(raw_data.nback);

%separate out the people who failed from the people who didn't 
tasknumlist = double(string(raw_data.tasknum));
tasknumlist = unique(tasknumlist(~isnan(tasknumlist)));
if sum(~isnan(tasknumlist))>0
    if tasknumlist(end) < (default_length-1)
        data.failed = true;
        data.subjnum = subjnum; %print some stuff about people who fail, eventually
        failed_data = dropout(raw_data,default_length); %pass it into a deeper function to figure out they dropped out
        return
    else
        failed = false;
        failed_data = table;
        demo = process_demo(raw_data); %we'll get this back
        data.age = demo.age; data.sex = demo.sex; data.hand = demo.hand; data.diffrating = demo.diffrating; data.blurs = demo.blurs; data.task_blurs = demo.task_blurs;
        data.session = unique(raw_data.session);
    end
else
    data.failed = true;
    data.subjnum = subjnum;
    failed_data = dropout(raw_data,default_length);
    return
end
fullscreen = false;

% pull out need for cognition scores
[NFC,SAPS,SAPSs,SAPSd] = grade_scales(raw_data);
%doesn't work for early pilots

rts = table2array(raw_data(:,rtidx));
logrts = log(rts);

overall_list = double(string(table2array(raw_data(:,overallidx)))); %overall accuracy
overall = unique(overall_list(~isnan(overall_list)));
total_length = str2num(char(table2array(raw_data(end,TOTidx))))/60000; %convert msec to minutes
task_list = table2array(raw_data(:,taskidx));
exp_version = table2array(raw_data(1,version_idx));

not_practice = double(string(raw_data.tasknum))>-1; % find practice trials and by extension, post-practice trials
real_task_idx = find(not_practice); real_task_idx = real_task_idx(1);
not_practice = false(height(raw_data),1); not_practice(real_task_idx:end) = true;
practice = ~not_practice; %a fun double negative
questionnaires = find(task_list==categorical({'NFC'}));
main_task = not_practice; 
if ~isempty(questionnaires)
    main_task(questionnaires(1):end) = false; %exclude questionnaire trials
else
    disp(['Subject ' num2str(raw_data.subjnum(1)) ' questionnaire data missing.'])
end

tasks = [categorical(cellstr('detection'));categorical(cellstr('combine'));categorical(cellstr('n-switch'));categorical(cellstr('n-back'))];
BDM_label = categorical(cellstr('BDM')); % find BDM trials in the sequence
BDMs = find(task_list == BDM_label);
BDM_rt = rts(BDMs);
displayed = raw_data.task_displayed(BDMs,:)';
displayed = double(displayed);

%first_trials = BDMs+4; 
%last_trials = BDMs-4;
%last_trials(1) = []; last_trials(end+1) = length(task_list) - 13; %pre-first BDM isn't a task trial
%debrief(debrief<last_trials(1)) = []; %only post-task debrief
debrief = find(task_list == categorical({'debrief'})&main_task);
task_trials = ismember(task_list,tasks)&main_task; %relevant trials. break into blocks by anything post-BDM and pre-debrief
for block = 1:default_length
    tasknum = block-1;
    relevant = raw_data(raw_data.tasknum == tasknum,:);
    first_trials(block) = find(raw_data.trial_index == relevant.trial_index(2));
    last_trials(block) = find(raw_data.trial_index == relevant.trial_index(end));
    TOT(block) = [relevant.time_elapsed(end)-relevant.time_elapsed(2)]/1000;
end
%first_trials = trials(1:15:end);
%last_trials = trials(15:15:end); %is this the solution? lol, no

% if sum(TOT>30|TOT<20)>0
%     idx = find(TOT>30|TOT<20)
%     TOT(idx)
%     subjnum
%     task_list(first_trials(idx))
%     task_list(last_trials(idx))
%     last_trials(idx)-first_trials(idx)
% end

if length(last_trials)<default_length %subject data got cut off for some unknown reason
    last_trials(end+1) = trials(end); %arbitrary last trial, since actual identity is unimportant
end

perf_list = table2array(raw_data(:,perfidx));
perf_by_block = perf_list(task_list==categorical({'debrief'})&main_task);
perf_by_block = str2num(char(perf_by_block));
if length(perf_by_block)<default_length
    disp(['check perf_by_block, subject ' string(subjnum) ' data may have stopped saving'])
    perf_by_block = [perf_by_block; 0];
end
if isempty(overall)
   overall = nanmean(perf_by_block);
end

%subject fair wages
value_list = double(string(table2array(raw_data(:,value))));
value_by_block = value_list(BDMs+1);
% BDM offers (random numbers)
offer_list = double(string(table2array(raw_data(:,offer))));
offer_by_block = offer_list(BDMs+1);

values = value_by_block; %points requested by block #
offers = offer_by_block;

%% get practice accuracy for an idea of 'baseline executive function'

max_practices = 14; 
acc = raw_data.practice_accuracy(~isnan(raw_data.practice_accuracy));
numpractices = max(double(string(raw_data.number_practice_hard)));
duplicates = [1; diff(acc(~isnan(acc)))]==0; %this isn't very precise but it's okay
%accuracy = acc(raw_data.task==categorical({'debrief'})); 
accuracy = acc(~isnan(acc)); accuracy(duplicates) = [];
practice_acc = NaN(1,max_practices); 
if numpractices > length(tasks) %only if they have to repeat the 4th task
    practice_acc(:,1:length(accuracy)) = accuracy;
else
    practice_acc(:,1:length(acc)) = acc;
end

% first block is detection, second block is 1-back, 3rd and onward is
% 2-back

%% compile rts by block and task
detectrts = NaN(1,default_length); nbackmeanrts = NaN(2,default_length); %defaults for concatenating failed subjects
for i=1:length(first_trials)
    if task_list(first_trials(i))==tasks(1)
        detectrts(:,i) = nanmean(rts(first_trials(i):last_trials(i)));
        nbackmeanrts(:,i) = NaN(2,1);
    elseif task_list(first_trials(i))==tasks(4)
        detectrts(:,i) = NaN;
        n = table2array(raw_data(first_trials(i),nidx));
        if n==1
            stacked = [nanmean(rts(first_trials(i):last_trials(i))); NaN];
        elseif n==2
            stacked = [NaN; nanmean(rts(first_trials(i):last_trials(i)))];
        end
        nbackmeanrts(:,i) = stacked;
    end
end

% starttime = raw_data(first_trials,timesteps); starttime = table2array(starttime);
% endtime = raw_data(last_trials,timesteps); endtime = table2array(endtime);
% TOT = (endtime-starttime)/1000;

task_progression = table2array(raw_data(first_trials,taskidx));

nbackacc = NaN(2,default_length); detectacc = NaN(1,default_length);
for i = 1:length(task_progression)
    task = task_progression(i);
    if task == tasks(1)
        detectacc(:,i) = perf_by_block(i);
        nbackacc(:,i) = NaN(2,1);
    elseif task == tasks(4)
        detectacc(:,i) = NaN;
        n = table2array(raw_data(first_trials(i),nidx));
        if n==1
            stacked = [perf_by_block(i); NaN];
            task_progression(i) = categorical({'n1'});
        elseif n==2
            stacked = [NaN; perf_by_block(i)];
            task_progression(i) = categorical({'n2'});
        end
        nbackacc(:,i) = stacked;
    end
end

%% look into characteristics of individual task completion (e.g. mean switch costs during a round of the combine)
init = NaN(1,default_length);
swrts = NaN(default_length,2); blrts = NaN(default_length,2);  %task switch rts, baseline rts for each task (detect, then nback)
%detectrts = [];
matchcount = NaN(1,default_length); missedmatches = NaN(1,default_length); nmatches = NaN(1,default_length); 
nswitches = NaN(1,default_length);
nswitchrts = NaN(2,default_length);
furthererrors = NaN(1,default_length);
% look into dual task cost in RT (so mean hit rts in combine - detect)
for block = 1:default_length
    %for trial = 1:length(first_trials)
    tasknum = block-1;
    relevant = raw_data(raw_data.tasknum==tasknum,:);
    %rts = table2array(raw_data(:,rtidx)); %initialize this properly
    %first = first_trials(trial);
    %last = last_trials(trial);
    %task_rts = rts(first:last);
    task_rts = relevant.rt;
    task_accuracy = relevant.new_correct;
    task = relevant.task(2);
    if task == tasks(1) %detection
        detect = relevant.detect == '1';
        rts = task_rts(detect);
    elseif task == tasks(2) %combine
        detect = table2array(raw_data(first:last,detectidx)) == '1'; %need to convert to logical
        nback = table2array(raw_data(first:last,nbackidx)) == '1'; %make sure that this is 1, or 'true' or something
        accurate = raw_data.new_correct(first:last);
        condition = init'; condition(nback) = 2; condition(detect) = 1; %always need to count detect trials as non-nbacks, no overlapping there
        rts = task_rts(~isnan(condition)&condition~=0); %rts(condition==0) = [];
        condition(isnan(condition)) = []; condition(condition==0) = []; %prune this, wasn't lining up before
        switches = [false; diff(condition)~=0];
        bl_detect = condition==1&~switches; bl_nback = condition==2&~switches; %bl short for "baseline"
        sw_detect = condition==1&switches; sw_nback = condition==2&switches;
        swrts(trial,:) = [nanmean(rts(sw_detect)) nanmean(rts(sw_nback))];
        blrts(trial,:) = [nanmean(rts(bl_detect)) nanmean(rts(bl_nback))];
        matchcount = [matchcount nansum(nback&accurate)];
    elseif task == tasks(3) %figure out switch costs in n-switch
        relevant = raw_data.task == tasks(3)&raw_data.new_correct&main_task; %only correct RTs, in main task, for n-switch task
        switches = table2array(raw_data(first:last,swidx))==categorical({'1'}) & relevant(first:last); %find switches, convert to logical
        notswitches = table2array(raw_data(first:last,swidx))==categorical({'0'}) & relevant(first:last); %find notswitches, convert to logical
        stacked = [nanmean(task_rts(notswitches)); nanmean(task_rts(switches))]; %let's not analyze first trials, since they're both not switches and not not switches
        nswitchrts(:,trial) = stacked;
        nswitches(trial) = sum(switches); %first trial is not a switch
        %[max(stacked) max(task_rts(switches)) min(task_rts(switches)) max(task_rts(notswitches)) min(task_rts(notswitches)) last-first]
        %print-out for debugging purposes, switch costs were all wrong here
    elseif task == tasks(4) %n-back (NOT combine)
        nback = (relevant.nback) == '1'; %make sure that this is 1, or 'true' or something
        accurate = relevant.new_correct;
        nmatches(block) = nansum(nback);
        missedmatches(block) = nansum(nback&~accurate);
        matchcount(block) = nansum(nback&accurate);
        if missedmatches(block)>0
            %if one error made, do you give up?
            matches = relevant(nback,:);
            errors = find(matches.new_correct==false);
            % errors made over matches still to go
            matchesleft = sum(find(matches.nback=='1')>errors(1));
            furthererrors(block) = sum(errors>errors(1))/(matchesleft);
            if matchesleft == 0
                furthererrors(block) = 0;
            end
            %so, error rate after 1 error was made
        end
    end
end

nbackmatches{1} = matchcount; %count number of matches, n

%% Calculate "effect" of n-back matches on BDM value (like set size effect)

%pick out individual tasks for match effect calculations
tasks = [categorical({'n1'}),categorical({'n2'})];
n1effect = []; n2effect = []; % keep track of numbers pulled out here for stats later
for task = 1:2 %cycle through 1- and 2-back
    idx = find(displayed == task); 
    completed = find(task_progression==tasks(task));
    for trial = 1:length(idx)
        now = idx(trial);
        if sum(completed<now)>0 %they've done the last once before
            last = completed(completed<now); last = last(end);
            delay = now-last;
            eval(['n' num2str(task) 'effect = [n' num2str(task) 'effect; matchcount(last) delay values(now) perf_by_block(now)];'])
        end
    end
end

range = 3:5;
n1matcheffect = NaN(1,length(range)); n2matcheffect = NaN(1,length(range));
for i = 1:length(range)
    match = range(i);
    idx = matchcount == match;
    idx(end) = []; idx = [false idx];
    idx2 = missedmatches == match-1;
    idx2(end) = []; idx2 = [false idx2];
    matcheffect(i) = nanmean(values(idx));
    missedeffect(i) = nanmean(values(idx2));
    existlogical = [~isempty(n1effect) ~isempty(n2effect)]; %sum(n1effect(:,1)==match)>0 sum(n2effect(:,1)==match)>0
    if existlogical(1)
        n1matcheffect(i) = nanmean(n1effect(n1effect(:,1)==match,3)); %pull out BDM in n1 when preceding matches was match
    end
    if existlogical(2)
        n2matcheffect(i) = nanmean(n2effect(n2effect(:,1)==match,3));
    end
end


%% look into frustration effects, changed answers, missed trials
% block by block
missed = NaN(1,default_length); lateresponse = missed; changedresponse = missed;
for block = 1:length(first_trials) % go by first-trial and last-trial indices
    first = first_trials(block); last = last_trials(block);
    ITIs = task_list(first:last) == categorical({'ITI'});
    fb = task_list(first:last) == categorical({'feedback'}); %it's 
    task = ~ITIs & ~fb;
    rts = raw_data.rt(first:last);
    missed(block) = sum(isnan(rts(task)));
    lateresponse(block) = sum(~isnan(rts(ITIs))) + sum(~isnan(rts(fb)));
    %changedresponse(block) = sum(~isnan(rts(task))&~isnan(rts(fb))); %responded both during task and during fb
end

%% make individual scores for BDM offer mimicry
y = values(2:end);
x = offers(1:end-1);
%[r,p] = corr(x,y);
%BDMmimicry = [r,p];

data.subjnum = subjnum;
data.version = exp_version;
data.fullscreen = fullscreen;
data.overall = overall;
data.practiceacc = practice_acc;
data.n1acc = nbackacc(1,:);
data.n2acc = nbackacc(2,:);
data.detectacc = detectacc;
data.values = values';
data.perf = perf_by_block';
data.task_progression = task_progression';
data.TOT = TOT;
data.n1rts = nbackmeanrts(1,:);
data.n2rts = nbackmeanrts(2,:);
data.detectrts = detectrts;
data.BDMrt = BDM_rt';
data.nbackmatches = nbackmatches;
data.nbackmisses = missedmatches;
data.intendednbackmatches = nmatches;
data.matcheffect = matcheffect;
data.n1matcheffect = n1matcheffect;
data.n2matcheffect = n2matcheffect;
data.missedeffect = missedeffect;
data.nbackfurthererrors = furthererrors;
data.offers = offers';
%data.BDMmimicry = BDMmimicry;
data.failed = failed;
data.missedtrials = missed;
data.lateresponse = lateresponse;
data.changedresponse = changedresponse;
data.task_displayed = displayed;
data.NFC = NFC; 
data.SAPS = SAPS;
data.SAPSs = SAPSs;
data.SAPSd = SAPSd;
%end