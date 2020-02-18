function [data,failed_data] = make_data_table(raw_data)

%filename = ['\data\' filename];
data = table; failed_counts = [];
colnames = raw_data.Properties.VariableNames;
default_length = 24;
subjnum = unique(raw_data.subjnum);

fullscreenidx = 0;
for i = 1:length(colnames)
    if strmatch(colnames{i},'TOT')
        TOTidx = i;
    end
    if strmatch(colnames{i},'overall')
        overallidx = i;
    end
    if strmatch(colnames{i},'task')
        taskidx = i;
    end
    if strmatch(colnames{i},'time_elapsed')
        timesteps = i;
    end
    if strmatch(colnames{i},'trial_index')
        trialnum = i;
    end
    if strmatch(colnames{i},'rt')
        rtidx = i;
    end
    if strmatch(colnames{i},'success')
        fullscreenidx = i;
    end
    if strmatch(colnames{i},'detect')
        detectidx = i;
    end
end

nbackidx = find(ismember(colnames,{'nback'}));
perfidx = find(ismember(colnames,{'perf'}));
value = find(ismember(colnames,{'value'}));
offer = find(ismember(colnames,{'offer'}));
version_idx = ismember(colnames,{'exp_version'});
pswindex = find(ismember(colnames,{'prob_switch'}));
swidx = find(ismember(colnames,{'switch1'}));
%accidx = find(ismember(colnames,{'correct'}))
% classic accuracy column doesn't work here since a withheld response is
% correct so let's do
raw_data.new_correct = [categorical(raw_data.key_press) == categorical(raw_data.correct_key)]; %compare these hits/false alarms/misses
raw_data.new_correct(raw_data.correct_key==90) = true;
% so this is the whole thing, which is unfortunate because it includes
% feedback and instructions and BDM and stuff
% can be trimmed as you go
accidx = find(ismember(colnames,{'new_correct'})); 

%clean the switch column up, it's messy because of the zero's in it become
%"undefined"s
switch1 = table2array(raw_data(:,swidx));
valid_switch_data = (switch1==categorical({'1'})|switch1==categorical({'NULL'}));
raw_data.switch1(~valid_switch_data) = categorical({'0'});
%
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
    if tasknumlist(end) < 23
        data.failed = true;
        data.subjnum = subjnum; %print some stuff about people who fail, eventually
        failed_data = dropout(raw_data); %pass it into a deeper function to figure out they dropped out
        return
    else
        failed = false;
        failed_data = table;
        %demo = process_demo(raw_data); %we'll get this back
    end
else
    data.failed = true;
    data.subjnum = subjnum;
    failed_data = dropout(raw_data);
    return
end
% see whether subjects were in fullscreen mode the entire experiment
fullscreen = false;
if fullscreenidx ~= 0
    fullscreen = table2array(raw_data(end,fullscreenidx));
    if strmatch(fullscreen,"true")
        fullscreen = true;
    elseif fullscreen == categorical({'1'})
        fullscreen = true;
    else
        fullscreen = false;
    end
end
raw_data(end,:) = []; %trim out fullscreen trial with no data

% pull out need for cognition scores
%[NFC,NFCmeaned] = gradeNFC(raw_data);
%doesn't work for early pilots

rts = table2array(raw_data(:,rtidx));
logrts = log(rts);

overall = str2num(char(table2array(raw_data(end,overallidx)))); %overall accuracy
total_length = str2num(char(table2array(raw_data(end,TOTidx))))/60000; %convert msec to minutes
task_list = table2array(raw_data(:,taskidx));
exp_version = table2array(raw_data(1,version_idx));

tasks = [categorical(cellstr('detection'));categorical(cellstr('combine'));categorical(cellstr('n-switch'))];
BDM_label = categorical(cellstr('BDM')); % find BDM trials in the sequence
BDMs = find(task_list == BDM_label);
BDM_rt = rts(BDMs);

%data_long = raw_data(task_start:end,:);
debrief = find(task_list == categorical({'debrief'}));
first_trials = BDMs+4; %start of individual tasks
last_trials = BDMs-2;
last_trials(1) = []; last_trials(end+1) = length(task_list) - 9; %pre-first BDM isn't a task trial
%end of individual tasks

not_practice = double(string(raw_data.tasknum))>-1; % find practice trials and by extension, post-practice trials
real_task_idx = find(not_practice); real_task_idx = real_task_idx(1);
not_practice = false(height(raw_data),1); not_practice(real_task_idx:end) = true;
practice = ~not_practice; %a fun double negative
main_task = not_practice; main_task(last_trials(end)+2:end) = false; %exclude questionnaire trials

perf_list = table2array(raw_data(:,perfidx));
perf_by_block = perf_list(perf_list~=categorical({'NULL'})&main_task);
perf_by_block = str2num(char(perf_by_block));

%subject fair wages
value_list = table2array(raw_data(:,value));
value_by_block = value_list(value_list~=categorical({'NULL'})&main_task);
% BDM offers (random numbers)
offer_list = table2array(raw_data(:,offer));
offer_by_block = offer_list(offer_list~=categorical({'NULL'})&main_task);
for i=1:length(value_by_block)
    if value_list(i) == "" %subj didn't respond in time to BDM prompt
        value_list(i) = "0";
    end
end
values = str2num(char(value_by_block)); %points requested by block #
offers = str2num(char(offer_by_block));
% prune out weird values, gee thanks jspsych!
values = values(values<=5); offers = offers(offers<=5);

%% compile rts by block and task
detectrts = NaN(1,default_length); nbackrts = NaN(1,default_length); nswitchmeanrts = NaN(2,default_length); %defaults for concatenating failed subjects
for i=1:length(first_trials)
    if task_list(first_trials(i))==tasks(1)
        detectrts(:,i) = nanmean(logrts(first_trials(i):last_trials(i)));
        nbackrts(:,i) = NaN;
        nswitchmeanrts(:,i) = NaN(2,1);
    elseif task_list(first_trials(i))==tasks(2)
        nbackrts(:,i) = nanmean(logrts(first_trials(i):last_trials(i)));
        detectrts(:,i) = NaN;
        nswitchmeanrts(:,i) = NaN(2,1);
    elseif task_list(first_trials(i))==tasks(3)
        nbackrts(:,i) = NaN;
        detectrts(:,i) = NaN;
        pswitch = table2array(raw_data(first_trials(i),pswindex));
        if pswitch==0.1
            stacked = [nanmean(logrts(first_trials(i):last_trials(i))); NaN];
        elseif pswitch==0.9
            stacked = [NaN; nanmean(logrts(first_trials(i):last_trials(i)))];
        end
        nswitchmeanrts(:,i) = stacked;
    end
end

starttime = raw_data(first_trials,timesteps); starttime = table2array(starttime);
endtime = raw_data(last_trials,timesteps); endtime = table2array(endtime);
TOT = (endtime-starttime)/1000;

task_progression = table2array(raw_data(first_trials,taskidx));

nbackacc = NaN(1,default_length); detectacc = NaN(1,default_length); nswitchacc = NaN(2,default_length);
for i = 1:length(task_progression)
    task = task_progression(i);
    if task == tasks(2)
        nbackacc(:,i) = perf_by_block(i); %row of values, same length for both
        detectacc(:,i) = NaN;
        nswitchacc(:,i) = NaN(2,1);
    elseif task == tasks(1)
        detectacc(:,i) = perf_by_block(i);
        nbackacc(:,i) = NaN;
        nswitchacc(:,i) = NaN(2,1);
    elseif task == tasks(3)
        nbackacc(:,i) = NaN; %row of values, same length for both
        detectacc(:,i) = NaN;
        pswitch = table2array(raw_data(first_trials(i),pswindex));
        if pswitch==0.1
            stacked = [perf_by_block(i); NaN];
            task_progression(i) = categorical({'pswitch0.1'});
        elseif pswitch==0.9
            stacked = [NaN; perf_by_block(i)];
            task_progression(i) = categorical({'pswitch0.9'});
        end
        nswitchacc(:,i) = stacked;
    end
end

%% look into characteristics of individual task completion (e.g. mean switch costs during a round of the combine)
init = NaN(1,length(task_progression));
swrts = NaN(default_length,2); blrts = NaN(default_length,2); %task switch rts, baseline rts for each task (detect, then nback)
%detectrts = [];
matchcount = []; nswitches = NaN(1,default_length);
nswitchrts = NaN(2,default_length);
% look into dual task cost in RT (so mean hit rts in combine - detect)
for trial = 1:length(first_trials)
    rts = table2array(raw_data(:,rtidx)); %initialize this properly
    first = first_trials(trial);
    last = last_trials(trial);
    task_rts = table2array(raw_data(first:last,rtidx));
    task_accuracy = table2array(raw_data(first:last,accidx));
    task = task_list(first);
    if task == tasks(1) %detection
        detect = table2array(raw_data(first:last,detectidx)) == 'true';
        condition = init'; condition(detect) = 1;
        rts = task_rts(~isnan(condition)&condition~=0);
        %detectrts = [detectrts; nanmean(rts)];
    elseif task == tasks(2)
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
        relevant = raw_data.task == tasks(3);% & raw_data.practice==categorical({'false'}); %prune out feedback and ITI trials
        switches = table2array(raw_data(:,swidx))==categorical({'1'}) & relevant; %find switches, convert to logical
        notswitches = table2array(raw_data(:,swidx))==categorical({'0'}) & relevant; %find switches, convert to logical
        stacked = [nanmean(rts(notswitches(first+1:last))); nanmean(rts(switches(first+1:last)))]; %let's not analyze first trials, since they're both not switches and not not switches
        nswitchrts(:,trial) = stacked;
        nswitches(trial) = sum(raw_data.new_correct(first+1:last)&switches(first+1:last)); %first trial is not a switch
    end
end

correct_rts = []; 
switchcosts{1} = swrts-blrts;
combineswitchrts{1} = nanmean(swrts-blrts,2); %across tasks, switch costs within each combine block
nbackmatches{1} = matchcount;
allnswitchrts{1} = nswitchrts;
allnswitchcosts{1} = diff(nswitchrts);

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
data.nbackacc = nbackacc;
data.detectacc = detectacc;
data.nswitch9acc = nswitchacc(2,:);
data.nswitch1acc = nswitchacc(1,:);
data.values = values';
data.perf = perf_by_block';
data.task_progression = task_progression';
data.TOT = TOT';
data.nbackrts = nbackrts;
data.detectrts = detectrts;
data.nswitch9rts = nswitchmeanrts(2,:);
data.nswitch1rts = nswitchmeanrts(1,:);
data.nswitchblrts = nswitchrts(1,:);
data.nswitchswrts = nswitchrts(2,:);
data.BDMrt = BDM_rt';
data.combineswitchcosts = switchcosts;
data.combineswitchrts = combineswitchrts;
data.nbackmatches = nbackmatches;
data.nswitches = nswitches;
data.offers = offers';
%data.BDMmimicry = BDMmimicry;
data.failed = failed;
data.allnswitchrts = allnswitchrts;
data.allnswitchcosts = allnswitchcosts;
data.missedtrials = missed;
data.lateresponse = lateresponse;
data.changedresponse = changedresponse;
%data.NFC = NFC;
% data.age = demo.age;
% data.sex = demo.sex;
% data.edu = demo.edu;

%end