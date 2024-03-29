function [data,failed_data] = make_data_table_v04(raw_data)

% when data comes in from the experiment, it's in a very long JSON file
% which includes data on every single event in the experiment, including
% stim presentation, feedback presentation, BDM ratings, etc.

% data is loaded up using data_loading_and_scoring/load_cost_data.m, then
% is fed into this function in the raw_data table. raw_data contains all
% trials' information for ONE subject at a time. if the subject passes
% exclusion criteria (basically, if they finished the experiment itself),
% then their data are formatted into a table called "data." if not,

% this function serves to make these raw events (raw_data), which have been
% preprocessed out of JSON format into MATLAB table format, into a usable
% data table with all relevant columns (RTs, BDM fair wage ratings (value),
% computer offers in the auction (offer), whether a trial is an n-back or
% detection match, subject identification number, which task was rated
% and which task was completed).

% as such, this function includes a lot of trial trimming, lots of datatype
% conversion, etc.

raw_data = sortrows(raw_data,'trial_index'); %re-sort scrambled data
% sometimes it comes in out of order, but the trial_index restores it

data = table;
% initialize clean MATLAB table

colnames = raw_data.Properties.VariableNames;
% what do we have saved in raw_data?

default_length = 32;
subjnum = unique(raw_data.subjnum);
NFC = NaN; SAPS = NaN;

perfidx = find(ismember(colnames,{'perf'}));

% classic accuracy column doesn't work here since a withheld response is
% correct so I re-score below:
raw_data.new_correct = categorical(raw_data.key_press) == categorical(raw_data.correct_key); %compare these hits/false alarms/misses
raw_data.new_correct(raw_data.correct_key==90) = true;

% so this is the entire experiment, including junk throw-away trials
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
tasknumlist = unique(tasknumlist(~isnan(tasknumlist)&tasknumlist>-1));
if sum(~isnan(tasknumlist))>0
    if length(tasknumlist) < default_length
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
% we started out with some fullscreen versions, but because implementation
% is not uniform across browsers, etc., we assume everyone's data come from
% fullscreen 
exp_version = unique(raw_data.exp_version);

% pull out need for cognition scores
[NFC,SAPS,SAPSs,SAPSd] = grade_scales(raw_data);
% pull out difficulty ratings
[task_ratings] = grade_difficulty_ratings(raw_data);

%rts = table2array(raw_data(:,ismember(colnames,{'rt'})));
rts = raw_data.rt;
overall_list = double(string(raw_data.overall)); %overall accuracy
overall = unique(overall_list(~isnan(overall_list)));

total_length = raw_data.TOT(end)/60000; %convert msec to minutes
% this is how long people spent on the entire experiment, including
% reading/filling out consent and data consent forms
totalTOT = [raw_data.time_elapsed(end)-raw_data.time_elapsed(3)]/60000;
% this is how many minutes subjects spent on the experiment POST-consent
% forms. given that workers often accept AMT HITs, open them up, and then
% sit with then open for a while actually starting the experiment, I think
% this is a more accurate measure of t

task_list = raw_data.task;
tasks = [categorical(cellstr('detection'));categorical(cellstr('combine'));categorical(cellstr('n-switch'));categorical(cellstr('n-back')); categorical(cellstr('ndetection'))];
BDM_label = categorical(cellstr('BDM')); % find BDM trials in the sequence

not_practice = double(string(raw_data.tasknum))>-1;
% find practice trials and by extension, post-practice trials
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

BDMs = find(task_list == BDM_label);
if length(BDMs)>default_length %more than 32 BDMs saved = there's a duplicate trial saved
    match = find(diff(raw_data.tasknum(BDMs))==0);
    BDMs(match) = []; %delete tag to duplicate trial
end
BDM_rt = NaN(1,default_length);
BDM_rt(1:length(BDMs)) = rts(BDMs);
displayed = raw_data.task_displayed(BDMs,:)';
d = default_length-length(displayed);
displayed = [double(string(displayed)) NaN(1,d)];

debrief = find(task_list == categorical({'debrief'})&main_task);
task_trials = ismember(task_list,tasks)&main_task; %relevant trials. break into blocks by anything post-BDM and pre-debrief
for block = 1:default_length
    tasknum = block-1;
    relevant = raw_data(raw_data.tasknum == tasknum,:);
    if height(relevant)>3
        try
            first_trials(block) = find(raw_data.trial_index == relevant.trial_index(2));
        catch
            random_duplicates = find(raw_data.trial_index == relevant.trial_index(2));
            first_trials(block) = random_duplicates(2);
        end
        last_trials(block) = find(raw_data.trial_index == relevant.trial_index(end));
        TOT(block) = [relevant.time_elapsed(end)-relevant.time_elapsed(2)]/1000;
    else
        first_trials(block) = NaN;
        last_trials(block) = NaN;
        TOT(block) = NaN;
    end
end
last_trials(isnan(first_trials)) = [];
first_trials(isnan(first_trials)) = [];

perf_list = raw_data.perf;
perf_by_block = perf_list(task_list==categorical({'debrief'})&main_task);
perf_by_block = str2num(char(perf_by_block));
if length(perf_by_block)<default_length
    disp(['check perf_by_block, subject ' string(subjnum) ' data may have stopped saving'])
    d = default_length-length(perf_by_block);
    perf_by_block = [perf_by_block; NaN(d,1)];
elseif length(perf_by_block)>default_length
    perf_by_block(1) = []; %practice block snuck in, seen this once when length(perf) = 33
end
if isempty(overall)
    overall = nanmean(perf_by_block);
end

%subject fair wages
value_list = double(string(table2array(raw_data(:,ismember(colnames,{'value'})))));
value_by_block = value_list(BDMs+1);
% BDM offers (random numbers)
offer_list = double(string(table2array(raw_data(:,ismember(colnames,{'offer'})))));
offer_by_block = offer_list(BDMs+1);
values = NaN(1,default_length);
values(1:length(value_by_block)) = value_by_block; %points requested by block
offers = NaN(1,default_length); %offers by block
offers(1:length(offer_by_block)) = offer_by_block; %determines actual points on offer

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

%% compile rts and accuracies by block and task
task_progression = repmat(categorical(cellstr('NULL')),1,default_length);
nbackacc = NaN(2,default_length); detectacc = NaN(1,default_length); ndetectacc = NaN(1,default_length);
detectrts = NaN(2,default_length); n1rts = NaN(2,default_length); n2rts = NaN(2,default_length); ndetectrts = NaN(2,default_length); %defaults for concatenating failed subjects
meanRTs = NaN(1,default_length);

for block = 1:length(first_trials)
    tasknum = block-1;
    relevant = raw_data(raw_data.tasknum==tasknum,:);
    relevant = relevant(2:end,:); %trim first row, which is the BDM trial preceding the block starts
    correct = relevant.rt(relevant.new_correct);
    incorrect = relevant.rt(~relevant.new_correct);
    currenttask = relevant.task(1);
    task_progression(block) = currenttask;
    switch currenttask
        case tasks(1)
            detectrts(:,block) = [nanmean(correct); nanmean(incorrect)];
            detectacc(:,block) = perf_by_block(block);
        case tasks(4)
            n = relevant.n(1);
            eval(['n' num2str(n) 'rts(:,block) = [nanmean(correct); nanmean(incorrect)];'])
            nbackacc(n,block) = perf_by_block(block);
            task_progression(block) = categorical(cellstr(['n' num2str(n)]));
        case tasks(5)
            ndetectrts(:,block) = [nanmean(correct); nanmean(incorrect)];
            ndetectacc(:,block) = perf_by_block(block);
    end
    meanRTs(block) = nanmean(correct);
end
%% look into characteristics of individual task completion (e.g. ER post-first error)
init = NaN(1,default_length);
swrts = NaN(default_length,2); blrts = NaN(default_length,2);  %task switch rts, baseline rts for each task (detect, then nback)
%detectrts = [];
matchcount = NaN(1,default_length); missedmatches = NaN(1,default_length); nmatches = NaN(1,default_length); distractors = zeros(1,default_length);
nswitches = NaN(1,default_length); ndetectcount = NaN(1,default_length); ndetectmisses = NaN(1,default_length);
nswitchrts = NaN(2,default_length);
furthererrors = NaN(1,default_length); posterrorRT = NaN(1,default_length);
% look into dual task cost in RT (so mean hit rts in combine - detect)
for block = 1:length(first_trials)
    %for trial = 1:length(first_trials)
    tasknum = block-1;
    relevant = raw_data(raw_data.tasknum==tasknum,:);
    task_rts = relevant.rt;
    task_accuracy = relevant.new_correct;
    task = relevant.task(2);
    if task == tasks(1) %detection
        detect = relevant.detect == '1';
        rts = task_rts(detect);
    elseif task == tasks(4) %n-back (NOT combine)
        nback = (relevant.nback) == '1'; %make sure that this is 1, or 'true' or something
        n = relevant.n(2); idxes = find(nback);
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
            elseif matchesleft > 0
                remaining = find(matches.nback=='1')>errors(1);
                posterrorRT(block) = nanmean(matches.rt(remaining));
            end
            %so, error rate after 1 error was made
        end
        for i = 1:nmatches(block)
            idx = idxes(i);
            try
                within = relevant(idx-n:idx,:);
                other = nansum(within.stimnum(1:n-1)~=within.stimnum(n)); %count stimuli that are not the same
                % for 3detect this should always be 0
                distractors(block) = distractors(block) + other;
            catch
                distractors(block) = distractors(block) + 0; % do nothing
            end
        end
    elseif task == tasks(5) %ndetect (2-back but continuous evidence accumulation)
        detect = relevant.detect == '1';
        ndetectcount(block) = nansum(detect&relevant.new_correct);
        ndetectmisses(block) = nansum(detect&~relevant.new_correct);
        if ndetectmisses(block)>0
            %if one error made, do you give up?
            matches = relevant(detect,:);
            errors = find(matches.new_correct==false);
            % errors made over matches still to go
            matchesleft = sum(find(matches.detect=='1')>errors(1));
            furthererrors(block) = sum(errors>errors(1))/(matchesleft);
            if matchesleft == 0
                furthererrors(block) = 0;
            elseif matchesleft > 0
                remaining = find(matches.detect=='1')>errors(1);
                posterrorRT(block) = nanmean(matches.rt(remaining));
            end
            % same variable for every task, can separate by block
            % #/task_progression
        end
    end
end

nbackmatches = matchcount; %count number of matches, n

%% Calculate "effect" of n-back matches on BDM value (like set size effect)

% %pick out individual tasks for match effect calculations
% tasks = [categorical({'n1'}),categorical({'n2'})];
% n1effect = []; n2effect = []; % keep track of numbers pulled out here for stats later
% for task = 1:2 %cycle through 1- and 2-back
%     idx = find(displayed == task);
%     completed = find(task_progression==tasks(task));
%     for trial = 1:length(idx)
%         now = idx(trial);
%         if sum(completed<now)>0 %they've done the last once before
%             last = completed(completed<now); last = last(end);
%             delay = now-last;
%             eval(['n' num2str(task) 'effect = [n' num2str(task) 'effect; matchcount(last) delay values(now) perf_by_block(now)];'])
%         end
%     end
% end
%
% range = 2:4;
% n1matcheffect = NaN(1,length(range)); n2matcheffect = NaN(1,length(range));
% for i = 1:length(range)
%     match = range(i);
%     idx = matchcount == match;
%     idx(end) = []; idx = [false idx];
%     idx2 = missedmatches == match-1;
%     idx2(end) = []; idx2 = [false idx2];
%     matcheffect(i) = nanmean(values(idx));
%     missedeffect(i) = nanmean(values(idx2));
%     existlogical = [~isempty(n1effect) ~isempty(n2effect)]; %sum(n1effect(:,1)==match)>0 sum(n2effect(:,1)==match)>0
%     if existlogical(1)
%         n1matcheffect(i) = nanmean(n1effect(n1effect(:,1)==match,3)); %pull out BDM in n1 when preceding matches was match
%     end
%     if existlogical(2)
%         n2matcheffect(i) = nanmean(n2effect(n2effect(:,1)==match,3));
%     end
% end


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

%% why does detection task take longer in some people than in others?

raw_data.time_elapsed(1);
detect = find(task_progression==categorical(cellstr('detection')));
phases = [categorical(cellstr('detection')) categorical(cellstr('feedback')) categorical(cellstr('ITI'))];
for block = 1:length(detect)
    tasknum = detect(block)-1;
    trials = raw_data.trial_index(raw_data.tasknum==tasknum);
    first = find(raw_data.trial_index==trials(1)); last = find(raw_data.trial_index==(trials(end)+2));
    relevant = raw_data(first:last,:);
    realfirst = find(relevant.task==categorical(cellstr('instructions')));
    relevant = relevant((realfirst+1):end,:);
    trialtimes = [NaN; diff(relevant.time_elapsed)];
    numresp = nansum(~isnan(relevant.rt));
    for phase = 1:length(phases)
        totaltimes(phase) = nansum(trialtimes(relevant.task==phases(phase)));
    end
    bookkeeping(block,:) = [subjnum totaltimes TOT(block) numresp];
end

detect_time_breakdown{1} = bookkeeping;

%% Individual relationships of BDMs with other values
% This measure, which is explored in run_supplementary_analyses, is meant
% to understand whether subjects were rating all tasks in the same way, or
% perhaps whether subjects are rating the 1-back and 3-detect the same way
% while the 2-back task is more of a "popout" costly task.

y = values(2:end);
x = offers(1:end-1);
%[r,p] = corr(x,y);
%BDMmimicry = [r,p];
tasknumbers = [0 1 2 7];
BDMautocorrs = NaN(1,length(tasknumbers)-1);
BDMacrosscorrs = NaN(3,2);
for task = 1:(length(tasknumbers)-1)
    onscreen = find(displayed==tasknumbers(task+1)); lag = onscreen(2:end);
    others = [1 2 7];
    others(others==task) = [];
    for second = 1:(length(tasknumbers)-2)
        othertask = find(displayed==others(second));
        columns{1} = othertask; columns{2} = onscreen;
        [~,longer] = max([length(othertask),length(onscreen)]);
        d = abs(length(othertask)-length(onscreen)); d = d-1;
        columns{longer}(end-d:end) = [];
        [r,p] = corr(values(columns{1})',values(columns{2})');
        BDMacrosscorrs(task,second) = r; %original task each row, second task is column (so 3x2)
    end
    onscreen(end) = [];
    [r,p] = corr(values(onscreen)',values(lag)'); %task autocorrelations
    BDMautocorrs(task) = r;
end
taskBDMcorrs{1} = BDMacrosscorrs;
taskBDMcorrs{2} = BDMautocorrs;

%% save it all in a compressed MATLAB table, called data
% this is a much cleaner data format
% some things have been already calculated here, too, like meanRTs, mean
% accuracy on each type of task, etc.

data.subjnum = subjnum;
data.version = exp_version;
data.fullscreen = fullscreen;
data.overall = overall;
data.totalTOT = totalTOT;
data.practiceacc = practice_acc;
data.n1acc = nbackacc(1,:);
data.n2acc = nbackacc(2,:);
data.detectacc = detectacc;
data.ndetectacc = ndetectacc;
data.values = values;
data.perf = perf_by_block';
data.meanRTs = meanRTs;
data.task_progression = task_progression;
data.TOT = TOT;
data.n1rts = n1rts(1,:);
data.n2rts = n2rts(1,:);
data.detectrts = detectrts(1,:);
data.ndetectrts = ndetectrts(1,:);
data.incorrectrts{1} = [detectrts(2,:); n1rts(2,:); n2rts(2,:); ndetectrts(2,:)];
data.BDMrt = BDM_rt;
data.nbackmatches = nbackmatches;
data.nbackmisses = missedmatches;
data.intendednbackmatches = nmatches;
% data.matcheffect = matcheffect;
% data.n1matcheffect = n1matcheffect;
% data.n2matcheffect = n2matcheffect;
% data.missedeffect = missedeffect;
data.posterrorER = furthererrors;
data.posterrorRT = posterrorRT;
data.ndetectmatches = ndetectcount;
data.ndetectmisses = ndetectmisses;
data.distractors = distractors;
data.offers = offers;
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
data.taskratings = task_ratings;
data.detect_time_breakdown = detect_time_breakdown;
data.taskBDMcorrs = taskBDMcorrs;

end % of function