function [data] = make_data_table(subjnum)
%subjnum = 203;
filename = ['.\pilotdata\Pilot' num2str(subjnum) '.txt'];
raw_data = import_jspsych_txt(filename);
if subjnum == 4
    raw_data = import_subj4_txt(filename);
elseif subjnum == 5
    raw_data = import_subj5_txt(filename);
elseif subjnum == 6
    raw_data = import_subj6_txt(filename);
elseif subjnum >= 7
    raw_data = load(['.\pilotdata\Pilot' num2str(subjnum) '.mat']);
    eval(['raw_data = raw_data.Pilot' num2str(subjnum) ';'])
end

colnames = raw_data.Properties.VariableNames;
default_length = 12;

to_delete = {'stimulus','internal_node_id'}; %clean things up for your sanity scrolling through all this shit
for var = 1:length(to_delete)
    column = to_delete(var);
    idx = find(ismember(colnames,column));
    raw_data(:,idx) = [];
    colnames(:,idx) = [];
end
%raw_data(:,find(ismember(colnames,{'stimulus'}))) = []; %remove 'Stimulus' column since it's bulky and not informative
%raw_data(:,find(ismember(colnames,{'internal_node_id'}))) = []; %remove 'internal_node_id' column since it's also bulky

relevant_vars = {'TOT','overall','performance_by_block','task','time_elapsed','trial_index','value_list','rt','success',...
'detect','nback','correct'}; %easier way to do the following colname match operation?
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
    if strmatch(colnames{i},'internal_node_id')
        internalnode = i;
    end
    if strmatch(colnames{i},'stimulus')
        stimulus = i;
    end
    if strmatch(colnames{i},'detect')
        detectidx = i;
    end
end

nbackidx = find(ismember(colnames,{'nback'}));
perfidx = find(ismember(colnames,{'performance_by_block'}));
value = find(ismember(colnames,{'value_list'}));
offer = find(ismember(colnames,{'offer_list'}));
version = find(ismember(colnames,{'exp_version'}));
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

%separate out the people who failed from the people who didn't 
tasknumlist = raw_data.tasknum(~isnan(raw_data.tasknum));
if tasknumlist(end) < 11
    failed = true;
else
    failed = false;
end
% see whether subjects were in fullscreen mode the entire experiment
fullscreen = false;
if fullscreenidx ~= 0
    fullscreen = table2array(raw_data(end,fullscreenidx));
    if strmatch(fullscreen,"true")
        fullscreen = true;
    else
        fullscreen = false;
    end
end
if subjnum >= 5 % to be removed when no longer analyzing pilot data
    raw_data(end,:) = []; %trim out fullscreen trial with no data
end

perf = split(table2array(raw_data(end,perfidx)),',');
perf_by_block = str2num(char(perf)); %performance by block #
value_list = split(table2array(raw_data(end,value)),',');
offer_list = split(table2array(raw_data(end,offer)),',');
for i=1:length(value_list)
    if value_list(i) == "" %subj didn't respond in time to BDM prompt
        value_list(i) = "0";
    end
end
values = str2num(char(value_list)); %points requested by block #
offers = str2num(char(offer_list));
if subjnum == 2 % you accidentally deleted this so put it here
    value_list = ["4.04", "3.08", "4.04", "4.04", "1.96", "1.80", "1.80", "1.48", "1.64", "1.16", "2.12", "1.64"];
    for i = 1:length(value_list)
        values(i) = str2double(char(value_list(i)));
    end
end

rts = table2array(raw_data(:,rtidx));

overall = str2num(char(table2array(raw_data(end,overallidx)))); %overall accuracy
total_length = str2num(char(table2array(raw_data(end,TOTidx))))/60000; %convert msec to minutes
task_list = table2array(raw_data(:,taskidx));
if ~isempty(version)
    exp_version = str2num(char(table2array(raw_data(end,version))));
else
    exp_version = 1;
end

tasks = [categorical(cellstr('detection'));categorical(cellstr('combine'));categorical(cellstr('n-switch'))];
BDM_label = categorical(cellstr('BDM')); % find BDM trials in the sequence
BDMs = find(task_list == BDM_label);
BDM_rt = rts(BDMs);

%data_long = raw_data(task_start:end,:);
first_trials = BDMs+4; %start of individual tasks
last_trials = BDMs-2; %pre-first BDM isn't a task trial
last_trials(1) = []; last_trials(end+1) = height(raw_data)-9;
%end of individual tasks
detectrts = NaN(1,default_length); nbackrts = NaN(1,default_length); nswitchmeanrts = NaN(2,default_length); %defaults for concatenating failed subjects
for i=1:length(first_trials)
    if task_list(first_trials(i))==tasks(1)
        detectrts(:,i) = nanmean(raw_data.rt(first_trials(i):last_trials(i)));
        nbackrts(:,i) = NaN;
        nswitchmeanrts(:,i) = NaN(2,1);
    elseif task_list(first_trials(i))==tasks(2)
        nbackrts(:,i) = nanmean(raw_data.rt(first_trials(i):last_trials(i)));
        detectrts(:,i) = NaN;
        nswitchmeanrts(:,i) = NaN(2,1);
    elseif task_list(first_trials(i))==tasks(3)
        nbackrts(:,i) = NaN;
        detectrts(:,i) = NaN;
        pswitch = table2array(raw_data(first_trials(i),pswindex));
        if pswitch==0.1
            stacked = [nanmean(rts(first_trials(i):last_trials(i))); NaN];
        elseif pswitch==0.9
            stacked = [NaN; nanmean(rts(first_trials(i):last_trials(i)))];
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
matchcount = [];
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
        detect = table2array(raw_data(first:last,detectidx)) == 'true'; %need to convert to logical
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
        relevant = raw_data.task == tasks(3); %prune out feedback and ITI trials
        switches = table2array(raw_data(:,swidx))==categorical({'true'}); %find switches, convert to logical
        stacked = [nanmean(rts(relevant(first:last)&~switches(first:last))); nanmean(rts(relevant(first:last)&switches(first:last)))];
        nswitchrts(:,trial) = stacked;
    end
end

correct_rts = []; 
switchcosts{1} = swrts-blrts;
combineswitchrts{1} = nanmean(swrts-blrts,2); %across tasks, switch costs within each combine block
nbackmatches{1} = matchcount;
allnswitchrts{1} = nswitchrts;
allnswitchcosts{1} = diff(nswitchrts);

%% make individual scores for BDM offer mimicry
y = values(2:end);
x = offers(1:end-1);
[r,p] = corr(x,y);
BDMmimicry = [r,p];

data = table;
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
data.offers = offers';
data.BDMmimicry = BDMmimicry;
data.failed = failed;
data.allnswitchrts = allnswitchrts;
data.allnswitchcosts = allnswitchcosts;

%end