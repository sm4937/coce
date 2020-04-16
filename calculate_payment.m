%% Pay subjects from MTurk
% go file by file
% not working yet

clear all
filename = '07.04.2020_trial';
subjs = [];
long_format = import_mturk_data(['./data/' filename '.csv']);
%raw_data = load('11.02.2020.mat','Untitled');
%long_format = raw_data.Untitled;
subjs = unique(string(long_format.worker_ID));
subjs(subjs=='0?assignmentId') = [];
subjs(subjs=='NULL') = []; %remove weird cases of data saving, debugging, etc

workers = [];total_points = [];TOTs = [];bonuses = [];codes =[];
for i = 1:length(subjs)
    subjnum = unique(long_format.subjnum(long_format.worker_ID==subjs(i)));
    if length(subjnum)>1
        repeatflag = true;
        subjnum = subjnum(2);
    else
        repeatflag = false;
    end
    indexs = find(long_format.subjnum == subjnum);
    index = indexs(end);
    if length(indexs)>5 %more than 5 trials completed
    subj_points = double(string(long_format.total_points(indexs)));
    points = nansum(subj_points);
    total_points = [total_points; points];
    TOT = [long_format.time_elapsed(indexs(end))-long_format.time_elapsed(indexs(3))]/60000;
    TOTs = [TOTs; TOT];
    worker =  subjs(i);
    worker = split(string(worker),'?');
    workers = [workers; worker(1)];
    codes = [codes; subjnum];
    bonus = (10*(TOT/60) + points/100) - 2.5;
    bonuses = [bonuses; bonus];
    disp(['Worker: ' worker(1) ' earned ' num2str(points) ' points and spent ' num2str(TOT) ' minutes on the task. Bonus worker $' num2str(bonus)])
    if repeatflag
        disp(['Worker ' worker(1) ' also restarted the task in the middle.'])
    end
    end %end more than 5 trials if statement
end

%to do: print this all to a CSV for record keeping
%here's a table for easy viewing

[workers(1:length(total_points)) codes(1:length(total_points)) TOTs(1:length(total_points)) total_points(1:length(total_points)) bonuses(1:length(total_points))]

%% why final screen bug on points screen?
broken_workers = {'A2GOYSTIL3LOV1','AA8PZKO9XGCKO','A12HZGOZQD5YK7','A1WZBK5PRHBW3H'};
for i = 1:length(broken_workers)
    indexes = find(contains(string(long_format.worker_ID),broken_workers{i}));
    if length(indexes)>100 %exclude subjects who did almost nothing or who are missing from the currently loaded file
    subj_points = double(string(long_format.points(indexes)));
    total = nansum(subj_points);
    TOT = [long_format.time_elapsed(indexes(end))-long_format.time_elapsed(indexes(3))]/60000;
    bonus = (10*(TOT/60) + total/100) - 2.5;
    bonuses = [bonuses; bonus];
    disp(['Worker: ' broken_workers{i} ' earned ' num2str(total) ' points and spent ' num2str(TOT) ' minutes on the task. Bonus worker $' num2str(bonus)])
    end
end

