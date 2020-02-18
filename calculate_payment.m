%% Pay subjects from MTurk
% go file by file
% not working yet

filename = '13.02.2020_trial';
subjs = [];
long_format = import_mturk_data(['./data/' filename '.csv']);
%raw_data = load('11.02.2020.mat','Untitled');
%long_format = raw_data.Untitled;
subjs = unique(long_format.worker_ID);
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
    subj_points = double(string(long_format.total_points(indexs)));
    points = unique(subj_points(~isnan(subj_points)));
    if sum(~isnan(subj_points))<1 %happens sometimes when pulling from trial and not from test
        total_points = [total_points; 0];
    else
        total_points = [total_points; points(end)];
    end
    TOT = [long_format.time_elapsed(indexs(end))-long_format.time_elapsed(indexs(3))]/60000;
    TOTs = [TOTs; TOT];
    worker =  subjs(i);
    worker = split(string(worker),'?');
    workers = [workers; worker(1)];
    codes = [codes; subjnum];
    bonus = (10*(TOT/60) + total_points(end)/100) - 2;
    bonuses = [bonuses; bonus];
    disp(['Worker: ' + worker(1) + ' earned ' + points + ' points and spent ' + TOT + ' minutes on the task. Bonus worker $' + bonus])
    if repeatflag
        disp(['Worker ' + worker(1) + ' also restarted the task in the middle.'])
    end
end

%to do: print this all to a CSV for record keeping
%here's a table for easy viewing

[workers(1:length(total_points)) codes(1:length(total_points)) TOTs(1:length(total_points)) total_points(1:length(total_points)) bonuses(1:length(total_points))]