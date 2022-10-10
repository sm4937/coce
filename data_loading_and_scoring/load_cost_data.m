function [data,excluded] = load_cost_data(files)
%loadcostdata
%takes in relevant files, loads subject data in and prints it as data
%table, one row per subject

plot_flag = false;

if contains(files,'example_data')
    prefix = '';
else
    prefix = 'data/';
end

long_format = [];
for i = 1:length(files)
    file = files{i};
    raw =  load([prefix file]);
    colnames = raw.Untitled.Properties.VariableNames;
    raw.Untitled.TOT = double(string(raw.Untitled.TOT));
    raw.Untitled.total_points = double(string(raw.Untitled.total_points));
    raw.Untitled.stimnum = double(string(raw.Untitled.stimnum));
    raw.Untitled.prob_switch = double(string(raw.Untitled.prob_switch));
    raw.Untitled.rule = double(string(raw.Untitled.rule));
    raw.Untitled.practice_accuracy = double(string(raw.Untitled.practice_accuracy));
    raw.Untitled.overall = double(string(raw.Untitled.overall));
    raw.Untitled.success = double(string(raw.Untitled.success));
    raw.Untitled.number_practice_hard = double(string(raw.Untitled.number_practice_hard));
    raw.Untitled.response = double(string(raw.Untitled.response));
    raw.Untitled.tasknum = double(string(raw.Untitled.tasknum));
    raw.Untitled.practice = double(string(raw.Untitled.practice));
    raw.Untitled.n = double(string(raw.Untitled.n));
    possiblymissing = {'VarName43','VarName44','VarName45'}; %patch up the fact that JsPsych saves the data DIFFERENTLY EVERY FREAKING TIME!!!
    for name = 1:length(possiblymissing)
        if(~contains(colnames,possiblymissing{name}))
            eval(['raw.Untitled.' possiblymissing{name} ' = NaN(height(raw.Untitled),1);'])
        end
    end
    raw.Untitled.session = i*ones(height(raw.Untitled),1); %add in a session factor for random effects model
    %replace problem variables with doubles
    long_format = [long_format; raw.Untitled];
    %pull all the data into one long table
end
subjs = unique(long_format.subjnum(~isnan(long_format.subjnum)));
workers = unique(long_format.worker_ID(~isnan(long_format.subjnum)));
% the real worker IDs have been scrubbed from example_subjs.mat
% and replaced with random numbers (preserving that different workers are
% distinct from each other)
disp([num2str(length(workers)) ' unique workers who started in this batch'])

% process individual subjects, exclude subjects who are kicked out early,
% combine datasets into one big table
group = table; excluded = table;
for i = 1:length(subjs)
    subj = subjs(i);
    raw_data = long_format(long_format.subjnum == subj,:);
    if unique(raw_data.exp_version)==3
        [single,failed_counts] = make_data_table_v03(raw_data); %cleaned some stuff up for version 3, ease of analysis
    elseif unique(raw_data.exp_version)==4
        [single,failed_counts] = make_data_table_v04(raw_data); % % CURRENT TASK VERSION
    else
        [single,failed_counts] = make_data_table(raw_data);
    end
    if single.failed %separate people who got kicked out early
        excluded = [excluded; failed_counts];
    else
        group = [group; single];
    end
end

n = height(group);
data = group;

% BASIC DATA QUALITY CHECKS %%
%sanity check - BDM differences across sessions?
if plot_flag
    figure
    subplot(3,2,1)
    for i = 1:length(unique(data.session))
        y = data.values(data.session==i,:);
        y = reshape(y,numel(y),1);
        jitter = rand(length(y),1)-0.5;
        x = i*ones(length(y),1)+jitter;
        scatter(x,y,'o')
        hold on
    end
    title('BDM points requested by session number')
    ylabel('BDM points')
    xlabel('Session (jittered for visibility)')
    ax = gca; fig = gcf;
    fig.Color = 'w';
end

%any subjects with high amounts of BDM variability?
% i.e. unusual BDM behavior?
% short answer, no: the standard deviations of the BDMs we observe are
% fairly uniform, there's a very normal progression upward, no outliers to
% speak of

n1 = data.task_progression == categorical({'n1'});
% where are they doing the 1-back?
n2 = data.task_progression == categorical({'n2'});
% 2-back?

%check out the stds of BDMs on these specific tasks
for i = 1:height(data)
    dot = nanstd(data.values(i,n1(i,:)));
    dot2 = nanstd(data.values(i,n2(i,:)));
    list(i,:) = [dot dot2 i];
end

listn1 = sortrows(list,1,'Ascend');
listn2 = sortrows(list,2,'Ascend');

if plot_flag
    subplot(3,2,2)
    scatter(1:height(data),listn1(:,1))
    title('STD of BDM points requested per subject (1-back)')
    ylabel('STD')
    xlabel('Subject')
    xticks([1:height(data)])
    xticklabels(num2str(listn1(:,3)))
    ax = gca; fig = gcf;
    fig.Color = 'w';
    
    subplot(3,2,4)
    scatter(1:height(data),listn2(:,2))
    title('STD of BDM points requested per subject (2-back)')
    ylabel('STD')
    xlabel('Subject')
    xticks([1:height(data)])
    xticklabels(num2str(listn2(:,3)))
    fig = gcf;
    fig.Color = 'w';
    
    subplot(3,2,3)
    scatter(list(:,1),list(:,2),'Filled')
    title('relationship of BDM variances across tasks')
    ylabel('2-back STD BDMs')
    xlabel('1-back STD BDMs')
    
    subplot(3,2,5)
    missedBDMs = sum(isnan(data.values),2);
    matrix = [missedBDMs [1:length(missedBDMs)]'];
    matrix = sortrows(matrix,1,'ascend');
    plot(1:length(matrix),matrix(:,1),'o')
    ylabel('# Missed BDMS')
    title('Missed BDM answers by subject')
    xlabel('Subject number')
    xticks(1:n)
    xticklabels(num2str(matrix(:,2)))
    
end %plotflag

end

