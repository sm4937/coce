%% Load in subject data for use in simulations
% this code stands alone
% it loads subject data all on its own, formats it, runs it through our
% simple WM/attention process model, and outputs the operations that the
% model uses to complete each task.
% these cost components are then fed into the modelling code in /modelling,
% through loading toanalyze.mat, a file that gets saved at the end of the
% second section of this script.

% in addition, there's some analysis in here of the overall distributions
% of the cost components being extracted in this script (seeing whether
% they are correlated, normally distributed, etc., etc.)
% this level of analysis is more important on the computational modeling
% side of things, and doesn't aim to understand subject behavior, hence it
% stands alone here & is not tied to behavioral analysis in
% paper_graphs_and_stats.m or run_supplementary_analyses.m.

clear all

if isdir('../data')
    % ALL experimental data, all 100 subjects live in 'data'
    load('../data/filenames.mat')
    load('../data/fullsubjnumbers.mat')
    data_directory = '../data/';
    % grab all subjects from those files, then
elseif isdir('../example_data')
    files{1} = 'example_subjs.mat';
    load('../example_data/fullsubjnumbers.mat')
    data_directory = '../example_data/';
end
% specify behavioral files collected on version 4 of WFW behavioral
% task


    
load([data_directory 'toanalyze.mat'],'toanalyze','trim')
    

%% Check out individual differences in measures,
% relationships between measures

subjnums = unique(toanalyze.subj);
nsubjs = length(subjnums);

measures(:,1) = toanalyze.nmisses; measures(:,2) = toanalyze.nlures;
measures(:,3) = toanalyze.nFAs; measures(:,4) = toanalyze.maintained;
names = {'misses','lures','FAs','maintenance'};
for col = 1:size(measures,2)
    figure; subplot(2,2,2)
    measure = measures(:,col);
    for misses = 0:4
        for s = 1:nsubjs
            subj = subjnums(s);
            nmisses(s) = max(measure(toanalyze.subj==subj));
            BDMs = toanalyze.BDM(toanalyze.subj==subj,:);
            nmissesall(:,s) = [measure(toanalyze.subj==subj); nan(32-sum(toanalyze.subj==subj),1)];
            misseffect(misses+1,:,s) = NaN(1,32); %initialize empty
            idx = find(nmissesall(:,s)==misses)+1; idx(idx>length(BDMs)) = []; %trim out edges
            misseffect(misses+1,idx,s) = BDMs(idx); %catalog BDM effects from misses
        end
        scatter(1:nsubjs,sum(nmissesall==misses),'Filled')
        hold on
    end
    legend({'0','1','2','3','4'})
    title(['All ' names{col} ' per Subject'])
    ylabel('Frequency'); xlabel('Subject')
    subplot(2,2,3)
    histogram(nmisses)
    title(['Max number of ' names{col} ' for each subject'])
    subplot(2,2,1)
    histogram(measure)
    title(['All ' names{col} ' values all subjects'])
    fig = gcf; fig.Color = 'w';
    
    subplot(2,2,4)
    for misses = 0:4
        scatter(ones(nsubjs,1)*misses,nanmean(misseffect(misses+1,:,:),2),'Filled')
        hold on
    end
    for subj = 1:nsubjs
        plot(0:4,nanmean(misseffect(:,:,subj),2),'k--')
        hold on
    end
    fig = gcf; fig.Color = 'w';
    title(['Mean BDM following each # of ' names{col}])
    
end

[r,p] = corr(measures,'type','Spearman'); % not normally distributed
% are these cost components correlated?
% across all subjects, yes, which is unsurprising overall. maintenance,
% lures, & errors all go up with 2-back.

for n = 1:nsubjs
    % are they related within-subject?
    onesubj = toanalyze(toanalyze.subj==n,:);
    
    [r,p] = corr(onesubj.maintained,onesubj.nlures,'type','Spearman');
    proportions(n,1) = p<0.05;
    [r,p] = corr(onesubj.maintained,onesubj.nFAs,'type','Spearman');
    proportions(n,2) = p<0.05;
    [r,p] = corr(onesubj.nlures,onesubj.nFAs,'type','Spearman');
    proportions(n,3) = p<0.05;
    % get percentage of subjects w/ significant correlation of all these
    % components
end

percentages = sum(proportions)/nsubjs;

%% understand effect of delay between task iterations on BDM values - is there anything there ?

BDMs = [0:4:100];
BDMs = 1+(BDMs./25);
nblocks = 32;
sim_list = [];
for p = 1:length(BDMs) %cycle through fixed values of BDMs - mirror consistency of our subjects
    value = BDMs(p);
    tasks = [ones(1,ceil(nblocks/3)) 2.*ones(1,ceil(nblocks/3)) 3.*ones(1,floor(nblocks/3))];
    tasks = tasks(randperm(nblocks));
    for block = 1:32
        offer = BDMs(ceil(rand*length(BDMs)));
        if value>offer
            tasks(block) = 0;
        end
    end
    sim_list = [sim_list; tasks value];
end

%so, given 1 fixed fair wage for all tasks, how does completion of the
%default task (0) change? how does completion of the other tasks change?
%does it look like our subjects have fixed BDM values?
figure
for subj = 1:length(BDMs)
    value = BDMs(subj);
    row = sim_list(:,end)==value;
    subplot(5,6,subj)
    histogram(sim_list(row,1:end-1))
    xticklabels({'0','1','2','3'})
    xlabel('Task')
    ylabel('Freq of completion given fixed BDM value for every task')
    title(['value = ' num2str(value)])
end
fig = gcf; fig.Color = 'w';

delay1 = []; delay2 = [];
for task = 1:2
    for subj = 1:size(sim_list,1)
        idx = find(sim_list(subj,1:end-1)==task);
        value = sim_list(subj,end);
        for block = 2:length(idx)
            now = idx(block);
            last = idx(block-1);
            eval(['delay' num2str(task) ' = [delay' num2str(task) '; now-last value];'])
        end
    end
end

figure
scatter(delay1(:,1),delay1(:,2),'b','Filled')
hold on
scatter(delay2(:,1),delay2(:,2),'r','Filled')
title('Delay btwn tasks by fixed BDM value (sim.)')
ylabel('Fixed BDM')
xlabel('Delay between task iterations')
legend({'1-back','2-back'})
fig = gcf; fig.Color = 'w';




