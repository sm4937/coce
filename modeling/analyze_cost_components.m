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

addpath(genpath('../plotting/'))

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
measures(:,5) = toanalyze.nupdates; measures(:,6) = toanalyze.nresponses;
names = {'misses','lures','FAs','maintenance','updates','responses'};
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
    proportions(n,1) = p < 0.05;
    [r,p] = corr(onesubj.maintained,onesubj.nFAs,'type','Spearman');
    proportions(n,2) = p < 0.05;
    [r,p] = corr(onesubj.nlures,onesubj.nFAs,'type','Spearman');
    proportions(n,3) = p < 0.05;
    [r,p] = corr(onesubj.maintained,onesubj.nupdates,'type','Spearman');
    proportions(n,4) = p < 0.05;
    [r,p] = corr(onesubj.nFAs,onesubj.nupdates,'type','Spearman');
    proportions(n,5) = p < 0.05;
    [r,p] = corr(onesubj.nlures,onesubj.nupdates,'type','Spearman');
    proportions(n,6) = p < 0.05;
    % get percentage of subjects w/ significant correlation of all these
    % components
end

percentages = sum(proportions)/nsubjs;

%% When are false alarm errors being made?
% During 2-back only? During other tasks too?

figure
subplot(2,2,1)
histogram(toanalyze.nFAs(toanalyze.task==1,:),'DisplayName','1-detect')
title('1-detect false alarms')
ylabel('Frequency')
clean_fig()

subplot(2,2,2)
histogram(toanalyze.nFAs(toanalyze.task==2,:),'DisplayName','1-back')
title('1-back false alarms')
clean_fig()
 
subplot(2,2,3)
histogram(toanalyze.nFAs(toanalyze.task==3,:),'DisplayName','2-back')
title('2-back false alarms')
ylabel('Frequency')
xlabel('# false alarms per round')
clean_fig()

subplot(2,2,4)
histogram(toanalyze.nFAs(toanalyze.task==4,:),'DisplayName','3-detect')
title('3-detect false alarms')
xlabel('# false alarms per round')
clean_fig()

totalFAs = sum(toanalyze.nFAs);
vector = [sum(toanalyze.nFAs(toanalyze.task==1,:)) sum(toanalyze.nFAs(toanalyze.task==2,:)) ...
    sum(toanalyze.nFAs(toanalyze.task==3,:)) sum(toanalyze.nFAs(toanalyze.task==4,:))];
percentage = vector./totalFAs;
% percentage of total FAs in 1-detect, 1-back, 2-back, and 3-detect ^^

%% Do lure trials elicit false alarm errors?
% Above analysis suggests yes - as they *both* only really occur in the
% 2-back task.

% How often per subject in each task?
for sii = 1:nsubjs
    
    % pull each subject's amalgamated data (from toanalyze)
    onesubj = toanalyze(toanalyze.subj==sii,:);
    
    three_detect = onesubj.nFAs(onesubj.task==4,:);
    % this subject successfully avoided this task ;)
    if isempty(three_detect)
        three_detect = [0 0];
    end
    ntrials_three_detect(sii) = sum(onesubj.task==4); %nrounds x trials per round
    per_round_three_detect(sii) = sum(three_detect)./ntrials_three_detect(sii);
    min_max_three_detect(sii,:) = [min(three_detect) max(three_detect)];
    
    two_back = onesubj.nFAs(onesubj.task==3,:);
    if isempty(two_back)
        two_back = [0 0];
    end
    ntrials_two_back(sii) = sum(onesubj.task==3);
    per_round_two_back(sii) = sum(two_back)./ntrials_two_back(sii);
    min_max_two_back(sii,:) = [min(two_back) max(two_back)];
    
    one_back = onesubj.nFAs(onesubj.task==3,:);
    if isempty(one_back)
        one_back = [0 0];
    end
    ntrials_one_back(sii) = sum(onesubj.task==3);
    per_round_one_back(sii) = sum(one_back)./ntrials_one_back(sii);
    min_max_one_back(sii,:) = [min(one_back) max(one_back)];
    
end

%% What about other possible types of lures? (Noted by reviewer #2)
% To look at that, we need more specific trial-by-trial data.

load([data_directory 'simdata.mat'])

% Goes through the tasks round by round to tally the lures and 
% the false alarm errors that occur in response to those lures.
% We are IGNORING the 3-detect because false alarm responses are so few
% during that task.

for subj = 1:length(unique(tosim.subj)) % loop over subjects captured in tosim
    
    % pull experimental data (trial-by-trial) for one subject at a time
    oneperson = tosim(tosim.subj==subj,:);

    for block = 0:max(oneperson.block)
        
        % now, for each subject, pull the data from each block of each task
        blockdata = oneperson(oneperson.block==block,:);
        bii = block+1; %can't 0 index in matlab!
        
        %initialize lure analysis variables
         
        % all false alarms this block
        nFAs_all = sum((blockdata.detect|blockdata.nback)==false & ~isnan(blockdata.keypress)); 
        
        nFAs_lure(bii,:) = 0;
        nFAs_3backlure(bii,:) = 0;
        nFAs_nolure(bii,:) = 0;
        
        nFAs_1back_2backlure(bii,:) = 0;
        nFAs_1back_3backlure(bii,:) = 0;
        nFAs_1back_nolure(bii,:) = 0;

        
        nlures_2back(bii,1) = 0;
        nlures_1back(bii,1) = 0;
        
        nFAs_per_lure_1back(bii,1) = 0;
        nFAs_per_lure_2back(bii,1) = 0;

        
        %keep track of storage over time, to get a sense of maintained info over time
        task = blockdata.task(2)+1;
        
        storage = NaN(1,3); %empty storage for n-back
        
        for trial = 2:height(blockdata)
            
            % stimulus on that trial
            stim = blockdata.stimnum(trial);
            resp = ~isnan(blockdata.keypress(trial));
            non_match = (blockdata.detect(trial)|blockdata.nback(trial))==false;
            
            % are they making false alarms on these potential lure
            % trials?
            if non_match
                
                if task==3 %only in 2-back can there be WM lures inside buffer
                    
                    % an item 1 trial ago matches the current one
                    if storage(3)==stim
                        nFAs_lure(bii,1) = nFAs_lure(bii,1) + resp; %not a real 2-back, but a lure trial
                        nlures_2back(bii,1) = nlures_2back(bii,1) + 1;                    
                    % an item 3 trials ago matches the current one
                    elseif storage(1)==stim
                        nFAs_3backlure(bii,1) = nFAs_3backlure(bii,1) + resp; %not a real 2-back, but a lure trial
                        nlures_2back(bii,1) = nlures_2back(bii,1) + 1;
                    else
                        nFAs_nolure(bii,1) = nFAs_nolure(bii,1) + resp; %not a real 2-back, but a lure trial
                    end
                    
                end
                
                if task==2 % there are also false alarms during the 1-back task, check those out
                    
                    % an item 2 trials ago matches the current one
                    if storage(2)==stim
                        nFAs_1back_2backlure(bii,1) = nFAs_1back_2backlure(bii,1) + resp; %not a real 2-back, but a lure trial
                        nlures_1back(bii,1) = nlures_1back(bii,1) + 1;                    
                    % an item 3 trials ago matches the current one
                    elseif storage(1)==stim
                        nFAs_1back_3backlure(bii,1) = nFAs_1back_3backlure(bii,1) + resp; %not a real 2-back, but a lure trial
                        nlures_1back(bii,1) = nlures_1back(bii,1) + 1;
                    else
                        nFAs_1back_nolure(bii,1) = nFAs_1back_nolure(bii,1) + resp;
                    end
                end
                
                 %clear out old info, add in new info, shift stuff over
                storage(1) = [];
                storage(end+1) = stim;
                
            end
            
        end
        
        % don't redundantly save, save block according to task identity
        % 2-back
        nFAs_all(nFAs_all == 0) = NaN;
        if task == 3
            
            nFAs_lure(bii,1) = nFAs_lure(bii,1)./nFAs_all;
            
            nFAs_3backlure(bii,1) = nFAs_3backlure(bii,1)./nFAs_all;
            
            nFAs_nolure(bii,1) = nFAs_nolure(bii,1)./nFAs_all;
            
            nFAs_per_lure_2back(bii) = nlures_2back(bii,1)./nFAs_all;
                
        % 1-back
        elseif task == 2
            
            nFAs_1back_2backlure(bii,1) = nFAs_1back_2backlure(bii,1)./nFAs_all;
            
            nFAs_1back_3backlure(bii,1) = nFAs_1back_3backlure(bii,1)./nFAs_all;
            
            nFAs_1back_nolure(bii,1) = nFAs_1back_nolure(bii,1)./nFAs_all;
            
            nFAs_per_lure_1back(bii,1) = nlures_1back(bii,1)./nFAs_all;
            
        end
        
    end %end of block by block loop
    
    avg_per_subj(subj,:) = nanmean([nFAs_lure(:,1) nFAs_3backlure(:,1) nFAs_nolure(:,1) nFAs_1back_2backlure(:,1) nFAs_1back_3backlure(:,1) nFAs_1back_nolure(:,1)]);
    nFAs_per_lure(subj,:) = nanmean([nFAs_per_lure_2back nFAs_per_lure_1back]);

end %of loop over subjects

% what's up with lure trials & how they evoke false alarm errors?
figure
hold on
%plot([0.75 2.5 3.25],repmat(nanmean(nFAs_per_lure(:,1)),1,3),'--b','LineWidth',3)
%plot([3.75 5.5 6.25],repmat(nanmean(nFAs_per_lure(:,2)),1,3),'--r','LineWidth',3)
errorbar([1:6],nanmean(avg_per_subj),nanstd(avg_per_subj)./sqrt(n),'*k','LineWidth',2)
legend({'Mean lures per FA in 2-back','Mean lures per FA in 1-back'})
xticklabels({'Traditional lures','3-back lures in 2-back task','No lure','2-back lures in 1-back task','3-back lures in 1-back task','No lure'})
xtickangle(45)
xticks([1:6])
ylabel('% of total false alarms per lure type')
clean_fig()

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




