%% Analysis of cost of control task data

% THIS SCRIPT DEPENDS ON VARIABLES PULLED FROM PAPER_GRAPHS_AND_STATS.M 
% IT CAN'T BE RUN ON ITS OWN, ONLY NESTED WITHIN PAPER_GRAPHS_AND_STATS.M

% What is this script for?
% Produces all sorts of supplementary figures, NOT the supplementary
% figures in the paper, but data quality checks and other simple analyses
% of task & BDM rating data.


%% Check data quality - how long each task takes, 
% & what sort of task completion rates are we looking at?

% first, double-check that all tasks take the same amount of time, as we
% designed them to
% there was a bug in very early data collection which makes it appear that
% the 1-detect takes really long in 3-4 subjects, but that was fixed
% such that the 1-detect was accurately timed in subjects after that point
figure
subplot(1,3,3)
for subj = 1:n
    y = data.TOT(subj,data.task_progression(subj,:)==tasks(1));
    x = data.TOT(subj,data.task_progression(subj,:)==tasks(2));
    z = data.TOT(subj,data.task_progression(subj,:)==tasks(3));
    a = data.TOT(subj,data.task_progression(subj,:)==tasks(4));
    scatter(ones(length(y),1),y,[],subjcolors(subj,:),'Filled')
    hold on
    scatter(2*ones(length(x),1),x,[],subjcolors(subj,:),'Filled')
    scatter(3*ones(length(z),1),z,[],subjcolors(subj,:),'Filled')
    scatter(4*ones(length(a),1),a,[],subjcolors(subj,:),'Filled')
end
% Fix outliers - where did those come from?
%[h,p] = ttest2(z,y);
xlim([0.75 4.25])
ylabel('Time on task (seconds)')
xticks(1:length(tasklabels))
xticklabels(tasklabels)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;

[h,p] = ttest2(data.TOT(data.task_progression==tasks(1)),data.TOT(data.task_progression==tasks(2)));
[h,p] = ttest2(data.TOT(data.task_progression==tasks(3)),data.TOT(data.task_progression==tasks(1)));
[h,p] = ttest2(data.TOT(data.task_progression==tasks(3)),data.TOT(data.task_progression==tasks(2)));
[h,p] = ttest2(data.TOT(data.task_progression==tasks(1)),data.TOT(data.task_progression==tasks(4)));

errorbar(1:length(tasklabels),[mean(y) mean(x) mean(z) mean(a)],[std(y)/sqrt(length(y)) std(x)/sqrt(length(x)) std(z)/sqrt(length(z)) std(a)/sqrt(length(a))],'k')

% second, how many times was each task completed, across all subjects? how much
% data are we really analyzing/modeling?

subplot(1,3,2)
y = sum(data.taskfreqs(:,1));
x = sum(data.taskfreqs(:,2));
w = sum(data.taskfreqs(:,3));
z = sum(data.taskfreqs(:,4));
%bar([y,x,w,z])
bar([y,x,w,z])
%xticklabels({'detection','combine','n-switch,p0.1','n-switchp0.9'})
xticklabels(tasklabels)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title('# times each task completed across subjects')
ylabel('Count')

% Where in the task is dropout occurring?
subplot(1,3,1)
bar([sum(excluded.values(excluded.exp_version==4,1:4))+n n])
labels = [excluded.labels(1,1:4) 'finished'];
xticklabels([labels])
title('# Subjects completed each phase')
ylabel('n')
ax = gca; fig = gcf;
ax.FontSize = 14;
fig.Color = 'w';


%% Differences in task accuracy, RT, and explicit difficulty ratings (collected at conclusion of expt.)

% Mean accuracy, explicit difficulty ratings, 
% how accuracy changes over subjects across tasks, and mean reaction times
% (RTs)
figure
subplot(2,2,1)
fig = gcf;
fig.Color = 'w';
bar([nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3)) nanmean(tasks_overall(:,4))])
hold on
errorbar(1:length(tasklabels),[nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3)) nanmean(tasks_overall(:,4))],[nanstd(tasks_overall(:,1)) nanstd(tasks_overall(:,2)) nanstd(tasks_overall(:,3)) nanstd(tasks_overall(:,4))]/sqrt(n),'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
ylim([50 100])
title(['Accuracy by task'])

subplot(2,2,2)
fig = gcf;
fig.Color = 'w';
bar([nanmean(data.taskratings(:,1)) nanmean(data.taskratings(:,2)) nanmean(data.taskratings(:,4)) nanmean(data.taskratings(:,3))])
hold on
errorbar(1:length(tasklabels),[nanmean(data.taskratings(:,1)) nanmean(data.taskratings(:,2)) nanmean(data.taskratings(:,4)) nanmean(data.taskratings(:,3))],[nanstd(data.taskratings(:,1)) nanstd(data.taskratings(:,2)) nanstd(data.taskratings(:,4)) nanstd(data.taskratings(:,3))]/sqrt(n),'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
title(['Difficulty rating by task'])

subplot(2,2,3)
allsubjs = [];
for subj = 1:n
    line = [nanmean(data.detectacc(subj,:)) nanmean(data.n1acc(subj,:)) nanmean(data.ndetectacc(subj,:)) nanmean(data.n2acc(subj,:))];
    allsubjs = [allsubjs; line];
    plot(1:length(line),line,'--')
    hold on
end
errorbar(nanmean(allsubjs),nanstd(allsubjs)/sqrt(n),'k','LineWidth',1.5)
title('Accuracy subj by subj for each task')
xticklabels(tasklabels)
xticks(1:length(tasklabels))
fig = gcf;
fig.Color = 'w';
xlim([0.75 4.25])
missingvals = sum(isnan(allsubjs),2)>0;
allsubjs(missingvals,:) = []; %prune for stats
[r,p] = corr(allsubjs(:,1),allsubjs(:,2));
[r,p] = corr(allsubjs(:,2),allsubjs(:,3));
[r,p] = corr(allsubjs(:,1),allsubjs(:,3));
[r,p] = corr(allsubjs(:,1),allsubjs(:,4));

subplot(2,2,4)
fig = gcf;
fig.Color = 'w';
bar([nanmean(tasks_rts(:,1)) nanmean(tasks_rts(:,2)) nanmean(tasks_rts(:,4)) nanmean(tasks_rts(:,3))])
hold on
errorbar(1:length(tasklabels),[nanmean(tasks_rts(:,1)) nanmean(tasks_rts(:,2)) nanmean(tasks_rts(:,4)) nanmean(tasks_rts(:,3))],[nanstd(tasks_rts(:,1)) nanstd(tasks_rts(:,2)) nanstd(tasks_rts(:,4)) nanstd(tasks_rts(:,3))]/sqrt(n),'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
title('Mean RT by task')

%% Examine overall trends in fair wage (BDM) ratings
% and in task performance over time

% do the BDM's evolve? do BDM RTs evolve?
% does task performance evolve?

%just check out subject BDM strategy
figure
subplot(3,2,1)
for i = 1:n
    scatter(1:default_length,data.values(i,:),'o','Filled')
    hold on
end
%errorbar(nanmean(data.values),nanstd(data.values)/sqrt(n),'k','LineWidth',1)
title('Mean fair wage per subject per block')
fig = gcf; ax = gca;
fig.Color = 'w';
legend({'Subj 1','Subj 2','...','Subj n'})
ylabel('BDM points requested')
xlabel('Block')

subplot(3,2,2)
n1 = data.task_displayed==tasknumbers(2); %1
n2 = data.task_displayed==tasknumbers(4); %2
ndetect = data.task_displayed==tasknumbers(3); %7
% bar([nanmean(data.values(n1)) nanmean(data.values(n2)) nanmean(data.values(ndetect))],'BarWidth',0.5)
% for subj = 1:n %plot subj by subj
%     errorbar([nanmean(data.values(subj,n1(subj,:))) nanmean(data.values(subj,n2(subj,:))) nanmean(data.values(subj,ndetect(subj,:)))],[nanstd(data.values(subj,n1(subj,:)))/sqrt(sum(n1(subj,:))) nanstd(data.values(subj,n2(subj,:)))/sqrt(sum(n2(subj,:))) nanstd(data.values(subj,ndetect(subj,:)))/sqrt(sum(ndetect(subj,:)))],'LineWidth',0.75)
%     hold on
% end
errorbar([nanmean(data.values(n1)) nanmean(data.values(ndetect)) nanmean(data.values(n2))],[nanstd(data.values(n1)) nanstd(data.values(ndetect)) nanstd(data.values(n2))]./sqrt(n),'k','LineWidth',1.25)
xticks(1:(length(tasklabels)-1))
xlim([0.75 3.25])
ylim([1 5])
xticklabels(tasklabels(2:end))
title('Mean fair wage by task')
% legend({'Task Means','Subj 1','Subj 2','...','Subj n'})
xlabel('Task')

[h,p] = ttest2(data.values(n1),data.values(n2))
disp ('t-test task 1 versus task 2 BDM values'); %task 1 values versus task 2 values
[h,p] = ttest2(data.values(ndetect),data.values(n2))
disp ('t-test task 2 versus task 3 BDM values'); %task 2 values versus task 3 values
[h,p] = ttest2(data.values(n1),data.values(ndetect))
disp ('t-test task 1 versus task 3 BDM values'); %task 1 values versus task 3 values
% basically it's task 2 versus the world

subplot(3,2,3)
n1subjlearning = NaN(n,default_length+1); n2subjlearning = NaN(n,default_length+1); n3subjlearning = NaN(n,default_length+1);
n1subjvalue = NaN(n,default_length+1); n2subjvalue = NaN(n,default_length+1); n3subjvalue = NaN(n,default_length+1);
n1subjrt = NaN(n,default_length+1); n2subjrt = NaN(n,default_length+1); n3subjrt = NaN(n,default_length+1);
for row = 1:n %cycle through subjects
    for task = 1:3 %cycle through tasks
        for trial = 2:default_length
            iters = sum(data.task_progression(row,1:(trial-1)) == tasks(task+1))+1;
            if data.task_progression(row,trial) == tasks(task+1)
                eval(['n' num2str(task) 'subjlearning(row,iters) = data.perf(row,trial);'])
                eval(['n' num2str(task) 'subjvalue(row,iters) = data.values(row,trial);'])
                eval(['n' num2str(task) 'subjrt(row,iters) = data.BDMrt(row,trial);'])
            end
        end
    end
end
errorbar(nanmean(n1subjlearning),nanstd(n1subjlearning)/sqrt(n),'Color',taskcolors(1,:),'LineWidth',2)
hold on
errorbar(nanmean(n2subjlearning),nanstd(n2subjlearning)/sqrt(n),'Color',taskcolors(2,:),'LineWidth',2)
errorbar(nanmean(n3subjlearning),nanstd(n3subjlearning)/sqrt(n),'Color',taskcolors(3,:),'LineWidth',2) %task 3 is actually ndetect
title('Task performance per task iteration')
legend(tasklabels(2:end))
fig = gcf; ax = gca;
fig.Color = 'w';
ylabel('Accuracy')
xlabel('Task iteration')
xticks([1:21])
xticklabels([0:20])

subplot(3,2,4)
errorbar(nanmean(n1subjvalue),nanstd(n1subjvalue)/sqrt(n),'Color',taskcolors(1,:),'LineWidth',2)
hold on
errorbar(nanmean(n2subjvalue),nanstd(n2subjvalue)/sqrt(n),'Color',taskcolors(2,:),'LineWidth',2)
errorbar(nanmean(n3subjvalue),nanstd(n3subjvalue)/sqrt(n),'Color',taskcolors(3,:),'LineWidth',2)
title('BDM by task iteration')
legend(tasklabels(2:end))
fig = gcf; ax = gca;
fig.Color = 'w';
ylabel('BDM points requested')
xlabel('Task iteration')
xticks([1:21])
xticklabels([0:20])

subplot(3,2,5)
% are people getting more decisive on how many BDM points they want?
errorbar(nanmean(data.BDMrt),nanstd(data.BDMrt)/sqrt(n),'LineWidth',1.5)
title('BDM RT by block')
ylabel('RT on BDM choice')
xlabel('Block')
fig = gcf; fig.Color = 'w';

subplot(3,2,6)
% are people getting more decisive on how many BDM points they want?
for task = 1:(length(unique(data.task_displayed(~isnan(data.task_displayed)))))
    matrix = [];
    for subj = 1:n
        eval(['rts = n' num2str(task) 'subjrt(subj,:);'])
        display = data.task_displayed(subj,:);
        init = NaN(1,default_length);
        init(~isnan(rts)) = rts(~isnan(rts));
        matrix = [matrix; init];
    end
    errorbar(nanmean(matrix),nanstd(matrix)/sqrt(n),'Color',taskcolors(task,:),'LineWidth',1.5)
    hold on
end
title('BDM RT by iteration by task displayed')
ylabel('RT on BDM choice')
xlabel('Task iteration')
xticks([1:21])
xticklabels([0:20])
legend(tasklabels(2:end))
fig = gcf; fig.Color = 'w';

%% MEAN BDM RATINGS NOT DIFFERENT BETWEEN 1-BACK and 3-DETECT
% Is that true within-subject, too, or just on the group-level>
% Seems like it's just on the group-level, and subjects are not treating
% these tasks exactly the same way.

%are there individual differences in n1 vs ndetect valuation?
%there aren't group differences
figure
subplot(1,2,1)
for subj = 1:n
    values = data.values(subj,:);
    n1 = data.task_displayed(subj,:) == 1; ndetect = data.task_displayed(subj,:) == 7;
    errorbar(1:2,[nanmean(values(n1)) nanmean(values(ndetect))],[nanstd(values(n1)) nanstd(values(ndetect))]./[sqrt(sum(n1)) sqrt(sum(ndetect))],'Color',subjcolors(subj,:),'LineWidth',1)
    hold on
end
xticks([1:2])
xticklabels({'1-back','3-detect'})
xlim([0.75 2.25])
ylabel('Mean BDM value')
title('Differences in cost of 1-back and 3-detect?')

bars = [];
for row = 1:length(data.taskBDMcorrs)
    bars = [bars; data.taskBDMcorrs{row,2}(1) data.taskBDMcorrs{row,2}(3) data.taskBDMcorrs{row,1}(1,2)];
end
subplot(1,2,2)
bar([nanmean(bars(:,1)) nanmean(bars(:,2)) nanmean(bars(:,3))])
hold on
errorbar([nanmean(bars(:,1)) nanmean(bars(:,2)) nanmean(bars(:,3))], [nanstd(bars(:,1)) nanstd(bars(:,2)) nanstd(bars(:,3))]/sqrt(n),'*k','LineWidth',1.5)
title('BDM value correlations')
xticklabels({'Autocorr 1-back','Autocorr 3-detect','Cross corr'})
xtickangle(45)
fig = gcf; fig.Color = 'w';

%% LEARNING CURVES
% plot task performance by block, any visible learning or decay?
% (No.)

figure
subplot(1,4,1)
ax = gca; fig = gcf;
for i = 1:n
    plot(data.detectacc(i,:),'ok')
    hold on
end
nblocks = length(data.detectacc(i,:));
errorbar(nanmean(data.detectacc),nanstd(data.detectacc)/sqrt(n),'k','LineWidth',1)
for block = 1:nblocks
    [h,p] = ttest(data.detectacc(:,block),100);
    if h == 1
        plot(block,75,'k*','LineWidth',1)
    end
end
title('Learning Curves (0-back)')
fig.Color = 'w';
ax.FontSize = 12;
ylim([30 109])
xlabel('Block #')
ylabel('Accuracy')

%test statistical deviation from ceiling on specific blocks
for i = 1:32
    [h,p] = ttest(100.*ones(n,1),data.detectacc(:,i));
    if p < 0.05
        disp(['block ' num2str(i) ', p = ' num2str(p)])
    end
end

subplot(1,4,2)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.n1acc(i,:),'o','Color',taskcolors(1,:))
end
errorbar(nanmean(data.n1acc),nanstd(data.n1acc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (1-back)')
fig.Color = 'w';
ax.FontSize = 12;
ylim([30 109])
xlabel('Block #')
ylabel('Accuracy')

subplot(1,4,3)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.ndetectacc(i,:),'o','Color',taskcolors(3,:))
end
errorbar(nanmean(data.ndetectacc),nanstd(data.ndetectacc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (3-detect)')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylim([30 109])
ylabel('Accuracy')

subplot(1,4,4)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.n2acc(i,:),'o','Color',taskcolors(2,:))
end
errorbar(nanmean(data.n2acc),nanstd(data.n2acc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (2-back)')
fig.Color = 'w';
ax.FontSize = 12;
ylim([30 109])
xlabel('Block #')
ylabel('Accuracy')
linkaxes

%% How do people rate these tasks right before they complete them?
% Based on task iteration, what's the mean BDM?
% Grouping subjects with similar numbers of task iterations

figure
for task = 1:3
    subplot(3,1,task)
    rating_idx = data.task_displayed==tasknumbers(task+1);
    eval(['ratings = n' num2str(task) 'subjvalue;'])
    hold on
    eval(['completions = sum(~isnan(n' num2str(task) 'subjlearning),2);'])
    ratings = [completions ratings]; ratings = sortrows(ratings,1);
    n_rounds = unique(completions);
    for i = 1:length(n_rounds)
        tomean = ratings(ratings(:,1)==n_rounds(i),2:end);
        toplot = nanmean(tomean,1);
        errorbar(toplot,nanstd(tomean,[],1)./sqrt(size(tomean,1)),'Color',taskcolors(task,:))
        hold on
    end
    xlabel('# task iters completed')
    ylabel('Mean BDM rating')
    title(tasklabels{task+1})
end
fig = gcf; fig.Color = 'w';


%% Plot stuff related to BDM gaming etc.
% are people tracking the offer at stake, and adjusting performance to it?
% seems like no.

figure
subplot(1,2,1)
for row = 1:n
    init = NaN(default_length,1);
    init(data.offers(row,:)>data.values(row,:)) = data.offers(row,data.offers(row,:)>data.values(row,:));
    init(data.values(row,:)>data.offers(row,:)) = 1;
    y = init;
    x = [NaN data.perf(row,1:end-1)]';
    matrix = sortrows([y,x],1);
    plot(matrix(:,1),matrix(:,2),'o')
    hold on
end
ax = gca; fig = gcf;
title('Performance by BDM offer')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('BDM points at stake')
ylabel('Accuracy')

% are subjects matching the previous offers to try and manipulate the
% points system? (it wouldn't work, but is a possible strategy if you think
% the computer isn't random & is maybe going lower on its bids)

%plot BDM request by previous offer
all_subj_matrix = [];

subplot(1,2,2)
for row = 1:n
    y = data.offers(row,:)';
    x = [NaN data.values(row,2:end)]';
    matrix = sortrows([y,x],1);
    plot(matrix(:,1),matrix(:,2),'o')
    hold on
    matrix(sum(isnan(matrix),2)>0,:) = [];
    [r,ps(row)] = corr(matrix(:,1),matrix(:,2),'type','Spearman');
    all_subj_matrix = [all_subj_matrix; matrix];
end
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title('BDM value by prev. computer offer')
ylim([1 5.1])
ylabel('BDM points')
xlabel('Last offer')

% big correlation
[r,p] = corr(all_subj_matrix(:,1),all_subj_matrix(:,2),'type','Spearman');

% how many correlations come out on individual level?
n_sig_individual_corr = sum(ps<0.05);

%% Plot late responses/changed responses by task

% are subjects more inattentive in one task than another?%
% not enough data on this

figure
for task = 1:length(tasks)
    bar(task,nanmean(data.lateresponse(data.task_progression==tasks(task))))
    scatter(task*ones(sum(sum(data.task_progression==tasks(task))),1),data.lateresponse(data.task_progression==tasks(task)))
    hold on
end
ax = gca; fig = gcf;
xticks(1:length(tasks))
xticklabels(tasklabels)
xtickangle(45)
ylabel('mean late responses')
ax.FontSize = 12;
fig.Color = 'w';


%% Understand n-back matches' effect on costs
% very noisy to plot, not very informative
% the correlation analyses are more useful, but essentially fruitless

% this analysis was a first stab at understanding which components were
% costly, pre-process model - just pull matches, and see whether BDM
% requests went up in response to them. the number of matches was set each
% round to be 3, 4, or 5.

%first, plot general nbackmatches effect for all subjects
figure
matrix = [];
for i=1:n
    task_list = data.task_progression(i,:);
    nback = find(task_list==tasks(2)|task_list==tasks(3)&~inattentive(i,:));
    next = nback+1;
    l = length(task_list);
    y = data.values(i,next(next~=(l+1)));
    x = data.nbackmatches(i,nback(nback~=l)); %exclude last trial for sizing reasons
    %costs = NFC_sigmoid(x,data.NFC(i));
    matrix = [matrix; y' x'];
    color = subjcolors(i,:);
    scatter(x,y,[],color,'o','Filled');
    hold on
end
title('BDM value by prev. # n-back matches')
ylabel('BDM points')
xlabel('# of n-back matches in last block')
ax = gca; fig = gcf;
legend({'Subj 1','Subj 2','...','Subj n'})
fig.Color = 'w';
ax.FontSize = 12;
ylim([1 5.1])

%is this significant?
x = matrix(:,1);
y = matrix(:,2);
notvalid = isnan(y) | isnan(x);
x(notvalid) = []; y(notvalid) = [];
[r,p] = corr(x,y,'Type','Spearman');
%disp(['Nswitches vs. BDM request. \ Spearman rho: ' num2str(r) ', p value: ' num2str(p)])

same_long = []; diff_long = [];
figure
for task = 1:(length(tasklabels)-1) %cycle through what is being displayed, display by what happened before
    subplot(2,3,task)
    for subj = 1:n %go subject by subject
        displayed = data.task_displayed(subj,:)+1;
        idx = displayed == (task+1);
        matches = [NaN data.nbackmatches(subj,:)]; matches(end) = []; %last trial not influential
        BDMs = data.values(subj,:);
        
        % take stock of which task was actually completed when (stochastic
        % task completion because of rating procedure means the rated task is
        % not always the completed task on that round)
        completed = NaN(1,default_length);
        completed(data.task_progression(subj,:)==tasks(2)) = 2;
        completed(data.task_progression(subj,:)==tasks(3)) = 3;
        completed(data.task_progression(subj,:)==tasks(4)) = 8;
        completed(:,end) = []; completed = [NaN completed];%delete last column, but conserve size & shift numbers down w/ a NaN
        
        same = completed==displayed;
        different = ~isnan(completed)&completed~=displayed; %different but not detection
        same_long = [same_long; matches(same&idx)' BDMs(same&idx)'];
        diff_long = [diff_long; matches(different&idx)' BDMs(different&idx)'];
    end
    toplot = count_entries(diff_long);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'c','Filled')
    hold on
    toplot = count_entries(same_long);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'m','Filled')
    legend({'Diff task','Same Task'},'Location','Best')
    ylabel('BDM points')
    title(['BDM for ' tasklabels(task+1) ' vs. matches in prev round'])
    xlabel('# of n-back matches in last block')
end
fig = gcf;
fig.Color = 'w';
ylim([1 5.1])

diff_stats = diff_long; diff_stats(isnan(diff_stats(:,1)),:) = []; %prune to run a correlation
diff_stats(isnan(diff_stats(:,2)),:) = [];
[r,p] = corr(diff_stats(:,1),diff_stats(:,2));
disp(['Infl of prev matches in different task on current task BDM - r: ' num2str(r) ', p: ' num2str(p)])

same_stats = same_long; same_stats(isnan(same_stats(:,1)),:) = []; %prune to run a correlation
same_stats(isnan(same_stats(:,2)),:) = [];
[r,p] = corr(same_stats(:,1),same_stats(:,2));
disp(['Infl of prev matches in same task on current task BDM - r: ' num2str(r) ', p: ' num2str(p)])

n1effect = []; n2effect = []; n3effect = []; % keep track of numbers pulled out here for stats later
for task = 1:(length(tasklabels)-1) %cycle through what is being displayed, display by what happened before
    subplot(2,3,task+3)
    for subj = 1:n %go subject by subject
        displayed = data.task_displayed(subj,:);
        idx = find(displayed == (tasknumbers(task+1)));
        matches = data.nbackmatches(subj,:);
        misses = data.nbackmisses(subj,:);
        if task == 3
            matches = data.ndetectmatches(subj,:);
            misses = data.ndetectmisses(subj,:);
        end
        perf = data.perf(subj,:);
        BDMs = data.values(subj,:);
        completed = find(data.task_progression(subj,:)==tasks(task+1));
        color = subjcolors(subj,:);
        for trial = 1:length(idx)
            now = idx(trial);
            if sum(completed<now)>0 %they've done the last once before
                last = completed(completed<now); last = last(end);
                delay = now-last;
                scatter(matches(last),BDMs(now),[],color,'Filled') %do this subject by subject instead
                eval(['n' num2str(task) 'effect = [n' num2str(task) 'effect; subj matches(last) misses(now) delay BDMs(now) perf(now)];'])
            end
            hold on
        end
    end
    ylabel('BDM points')
    title(['BDM for ' tasklabels(task+1)])
    xlabel('matches in last block')
end
fig = gcf;
fig.Color = 'w';
ylim([1 5.1])

%% How does delay between task iterations affect behavior, or BDM ratings?

% are people forgetting their previous task experience?
% i.e. is there decay of BDM values when the time between task rating
% iterations is increased (measured by proxy, looking at the number of
% intervening task ratings)?
% what about task accuracy?
% or task RT?
figure
subplot(3,1,1)
toplot = count_entries([n2effect(:,4) n2effect(:,6)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
hold on
toplot = count_entries([n1effect(:,4) n1effect(:,6)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
legend({'2-back','1-back'})
ylabel('Accuracy')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';

subplot(3,1,2)
toplot = count_entries([n2effect(:,4:5)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
hold on
toplot = count_entries([n1effect(:,4:5)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
legend({'2-back','1-back'})
ylabel('BDM request')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';

subplot(3,1,3)
%effect structures are 1. subj 2. matches 3. misses 4. delay 5. BDM 6. perf
toplot = count_entries([n2effect(:,3:4)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
hold on
toplot = count_entries([n1effect(:,3:4)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
legend({'2-back','1-back'})
ylabel('Missed matches')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';


%% Some interesting plots dealing with individual differences measures,
% Need for Cognition (NFC) and Short Almost Perfect Scale (SAPS) responses,
% already scored.

figure
%tertile split SAPS scores
split = tertile_split(data.SAPS);
%plot accuracy by SAPS score

% are highly perfectionistic subjects actually better at the tasks?
% or does perfectionism not translate to better task performance?

subplot(1,3,1)
bars = [nanmean(tasks_overall(split==1,1)) nanmean(tasks_overall(split==2,1)) nanmean(tasks_overall(split==3,1)); ...
    nanmean(tasks_overall(split==1,2)) nanmean(tasks_overall(split==2,2)) nanmean(tasks_overall(split==3,2)); ...
    nanmean(tasks_overall(split==1,3)) nanmean(tasks_overall(split==2,3)) nanmean(tasks_overall(split==3,3)); ...
    nanmean(tasks_overall(split==1,4)) nanmean(tasks_overall(split==2,4)) nanmean(tasks_overall(split==3,4))];
E = [nanstd(tasks_overall(split==1,1)) nanstd(tasks_overall(split==2,1)) nanstd(tasks_overall(split==3,1)); ...
    nanstd(tasks_overall(split==1,2)) nanstd(tasks_overall(split==2,2)) nanstd(tasks_overall(split==3,2)); ...
    nanstd(tasks_overall(split==1,3)) nanstd(tasks_overall(split==2,3)) nanstd(tasks_overall(split==3,3)); ...
    nanstd(tasks_overall(split==1,4)) nanstd(tasks_overall(split==2,4)) nanstd(tasks_overall(split==3,4))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:12,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:12,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:12,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
title('Accuracy by SAPS group')
xticks([2:3:12])
legend({'Low SAPS','Mid SAPS','High SAPS'})
xticklabels(tasklabels)
ylim([50 100])
ylabel('Mean Accuracy')
fig = gcf; fig.Color = 'w';
%diff in accuracy in high and low SAPS groups?
[h,p] = ttest2(data.overall(split==1),data.overall(split==3)); % no
clean_fig();


% looks like highly perfectionist subjects are LESS accurate in their task
% performance, especially on 2-back task. how does that effect their BDM
% fair wage ratings?

% maybe they ask for higher wages after errors, because the errors are so
% aversive to them?

subplot(1,3,2)
low = []; mid = []; high = [];
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2)|data.task_progression(row,:)==tasks(4));
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.perf(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==1
        low = [low; x' y'];
    elseif split(row)==2
        mid = [mid; x' y'];
    elseif split(row)==3
        high = [high; x' y'];
    end
end
scatter(low(:,1),low(:,2),[],SAPScolors(1,:),'Filled')
hold on
scatter(mid(:,1),mid(:,2),[],SAPScolors(2,:),'Filled')
scatter(high(:,1),high(:,2),[],SAPScolors(3,:),'Filled')
low(isnan(low(:,2)),:) = []; mid(isnan(mid(:,2)),:) = []; high(isnan(high(:,2)),:) = [];
[r,p] = corr(low(:,1),low(:,2),'Type','Spearman');
[r,p] = corr(mid(:,1),mid(:,2),'Type','Spearman');
[r,p] = corr(high(:,1),high(:,2),'Type','Spearman');
legend({'Low SAPS','Mid SAPS','High SAPS'},'Location','Best')
legend('boxoff')
xlabel('Accuracy')
ylabel('BDM value')
xlim([50 100])
clean_fig();

% same logic as above for the following analysis, but now looking at the
% frequency of ommission errors specifically
% it kind of looks vaguely like they're asking for fewer points when their
% performance is lower, almost apologetically

subplot(1,3,3)
low = []; mid = []; high = [];
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2)|data.task_progression(row,:)==tasks(4));
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.nbackmisses(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==1
        low = [low; x' y'];
    elseif split(row)==2
        mid = [mid; x' y'];
    elseif split(row)==3
        high = [high; x' y'];
    end
end
scatter(low(:,1),low(:,2),[],SAPScolors(1,:),'Filled')
hold on
scatter(mid(:,1),mid(:,2),[],SAPScolors(2,:),'Filled')
scatter(high(:,1),high(:,2),[],SAPScolors(3,:),'Filled')
low(isnan(low(:,2)),:) = []; mid(isnan(mid(:,2)),:) = []; high(isnan(high(:,2)),:) = [];
[r,p] = corr(low(:,1),low(:,2),'Type','Spearman');
[r,p] = corr(mid(:,1),mid(:,2),'Type','Spearman');
[r,p] = corr(high(:,1),high(:,2),'Type','Spearman');
legend({'Low SAPS','Mid SAPS','High SAPS'},'Location','Best')
legend('boxoff')
xlabel('Missed Matches')
ylabel('BDM value')
clean_fig();


%% What does behavior look like after errors are made?
% how do perfectionistic subjects bounce back? do they get back on track,
% or does their overall error frequency go up?
% they seem to slow down, but not sure whether they get more or less
% accurate

%make this figure for mean post-error ER and RT
ys{1} = data.posterrorER;
ys{2} = data.posterrorRT;
datalabels = {'ER','RT'};
for col = 1:length(ys)
    figure
    for task = 1:(length(tasklabels)-1)
        y = ys{col};
        relevant = y(data.task_progression==tasks(task+1));
        subplot(3,3,task)
        histogram(relevant)
        title(['Distribution of ' datalabels{col} ' after first error - ' tasklabels(task+1)])
        xlabel(['Post-Error' datalabels{col}])
        subplot(3,3,task+3)
        for subj = 1:n
            specific = y(subj,:);
            relevant = specific(data.task_progression(subj,:)==tasks(task+1));
            scatter(data.SAPS(subj),nanmean(relevant),[],SAPScolors(2,:),'Filled')
            hold on
        end
        title([ datalabels{col} ' after first error - ' tasklabels(task+1)])
        ylabel(['Mean Post-Error ' datalabels{col}])
        xlabel('Perfectionism')
        clean_fig();
        subplot(3,3,task+6)
        for subj = 1:n
            specific = y(subj,:);
            relevant = specific(data.task_progression(subj,:)==tasks(task+1));
            scatter(data.NFC(subj),nanmean(relevant),[],NFCcolors(2,:),'Filled')
            hold on
        end
        title([ datalabels{col} ' after first error - ' tasklabels(task+1)])
        ylabel(['Mean Post-Error ' datalabels{col}])
        xlabel('NFC')
        clean_fig();
    end
end

%% Perfectionism and mean RT on each task
% Maybe perfectionist subjects are more careful during task completion,
% resulting in overall slower RTs?

figure
bars = [nanmean(tasks_rts(split==1,1)) nanmean(tasks_rts(split==2,1)) nanmean(tasks_rts(split==3,1)); ...
    nanmean(tasks_rts(split==1,2)) nanmean(tasks_rts(split==2,2)) nanmean(tasks_rts(split==3,2)); ...
    nanmean(tasks_rts(split==1,3)) nanmean(tasks_rts(split==2,3)) nanmean(tasks_rts(split==3,3)); ...
    nanmean(tasks_rts(split==1,4)) nanmean(tasks_rts(split==2,4)) nanmean(tasks_rts(split==3,4))];
E = [nanstd(tasks_rts(split==1,1)) nanstd(tasks_rts(split==2,1)) nanstd(tasks_rts(split==3,1)); ...
    nanstd(tasks_rts(split==1,2)) nanstd(tasks_rts(split==2,2)) nanstd(tasks_rts(split==3,2)); ...
    nanstd(tasks_rts(split==1,3)) nanstd(tasks_rts(split==2,3)) nanstd(tasks_rts(split==3,3)); ...
    nanstd(tasks_rts(split==1,4)) nanstd(tasks_rts(split==2,4)) nanstd(tasks_rts(split==3,4))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:12,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:12,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:12,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticks([2:3:12])
title('Mean RTs by SAPS group')
legend({'Low SAPS','Mid SAPS','High SAPS'})
xticklabels(tasklabels)
clean_fig();

%% Are perfectionistic subjects better from the start of the task (practice rounds)?
% Are they more confident? Do they find certain tasks more difficult
% (explicit difficulty ratings)?

% Plot perfectionism tertiles by practice accuracy on 2-back task
figure
subplot(1,3,1)
scatter(data.practiceacc(:,4),data.SAPS,[],SAPScolors(2,:),'Filled')
[r,p] = corr(data.practiceacc(~isnan(data.SAPS),4),data.SAPS(~isnan(data.SAPS)));
xlabel('Practice Accuracy on 2-back')
ylabel('Perfectionism')
ax = gca; fig = gcf;
fig.Color = 'w'; ax.FontSize = 12;

% look into self-rated competence i.e. how "easy" or "hard" they rated the
% whole experiment to be
subplot(1,3,2)
E = []; means = [];
for group = 1:3
    relevant = data.diffrating(split==group,:);
    means(group) = nanmean(relevant); E(group) = nanstd(relevant)./sqrt(sum(split==group));
    bar(group,nanmean(relevant),'FaceColor',SAPScolors(group,:))
    hold on
    title(['Overall difficulty rating by SAPS group'])
    ylabel('Mean rating (1 lowest, 5 highest)')
    xticks([1:3])
    xticklabels({'Low SAPS','Mid SAPS','High SAPS'})
end
errorbar(means,E,'k','LineWidth',1.5)
fig = gcf; fig.Color = 'w';

% now break into sub-tasks
subplot(1,3,3)
E = []; means = [];
for task = 1:length(tasklabels)
    for group = 1:3
        relevant = data.taskratings(split==group,:);
        means(group,:) = nanmean(relevant); E(group,:) = nanstd(relevant)./sqrt(sum(split==group));
    end
end
bar(1:3:12,means(1,:),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:12,means(2,:),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:12,means(3,:),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
means = reshape(means,numel(means),1); E = reshape(E,numel(means),1);
errorbar(means,E,'*k')
title('Overall difficulty rating by SAPS group')
ylabel('Mean rating (1 lowest, 5 highest)')
legend({'Low SAPS','Mid SAPS','High SAPS'})
xticklabels(tasklabels)
xticks([2:3:12])
fig = gcf; fig.Color = 'w';

%% Now let's dive into NFC, instead of perfectionism

% are subjects with high NFC more responsive to the number of matches per
% round?

% individual differences in NFC and cognition
split = tertile_split(data.NFC);
%plot nswitches vs. BDM as function of NFC score
low = []; mid = []; high = [];
figure
subplot(1,3,1)
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
    nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.nbackmatches(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==1
        scatter(x,y,'o','Filled')
        low = [low; x' y'];
    end
    hold on
end
title('Low NFC subjects')
ax = gca; fig = gcf;
ax.FontSize = 12;
fig.Color = 'w';
xlabel('N-back Matches')
ylabel('BDM value')
ylim([1 5.1])
[r,p] = corr(low(~isnan(low(:,2)),1),low(~isnan(low(:,2)),2),'Type','Spearman'); %low NFC has relationship p = 0.03 with n = 37
if p<0.05
    y = r.*(unique(low(:,1)))+1;
    plot(unique(low(:,1)),y,'k','LineWidth',1.5)
end

subplot(1,3,2)
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
    nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.nbackmatches(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==2
        scatter(x,y,'o','Filled')
        mid = [mid; x' y'];
    end
    hold on
end
title('Mid NFC subjects')
ax = gca; fig = gcf;
ax.FontSize = 12;
fig.Color = 'w';
xlabel('N-back Matches')
ylabel('BDM value')
ylim([1 5.1])
[r,p] = corr(mid(~isnan(mid(:,2)),1),mid(~isnan(mid(:,2)),2),'Type','Spearman'); %low NFC has relationship p = 0.03 with n = 37
if p<0.05
    y = r.*(unique(mid(:,1)))+2;
    plot(unique(mid(:,1)),y,'k','LineWidth',1.5)
end

subplot(1,3,3)
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
    nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.nbackmatches(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==3
        scatter(x,y,'o','Filled')
        high = [high; x' y'];
    end
    hold on
end
title('High NFC subjects')
ax = gca; fig = gcf;
ax.FontSize = 12;
fig.Color = 'w';
xlabel('N-back Matches')
ylabel('BDM value')
ylim([1 5.1])
[r,p] = corr(high(~isnan(high(:,2)),1),high(~isnan(high(:,2)),2),'Type','Spearman'); %high NFC has no relationship of n switches vs. value
if p<0.05
    y = r.*(unique(high(:,1)))+2.5;
    plot(unique(high(:,1)),y,'k','LineWidth',1.5)
end

%% Baseline executive function by NFC group
% Are high NFC people showing more baseline EF, as has been seen in
% previous work?
% Plot NFC group by practice round accuracy

split = tertile_split(data.NFC);

figure
subplot(1,2,1)
% accuracy by group
bars = [nanmean(tasks_overall(split==1,1)) nanmean(tasks_overall(split==2,1)) nanmean(tasks_overall(split==3,1)); ...
    nanmean(tasks_overall(split==1,2)) nanmean(tasks_overall(split==2,2)) nanmean(tasks_overall(split==3,2)); ...
    nanmean(tasks_overall(split==1,3)) nanmean(tasks_overall(split==2,3)) nanmean(tasks_overall(split==3,3)); ...
    nanmean(tasks_overall(split==1,4)) nanmean(tasks_overall(split==2,4)) nanmean(tasks_overall(split==3,4))];
E = [nanstd(tasks_overall(split==1,1)) nanstd(tasks_overall(split==2,1)) nanstd(tasks_overall(split==3,1)); ...
    nanstd(tasks_overall(split==1,2)) nanstd(tasks_overall(split==2,2)) nanstd(tasks_overall(split==3,2)); ...
    nanstd(tasks_overall(split==1,3)) nanstd(tasks_overall(split==2,3)) nanstd(tasks_overall(split==3,3)); ...
    nanstd(tasks_overall(split==1,4)) nanstd(tasks_overall(split==2,4)) nanstd(tasks_overall(split==3,4))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:12,bars(:,1),'FaceColor',NFCcolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:12,bars(:,2),'FaceColor',NFCcolors(2,:),'BarWidth',0.2)
bar(3:3:12,bars(:,3),'FaceColor',NFCcolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticks([2:3:12])
xticklabels(tasklabels)
title('Accuracy by NFC group')
legend({'Low NFC','Mid NFC','High NFC'})
xticklabels(tasklabels)
clean_fig();

% are high NFC-ers more accurate than everyone else on the 2-back?
[h,p] = ttest2(tasks_overall(split==3,3),tasks_overall(split==2,3));
[h,p] = ttest2(tasks_overall(split==3,3),tasks_overall(split==1,3));


% how does NFC interact with 'baseline EF' i.e. practice accuracy?
subplot(1,2,2)
pracn0acc = data.practiceacc(:,1);
prac3dacc = data.practiceacc(:,2);
pracn1acc = data.practiceacc(:,3); %first practice round for each task
pracn2acc = data.practiceacc(:,4);
bars = [nanmean(pracn0acc(split==1,:)) nanmean(pracn0acc(split==2,:)) nanmean(pracn0acc(split==3,:)); ...
    nanmean(pracn1acc(split==1,:)) nanmean(pracn1acc(split==2,:)) nanmean(pracn1acc(split==3,:)); ...
    nanmean(pracn2acc(split==1,:)) nanmean(pracn2acc(split==2,:)) nanmean(pracn2acc(split==3,:)); ...
    nanmean(prac3dacc(split==1,:)) nanmean(prac3dacc(split==2,:)) nanmean(prac3dacc(split==3,:))];
E = [nanstd(pracn0acc(split==1,:)) nanstd(pracn0acc(split==2,:)) nanstd(pracn0acc(split==3,:)); ...
    nanstd(pracn1acc(split==1,:)) nanstd(pracn1acc(split==2,:)) nanstd(pracn1acc(split==3,:)); ...
    nanstd(pracn2acc(split==1,:)) nanstd(pracn2acc(split==2,:)) nanstd(pracn2acc(split==3,:)); ...
    nanstd(prac3dacc(split==1,:)) nanstd(prac3dacc(split==2,:)) nanstd(prac3dacc(split==3,:))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
% bar(1:3:12,bars(:,1),'FaceColor',NFCcolors(1,:),'BarWidth',0.2)
% hold on
% bar(2:3:12,bars(:,2),'FaceColor',NFCcolors(2,:),'BarWidth',0.2)
% bar(3:3:12,bars(:,3),'FaceColor',NFCcolors(3,:),'BarWidth',0.2)
% xticks([2:3:12])
% bars = reshape(bars',1,numel(bars));
% E = reshape(E',1,numel(E));
% errorbar(bars,E,'*k')
% xticklabels(tasklabels)
% legend({'Low NFC','Mid NFC','High NFC'})
% title('Interaction of NFC group with practice accuracy')
scatter(data.practiceacc(:,4),data.NFC,[],NFCcolors(2,:),'Filled')
[r,p] = corr(data.practiceacc(~isnan(data.NFC),4),data.NFC(~isnan(data.NFC)));
xlabel('Practice Accuracy 2-back')
ylabel('NFC')
ax = gca; fig = gcf;
fig.Color = 'w'; ax.FontSize = 12;

%run ANOVA over PRACTICE accuracy by task and NFC group
maxlength = min([sum(split==1) sum(split==2) sum(split==3)]);
distances = [sum(split==1) sum(split==2) sum(split==3)]-maxlength;
trim = find(distances>0); matrix = [pracn0acc pracn1acc pracn2acc prac3dacc split];
idxes = [];
for t = 1:length(trim)
    temp = find(split==trim(t)); %trim larger NFC groups, need standard size for ANOVA
    temp = temp(randperm(length(temp)));
    idxes = [idxes; temp(1:distances(trim(t)))];
end
matrix(idxes,:) = [];
matrix = sortrows(matrix,5); matrix(isnan(matrix(:,5)),:) = [];
[~,~,stats] = anova2(matrix(:,1:4),maxlength,'off');
% no effect of NFC group on practice accuracy when all tasks taken into
% account

%2-way ANOVA on MAIN EXPERIMENT accuracy by task and NFC group
tasklist = [repmat(1,n,1); repmat(2,n,1); repmat(3,n,1); repmat(4,n,1)];
matrix = [reshape(tasks_overall,numel(tasks_overall),1) repmat(split,4,1) tasklist];
[~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{'NFC','task'},'Display','off');

[h,p] = ttest2(data.overall(split==1),data.overall(split==3));
[h,p] = ttest2(data.overall(split==2),data.overall(split==3));
[h,p] = ttest2(data.overall(split==1),data.overall(split==2));

% short description of what's in here:
% there is a clear effect of NFC group on overall accuracy, primarily
% driven by the difference between overall accuracy in mid NFC versus high
% NFC subjects
% there's also a very clear effect of task, where 2-back performance is
% much worse across the board

%% What's driving differences in 1- and 2-back ratings?
% Is it NFC?
% Is it just rating noise?

split = tertile_split(data.NFC);

n1 = data.task_displayed == 1;
n2 = data.task_displayed == 2;
ndetect = data.task_displayed==7;
figure
subplot(1,2,1)
for subj = 1:n
    diff(subj) = [nanmean(data.values(subj,n2(subj,:)))-nanmean(data.values(subj,n1(subj,:)))];
    scatter(diff(subj),data.NFC(subj))
    hold on
end
[r,p] = corr(diff(~isnan(data.NFC))',data.NFC(~isnan(data.NFC)),'Type','Spearman');
title([num2str(r) ' : relationship of NFC and mean diff between task values'])
ylabel('NFC')
xlabel('2-back minus 1-back fair wages')
clean_fig();

subplot(1,2,2)
for subj = 1:n
    scatter(diff(subj),nanstd(data.values(subj,:)))
    hold on
end
title('Relationship of STD of BDMs and mean diff between task values')
ylabel('STD of BDM fair wages (all from 1 subj)')
xlabel('2-back minus 1-back fair wages')
clean_fig();

%% What's up with NFC tertiles and difficulty ratings?
% Also, are NFC subjects better at paying attention to the task as a whole?

% look into self-rated competence i.e. how "easy" or "hard" they rated the
% task to be
figure
subplot(1,2,1)
means = []; E = [];
for group = 1:3
    relevant = data.diffrating(split==group,:);
    means(group) = nanmean(relevant); E(group) = nanstd(relevant)/sqrt(sum(split==group));
    bar(group,nanmean(relevant),'FaceColor',NFCcolors(group,:))
    hold on
    title(['Overall difficulty rating by NFC group'])
    ylabel('Mean rating (1 lowest, 5 highest)')
    xticks([1:3])
    xticklabels({'Low NFC','Mid NFC','High NFC'})
end
errorbar(means,E,'k','LineWidth',1.5)
fig = gcf; fig.Color = 'w';

subplot(1,2,2)
for group = 1:3
    relevant = data.blurs(split==group,:);
    means(group) = nanmean(relevant); E(group) = nanstd(relevant)/sqrt(sum(split==group));
    bar(group,nanmean(relevant),'FaceColor',NFCcolors(group,:))
    hold on
    scatter(group.*ones(1,length(relevant)),relevant,'k','Filled')
    title(['Mean # of times subject clicked away from task'])
    xticks([1:3])
    xticklabels({'Low NFC','Mid NFC','High NFC'})
end
errorbar(means,E,'k','LineWidth',1.5)
fig = gcf; fig.Color = 'w';

%% demographic data analysis
% does NFC change with sex or age?
% how about SAPS?
% how about reaction time and accuracy on the tasks?

figure
subplot(3,2,1)
bar([nanmean(data.NFC(data.sex==1)) nanmean(data.NFC(data.sex==2))])
hold on
xticklabels({'Men','Women'})
errorbar([nanmean(data.NFC(data.sex==1)) nanmean(data.NFC(data.sex==2))],[nanstd(data.NFC(data.sex==1)) nanstd(data.NFC(data.sex==2))]./[sqrt(sum(data.sex==1)) sqrt(sum(data.sex==2))],'k','LineWidth',1.5)
ylabel('Mean NFC')

subplot(3,2,2)
scatter(data.age,data.NFC,'Filled')
title('Age by NFC')
xlabel('Age')
ylabel('NFC')
[r,p] = corr(data.age,data.NFC);
if p<0.05
    lsline
end

subplot(3,2,3)
bar([nanmean(data.SAPS(data.sex==1)) nanmean(data.SAPS(data.sex==2))])
hold on
xticklabels({'Men','Women'})
errorbar([nanmean(data.SAPS(data.sex==1)) nanmean(data.SAPS(data.sex==2))],[nanstd(data.SAPS(data.sex==1)) nanstd(data.SAPS(data.sex==2))]./[sqrt(sum(data.sex==1)) sqrt(sum(data.sex==2))],'k','LineWidth',1.5)
ylabel('Mean SAPS')

subplot(3,2,4)
scatter(data.age,data.SAPS,'Filled')
title('Age by SAPS')
xlabel('Age')
ylabel('SAPS')
[r,p] = corr(data.age,data.SAPS);
if p<0.05
    lsline
end

subplot(3,2,5)
scatter(data.age,tasks_rts(:,2),[],taskcolors(1,:),'Filled')
hold on
scatter(data.age,tasks_rts(:,3),[],taskcolors(2,:),'Filled')
scatter(data.age,tasks_rts(:,4),[],taskcolors(3,:),'Filled')
legend({'1-back','2-back','3-detect'})
xlabel('Age')
ylabel('Mean RT')
title('Age by Mean RT')
[r,p] = corr(data.age,tasks_rts(:,3));
if p < 0.05
    lsline
end

subplot(3,2,6)
scatter(data.age,data.overall,'Filled')
xlabel('Age')
ylabel('Overall Acc')
title('Age by Overall Accuracy')
[r,p] = corr(data.age,data.overall);
if p<0.05
    lsline
end
fig = gcf; fig.Color = 'w';

%% Are subjects more accurate at the tasks when more BDM points are on the line?
% Does this differ across NFC tertiles?

% basic diff rating versus accuracy stuff - how well do people gauge their
% own abilities?

figure
subplot(2,2,1)
for task = 1:length(tasklabels)
    scatter(tasks_overall(:,task),data.taskratings(:,task),[],taskcolors(task,:),'Filled')
    hold on
    prune = [tasks_overall(:,task) data.taskratings(:,task)];
    prune(isnan(tasks_overall(:,task))|isnan(data.taskratings(:,task)),:) = [];
    [r,p] = corr(prune(:,1),prune(:,2));
    disp(['relationship between ' tasklabels{task} ' overall acc and diff rating r = ' num2str(r) ' p = ' num2str(p)])
end
title('Difficulty rating by overall accuracy')
ylabel('Rating')
xlabel('Overall Accuracy')
legend(tasklabels)
fig = gcf; fig.Color = 'w';

subplot(2,2,2)
columns = {'detectacc','n1acc','n2acc','ndetectacc'};
for task = 1:size(columns,2)
    eval(['column = data.' columns{task} ';'])
    for subj = 1:n
        endacc = column(subj,~isnan(column(subj,:)));
        if length(endacc)>0
            endacc = endacc(end);
        else
            endacc = NaN;
        end
        endaccs(subj,task) = endacc;
    end
end
for task = 1:length(tasklabels)
    scatter(endaccs(:,task),data.taskratings(:,task),[],taskcolors(task,:),'Filled')
    hold on
    prune = [endaccs(:,task) data.taskratings(:,task)];
    prune(isnan(endaccs(:,task))|isnan(data.taskratings(:,task)),:) = [];
    [r,p] = corr(prune(:,1),prune(:,2));
    if p<0.05
        disp(['relationship between ' tasklabels{task} ' final acc and diff rating r = ' num2str(r) ' p = ' num2str(p)])
    end
end
title('Difficulty rating by final accuracy')
ylabel('Rating')
xlabel('Final Accuracy')
legend(tasklabels)
fig = gcf; fig.Color = 'w';

subplot(2,2,3)
%plot difficulty rating by practice attempts on 2-back
backattempts = data.practiceacc(:,4:end); backattempts = sum(~isnan(backattempts),2);
scatter(backattempts, data.taskratings(:,3),[],taskcolors(3,:),'Filled')
title('# 2-back practices vs. diff rating')
xlabel('Number of attempts at the 2-back in practice')
ylabel('Difficulty rating on 2-back')
backattempts = backattempts(~isnan(data.taskratings(:,1))); taskratings = data.taskratings(~isnan(data.taskratings(:,1)),:);
[r,p] = corr(backattempts,taskratings(:,3));
if p < 0.05
    lsline
end

subplot(2,2,1)
for subj = 1:n
    points = NaN(1,default_length);
    idx = data.values(subj,:)<data.offers(subj,:);
    points(:,idx) = data.offers(subj,idx);
    idx = data.values(subj,:)>data.offers(subj,:);
    points(:,idx) = ones(1,sum(idx));
    BDMs = points;
    perf = data.perf(subj,:);
    scatter(perf, BDMs, 'Filled')
    hold on
end
title('Points on Offer by Performance')
xlabel('Accuracy')
ylabel('Points to Win')

split = tertile_split(data.NFC);
NFCdata{1} = []; NFCdata{2} = []; NFCdata{3} = [];
for subj = 1:n
    points = NaN(1,default_length);
    idx = data.values(subj,:)<data.offers(subj,:);
    points(:,idx) = data.offers(subj,idx);
    idx = data.values(subj,:)>data.offers(subj,:);
    points(:,idx) = ones(1,sum(idx));
    BDMs = points;
    perf = data.perf(subj,:);
    group = split(subj);
    if ~isnan(group)
        NFCdata{group} = [NFCdata{group}; perf' BDMs'];
        %scatter(perf, BDMs, [], NFCcolors(group,:),'Filled')
    end
end
labels = {'Low NFC','Mid NFC','High NFC'};
for group = 1:max(split)
    subplot(2,2,1+group)
    toplot = count_entries(NFCdata{group});
    scatter(toplot(:,1),toplot(:,2),[],NFCcolors(group,:),'Filled')
    [r,p] = corr(toplot(:,1),toplot(:,2));
    if p < 0.05
        lsline
    end
    title('Points on Offer by Performance')
    legend(labels{group})
    ylabel('Points')
    xlabel('Accuracy')
end

%% modeling results (parameters) versus behavioral results (NFC and SAPS)
% this is a little more extensive than what I have going on in
% paper_graphs_and_stats.m

% more params, more groups plotted (as in, people also best fit by 2nd and
% 3rd winning models, not just people best fit by 1st winning model)

modelfits = load('modeling/HBI/HBI_modelStruct_2023.mat');
no_fits = load('simdata/toanalyze.mat','trim'); no_fits = no_fits.trim;
% remove subjects who weren't fit to models
% grab best model
best_model = modelfits.best_model;
responsibility = modelfits.best_model.overallfit.responsibility;
models_at_play = find(sum(responsibility>=0.01,1)>0);
assignments(responsibility(:,models_at_play(1))>responsibility(:,models_at_play(2))) = models_at_play(1);
assignments(responsibility(:,models_at_play(2))>responsibility(:,models_at_play(1))) = models_at_play(2); %split into best model groups

modelnames = best_model.overallfit.fitmodels;
nparams = best_model.nparams;
params = applyTrans_parameters(best_model,best_model.lowparams); paramnames = best_model.paramnames;
measures = [data.NFC data.SAPS]; measures(no_fits,:) = [];
%remove subject NFC and SAPS who don't have model results

names = {'NFC','SAPS'};
for measure = 1:size(measures,2)
    split = tertile_split(measures(:,measure));
    column = measures(:,measure);
    eval(['colors =' names{measure} 'colors;']);
    figure
    for p = 1:nparams
        subplot(3,2,p)
        scatter(column,params(:,p),[],colors(2,:),'Filled')
        [r,pval] = corr(column(~isnan(column)),params(~isnan(column),p));
        if pval < 0.05
            lsline
        end
        xlabel(names{measure})
        ylabel(paramnames(p))
        disp([names{measure} ' vs ' paramnames(p) ' r = ' num2str(r) ', p = ' num2str(pval)])
    end
    fig = gcf; fig.Color = 'w';
    
    
    figure
    for p = 1:nparams
        subplot(3,2,p)
        for group = 1:3
            relevant = params(split==group,p);
            bar(group,nanmean(relevant),'FaceColor',colors(group,:))
            hold on
        end
        errorbar([nanmean(params(split==1,p)) nanmean(params(split==2,p)) nanmean(params(split==3,p))], ...
            [nanstd(params(split==1,p)) nanstd(params(split==2,p)) nanstd(params(split==3,p))]./sqrt([sum(split==1) sum(split==2) sum(split==3)]),'*k')
        title(paramnames(p))
        xticks([1:group])
        xticklabels({['Low ' names{measure}],['Mid ' names{measure}],['High ' names{measure}]})
        % run group-level t-tests
        [h,pval] = ttest2(params(split==1,p),params(split==2,p));
        disp([paramnames{p} ': ' num2str(pval) ' difference between low and mid ' names{measure}])
        [h,pval] = ttest2(params(split==3,p),params(split==2,p));
        disp([paramnames{p} ': ' num2str(pval) ' difference between high and mid ' names{measure}])
        [h,pval] = ttest2(params(split==1,p),params(split==3,p));
        disp([paramnames{p} ': ' num2str(pval) ' difference between low and high ' names{measure}])
    end
    fig = gcf; fig.Color = 'w';
    
    figure
    bar([nanmean(column(assignments==models_at_play(1))); nanmean(column(assignments==models_at_play(2)))],'FaceColor','w')
    hold on
    errorbar([nanmean(column(assignments==models_at_play(1))); nanmean(column(assignments==models_at_play(2)))],[nanstd(column(assignments==models_at_play(1))); nanstd(column(assignments==models_at_play(2)))]./sqrt(n),'*k')
    scatter(ones(sum(assignments==models_at_play(1)),1),column(assignments==models_at_play(1)),[],colors(2,:),'Filled')
    scatter(2*ones(sum(assignments==models_at_play(2)),1),column(assignments==models_at_play(2)),[],colors(2,:),'Filled')
    ylabel(names{measure})
    xticklabels(modelnames(models_at_play))
    xtickangle(45)
    fig = gcf; fig.Color = 'w';
    title(['Mean ' names{measure} ' by model class'])
end

% Assign subjects to their model, plot parameters that way
all_params = best_model.overallfit.parameters;
for measure = 1:2
    if measure == 1
        split = tertile_split(data.NFC); colors = NFCcolors;
        % Also, run multiple linear regressions
        X = [ones(n,1) data.NFC data.NFC.^2];
    elseif measure == 2
        split = tertile_split(data.SAPS); colors = SAPScolors;
        X = [ones(n,1) data.SAPS data.SAPS.^2];
    end
    invalid = sum(isnan(X),2)>0; X(invalid,:) = [];
    
    % initialize plotting variables
    plotted = []; count = 0;
    figure; fig = gcf; fig.Color = 'w';
    
    for m = 1:length(models_at_play)
        modelnum = models_at_play(m);
        nparams = size(all_params{modelnum},2);
        name = best_model.overallfit.fitmodels{modelnum};
        paramnames = strsplit(strrep(name,'c',' cost'),'-');
        paramnames = strrep(paramnames,'u ','update '); 
        paramnames = strrep(paramnames,'main ','maintenance ');
        % add clearer labels to param plots
        costs = find(contains(paramnames,'cost'));
        values = all_params{modelnum}(:,costs);
        for c = 1:size(costs,2)
            count = count+1;
            
            % run stats on measurement groups
            % [p,~,stats] = anova1(values(:,1),split,'off');
            %if p < 0.05; [compared] = multcompare(stats,'display','on'); end
            % run group-based & continuous stats
            [h,pval] = ttest2(values(split==1,c),values(split==2,c));
            disp([paramnames{costs(c)} ': ' num2str(pval) ' difference between low and mid ' names{measure}])
            [h,pval] = ttest2(values(split==3,c),values(split==2,c));
            disp([paramnames{costs(c)} ': ' num2str(pval) ' difference between high and mid ' names{measure}])
            [h,pval] = ttest2(values(split==1,c),values(split==3,c));
            disp([paramnames{costs(c)} ': ' num2str(pval) ' difference between low and high ' names{measure}])
            
            Y = nanmean(values,2); Y(invalid,:) = [];
            [betas,BINV,~,~,stats] = regress(Y,X);
            % get betas for quadratic term
            predicted = X*betas;
            distance = predicted-Y; MSE = distance'*distance; %squared distance
            disp([names{measure} ' betas on ' paramnames{costs(c)} ', linear & quadratic: ' num2str(betas')])
            
            % plot effect of measure group on parameter values
            subplot(8,2,count)
            for group = 1:3
                relevant = values(split==group,c);
                bar(group,nanmean(relevant),'FaceColor',colors(group,:))
                hold on
            end
            errorbar([nanmean(values(split==1,c)) nanmean(values(split==2,c)) nanmean(values(split==3,c))], ...
                [nanstd(values(split==1,c)) nanstd(values(split==2,c)) nanstd(values(split==3,c))]./sqrt([sum(split==1) sum(split==2) sum(split==3)]),'*k')
            title(paramnames{costs(c)})
            xticks([1:group])
            xticklabels({['Low ' names{measure}],['Mid ' names{measure}],['High ' names{measure}]})
            
        end %of cycling over each cost
        
        plotted = [plotted paramnames(costs)];
    end %of cycling over each model
    
end %of cycling over SAPS & NFC

%% RELATE MORE MODEL-AGNOSTIC FINDINGS TO MORE MODEL-BASED FINDINGS 

% is error cost related to perfectionism?

%1) Standard correlation of error cost parameter magnitude and
%perfectionism scores (SAPS)

% where lure cost parameter in output?
modelstofit = best_model.overallfit.fitmodels;
model = coc_createModels(modelstofit{5});
column = find(contains(model.paramnames,'fac'));
cbm = best_model.cbm;
[~,assignments] = max(cbm.output.responsibility,[],2);

% interference costs from second-best model (number 3)
fa_costs = cbm.output.parameters{5}(:,column);

[r,p] = corr(fa_costs(~isnan(data.SAPS)),data.SAPS(~isnan(data.SAPS)));

%2) Tertile split, then examine SAPS tertiles' FA costs
% doesn't really make sense because... FA costs only explain 8 subjects'
% data, so not the best way to look at stuff, but I'm curious.
count = 0;
split = tertile_split(data.SAPS);
low = fa_costs(split==1); mid = fa_costs(split==2); high = fa_costs(split==3);

[h,pval1] = ttest2(low,mid);
[h,pval2] = ttest2(high,mid);
[h,pval3] = ttest2(low,high);
ps = [NaN pval1 pval3; pval1 NaN pval2; pval3 pval2 NaN];

X = [ones(length(data.SAPS),1) data.SAPS]; invalid = isnan(data.SAPS);
X(invalid,:) = [];
Y = fa_costs; Y(invalid,:) = [];
[betas,BINV,~,~,stats] = regress(Y,X);
% get betas for quadratic term
predicted = X*betas;
distance = predicted-Y; MSE = distance'*distance; %squared distance
if stats(3)<0.05
    disp(['SAPS significant betas on FA costs, linear & quadratic: ' num2str(betas')])
    disp(['p value = ' num2str(stats(3))])
end

figure;
% plot effect of SAPS group on FA cost values
superbar([nanmean(low) nanmean(mid) nanmean(high)], ...
    'E',[nanstd(low) nanstd(mid) nanstd(high)]./sqrt(n), ...
    'P',ps,'BarFaceColor',SAPScolors,'PStarShowNS',false,'PStarBackgroundColor','None');
ylabel('False alarm costs')
%title(model_labels{modelnum})
xticks([1:3])
xtickangle(30)
xticklabels({['Low'],['Mid'],['High']})
xlabel(['SAPS group'])
ylabel('Mean false alarm cost')
clean_fig();

% 3) are subjects in FA cost model more perfectionist than other subjects?
FA_subjects = assignments==5 | assignments==16; others = assignments~=5;
figure;
bar([nanmean(data.SAPS(FA_subjects)) nanmean(data.SAPS(others))],'FaceColor',SAPScolors(2,:))
hold on
errorbar([nanmean(data.SAPS(FA_subjects)) nanmean(data.SAPS(others))],...
    [nanstd(data.SAPS(FA_subjects)) nanstd(data.SAPS(others))]./[sqrt(sum(FA_subjects)) sqrt(sum(others))],...
    'k*')
xticklabels({'Subjs w/ FA costs','Subjs w/o'})
ylabel('Mean SAPS score')
clean_fig();


% % Another model-agnostic to model-dependent finding:
% are subjects with more than 1 cost more effort-approaching than others?
% no, but subject numbers are way off here, anyway
more_than_one_cost = assignments>6;
figure
bar([nanmean(data.NFC(more_than_one_cost)) nanmean(data.NFC(~more_than_one_cost))], ...
    'FaceColor',NFCcolors(2,:))
hold on
errorbar([nanmean(data.NFC(more_than_one_cost)) nanmean(data.NFC(~more_than_one_cost))], ...
    [nanstd(data.NFC(more_than_one_cost))./sqrt(sum(more_than_one_cost)) nanstd(data.NFC(~more_than_one_cost))./sqrt(sum(~more_than_one_cost))], ...
    '*k','LineWidth',2)
ylabel('NFC score')
xlabel('Group')
xticklabels({'Subjects w more than 1 cost','Subjects w 1 cost in best model'})
xtickangle(45)
clean_fig()
[h,p] = ttest2(data.NFC(more_than_one_cost),data.NFC(~more_than_one_cost));

% % Lastly, subjects who randomly respond to NFC questionnaire, as well as
% respond randomly on the BDM, will have 1) mid NFC score and 2) poor model
% fit as exhibited by high sigma parameter. Is that intuition correct?

%  Some subjects are poorly described even by their best model
% Identify outliers in sigma to identify these subjects
% They may be influencing other analyses in an outsized way (e.g. NFC
% analyses)

for sii = 1:n
    
    model = assignments(sii);
    model_structure = coc_createModels(modelstofit{model});
    sigma_index = contains(model_structure.paramnames,'epsilon');
    % yes, confusingly, I misnamed the std (sigma) parameter at first, and just
    % left it that way to avoid any strange bugs that could arise from
    % trying to re-name everything
    sigma_sii = cbm.output.parameters{model}(sii,sigma_index);
    sigmas_from_best_models(sii) = exp(sigma_sii);
    
end

%visually identify outliers
figure
scatter(1:n,sort(sigmas_from_best_models))
ylabel('\sigma value from best model')
xlabel('Subject # (sorted)')
clean_fig();

% cutoff, visually, appears to be sigma values of/above 1 
bad_fits = find(sigmas_from_best_models>1);
good_fits = 1:n; good_fits(bad_fits) = [];

figure
errorbar([nanmean(data.NFC(good_fits)) nanmean(data.NFC(bad_fits))], ...
    [nanstd(data.NFC(good_fits)) nanstd(data.NFC(bad_fits))]./[sqrt(length(good_fits)) sqrt(length(bad_fits))], ...
    'k','LineWidth',2)
set(gcf,'color','w')
ylabel('Mean NFC score')
xlabel('Group')
xticks([1 2]); xlim([0.75 2.25])
xticklabels({'Good model fits','Bad fits'})

