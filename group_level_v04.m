%% Analysis of cost of control task data 
% task written in JsPsych, data printed in CSV form
% Analysis script for version 3 of the task, the n-back version
%% Process raw jspsych text file data
clear all

% files = string(ls('./data')); files(1:2) = []; %trim directory entries
% files = files(contains(files,'.mat'));
% files = {'11.02.2020.mat','12.02.2020_trial.mat','13.02.2020_trial.mat'};
% version 1 ^
% files = {'18.02.2020.mat','19.02.2020.mat','24.02.2020.mat','25.02.2020.mat'};
% version 2
%files = {'24.03.2020.mat','26.03.2020.mat','31.03.2020.mat','07.04.2020.mat'}; %% %wrong # of variables here %'06.04.2020.mat'
% version 3
files = {'29.04.2020.mat','07.05.2020.mat','12.05.2020.mat','13.05.2020.mat','18.05.2020.mat'};
% version 4
[data,excluded] = loadcostdata(files); %load data in
n = height(data);

tasks = [categorical(cellstr('detection'));categorical(cellstr('n1')); categorical(cellstr('ndetection')); categorical(cellstr('n2'))];
tasklabels = {'1-detect','1-back','3-detect','2-back'};
tasknumbers = [0,1,7,2];
default_length = 32;

inattentive = data.perf<60;

% save important variables for later
for subj = 1:n
    for task = 1:length(tasks)
        list = find(data.task_progression(subj,:)==tasks(task));
        tasks_overall(subj,task) = nanmean(data.perf(subj,list));
    end
end
tasks_rts = [nanmean(data.detectrts,2) nanmean(data.n1rts,2) nanmean(data.n2rts,2) nanmean(data.ndetectrts,2)];

for task = 2
    data = run_stats(data,task);
end %run stats separately on 1 and 2 back

data.taskfreqs = [sum(data.task_progression==tasks(1),2) sum(data.task_progression==tasks(2),2) sum(data.task_progression==tasks(3),2) sum(data.task_progression==tasks(4),2)];
meanBDM = nanmean(data.values(:,2:end),2);

%display variables
subjcolors = rand(n);subjcolors(:,4:end) = []; %delete unnecessary columns
taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
%style variables
NFCcolors = [.5 0 .5; .75 0 .75; .95 0 .95];
SAPScolors = [0 .25 .25; 0 .5 .5; 0 .75 .75];

% list = data.subjnum;
% save('fullsubjnumbers.mat','list');

%% Check data quality - how long each task takes, how many times they've been completed

figure
subplot(1,2,1)
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

subplot(1,2,2)
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
%% accuracy, rt, important BDM plots

%Mean accuracy by task, mean RTs, make sure "harder" tasks are actually
%harder
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

% subplot(2,2,4)
% bar(nanmean(data.task_blurs))
% hold on
% errorbar(nanmean(data.task_blurs),nanstd(data.task_blurs)/sqrt(n),'LineWidth',1.5)
% title('Blur events by task')
% xticklabels(tasklabels)
% this is nothing, no one tunes out during tasks I guess?

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

% figure
% %what does this bdm x iter thing look like subj by subj?
% for task = 1:2
%     subplot(1,2,task)
%     eval(['toplot = n' num2str(task) 'subjvalue;'])
%     for subj = 1:n
%         plot(1:default_length+1,toplot(subj,:))
%         hold on
%     end
%     title(['BDM by iteration: ' num2str(task) '-back'])
% end

%plot task avoidance, well, sort of
% for i = 1:n
%     plot(1:3,data.taskfreqs(i,:),'LineWidth',1)
%     hold on
% end
% ylabel('# of times task completed')
% title('Task avoidance by subject')
% xlabel('Task')
% xticks([1:3])
% xticklabels(tasklabels)
% ax = gca; fig = gcf;
% ax.FontSize = 12; fig.Color = 'w';

% plot task performance by block, any learning or decay?
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
[h,p] = ttest(100.*ones(30,1),data.detectacc(:,i));
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

%% Plot first rating versus total iterations

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
ignore_flag = true;
if ~ignore_flag
figure
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

%plot BDM request by previous offer
figure
subplot(1,2,1)
for row = 1:n
    y = data.offers(row,:)';
    x = [NaN data.values(row,2:end)]';
    matrix = sortrows([y,x],1);
    plot(matrix(:,1),matrix(:,2),'o')
    hold on
end
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title('BDM value by prev. computer offer')
ylim([1 5.1])
ylabel('BDM points')
xlabel('Last offer')

%plot late responses/changed responses by task
subplot(1,2,2)
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

figure
rts = data.incorrectrts;
% just kidding, there aren't any incorrect rts yet...
for task = 1:length(tasklabels)
    subplot(2,2,task)
    histogram(rts{task})
    title(['Dist. of Incorrect RTs for ' tasklabels(task)])
    fig = gcf; fig.Color = 'w';
end
end %of plot_flag if

%% Understand n-back matches' effect on costs

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
        completed = NaN(1,default_length);completed(data.task_progression(subj,:)==tasks(2)) = 2;
        completed(data.task_progression(subj,:)==tasks(3)) = 3;completed(data.task_progression(subj,:)==tasks(4)) = 8;
        completed(:,end) = []; completed = [NaN completed];%delete last column
        same = completed==displayed; 
        different = ~isnan(completed)&completed~=displayed; %different but not detection
        same_long = [same_long; matches(same&idx)' BDMs(same&idx)'];
        diff_long = [diff_long; matches(different&idx)' BDMs(different&idx)'];
    end
    toplot = countEntries(diff_long);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'c','Filled')
    hold on
    toplot = countEntries(same_long);
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

for task = 1:3
    eval(['columns = [n' num2str(task) 'effect(:,2) n' num2str(task) 'effect(:,5)];']);
    trim = isnan(columns(:,1))|isnan(columns(:,2));
    columns(trim,:) = [];
    [r,p] = corr(columns(:,1),columns(:,2),'Type','Spearman');
    disp(['effect of prev matches on BDM, ' tasklabels{task+1} ': r = ' num2str(r) ' p = ' num2str(p)])
end

% for task = 1:2
%     subplot(3,3,task+6)
%     for subj = 1:n
%         eval(['y = data.n' num2str(task) 'matcheffect(subj,:);'])
%         plot(1:3,y,'LineWidth',1.5)
%         hold on
%     end
%     title(['Match "effect" on BDM from last ' num2str(task) '-back'])
%     eval(['y = data.n' num2str(task) 'matcheffect;'])
%     errorbar(nanmean(y),nanstd(y)/sqrt(n),'k','LineWidth',1.5)
%     xticklabels({'3','4','5'})
%     xlabel('Matches in Last Round of Task')
%     legend({'Subj 1','Subj 2','...','Subj n'})
% end

figure
subplot(2,2,1)
toplot = countEntries([n2effect(:,4) n2effect(:,6)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
hold on
toplot = countEntries([n1effect(:,4) n1effect(:,6)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
legend({'2-back','1-back'})
ylabel('Accuracy')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';

subplot(2,2,2)
toplot = countEntries([n2effect(:,4:5)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
hold on
toplot = countEntries([n1effect(:,4:5)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
legend({'2-back','1-back'})
ylabel('BDM request')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';

subplot(2,2,3)
%effect structures are 1. subj 2. matches 3. misses 4. delay 5. BDM 6. perf
toplot = countEntries([n2effect(:,3:4)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'b','Filled')
hold on
toplot = countEntries([n1effect(:,3:4)]);
scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1).*5,'r','Filled')
legend({'2-back','1-back'})
ylabel('Missed matches')
title('Effect of delay since last time task completed')
ax = gca; ax.FontSize = 12; fig = gcf; fig.Color = 'w';

figure
subplot(2,2,1)
for subj = 1:n
   effect = data.matcheffect(subj,:); 
   plot(1:3,effect,'LineWidth',1.5)
   hold on
end
errorbar(1:3,nanmean(data.matcheffect),nanstd(data.matcheffect)/sqrt(n),'k','LineWidth',2)
legend({'Subj 1','Subj 2','Subj ...','Subj N'},'Location','Best')
title('"Effect" of N-Back Matches for Each Subject, irrespective of task')
xticks([1:3])
xticklabels({'3 matches','4 matches','5 matches'})
xlabel('Correctly answered n-back matches')
ylabel('mean BDM value following that number of matches')
fig = gcf; fig.Color = 'w';
xlim([0.5 3.5])

subplot(2,2,2)
for subj = 1:n
   effect = data.missedeffect(subj,:); 
   plot(1:3,effect,'LineWidth',1.5)
   hold on
end
errorbar(1:3,nanmean(data.missedeffect),nanstd(data.missedeffect)/sqrt(n),'k','LineWidth',2)
legend({'Subj 1','Subj 2','Subj ...','Subj N'},'Location','Best')
title('"Effect" of Missed N-Back Matches for Each Subject, irrespective of task')
xticks([1:3])
xticklabels({'2 missed','3 missed','4 missed'})
xlabel('Missed n-back matches')
ylabel('mean BDM value')
fig = gcf; fig.Color = 'w';
xlim([0.5 3.5])

subplot(2,2,4)
for subj = 1:n
    nback = find((data.task_progression(subj,:)==tasks(2))|(data.task_progression(subj,:)==tasks(3)));
    next = (nback)+1; nback(next==default_length+1) = []; next(next==default_length+1) = []; 
    BDMs = data.values(subj,next);
    missed = data.nbackmisses(subj,nback);
    scatter(missed,BDMs,'Filled')
    hold on
end
xlabel('Number of Missed N-Back Matches')
ylabel('BDM value')
title('"Effect" of Missed N-Back Matches for Each Subject')

subplot(2,2,3)
for subj = 1:n
    nback = find((data.task_progression(subj,:)==tasks(2))|(data.task_progression(subj,:)==tasks(3)));
    next = (nback)+1; nback(next==default_length+1) = []; next(next==default_length+1) = []; 
    BDMs = data.values(subj,next);
    missed = data.intendednbackmatches(subj,nback);
    scatter(missed,BDMs,'Filled')
    hold on
end
xlabel('Number of N-Back Matches, missed or correct')
ylabel('BDM value')

% look at survivor analysis
% if they made one error on either task, how likely were they to just give
% up?

%% deal with individual differences in perfectionism
%style variables
SAPScolors = [0 .25 .25; 0 .5 .5; 0 .75 .75];

figure
subplot(1,3,1)
histogram(data.NFC)
title('Hist of NFC scores (normalized)')
subplot(1,3,2)
histogram(data.SAPS)
title('Hist of SAPS scores (normalized)')
meanBDM = nanmean(data.values(:,2:end),2);
subplot(1,3,3)
scatter(data.NFC,data.SAPS,'Filled')
lsline
title('NFC by SAPS score')
xlabel('NFC score')
ylabel('SAPS score')
[r,p] = corr(data.NFC(~isnan(data.NFC)),data.SAPS(~isnan(data.NFC)),'Type','Spearman');
if p<0.05
    hold on
    lsline
end
fig = gcf; fig.Color = 'w';
figure
%tertile split SAPS scores
split = tertileSplit(data.SAPS);
%plot accuracy by SAPS score

subplot(1,3,1)
%bar([nanmean(data.perf(split==1)) nanmean(data.perf(split==2)) nanmean(data.perf(split==3))])
% bar(1,[nanmean(data.overall(split==1))],'c')
% hold on
% bar(2,[nanmean(data.overall(split==2))],'b')
% bar(3,[nanmean(data.overall(split==3))],'r')
% errorbar([nanmean(data.overall(split==1)) nanmean(data.overall(split==2)) nanmean(data.overall(split==3))],[nanstd(data.overall(split==1)) nanstd(data.overall(split==2)) nanstd(data.overall(split==3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))],'*k','LineWidth',1.5)
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
fig = gcf; fig.Color = 'w';

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
    fig = gcf; fig.Color = 'w';
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
    fig = gcf; fig.Color = 'w';
end
end

figure
%slower performance, higher perfectionism?
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

figure
subplot(1,3,1)
pracn0acc = data.practiceacc(:,1);
prac3dacc = data.practiceacc(:,2);
pracn1acc = data.practiceacc(:,3);
pracn2acc = data.practiceacc(:,4); %first practice round for each task
bars = [nanmean(pracn0acc(split==1,:)) nanmean(pracn0acc(split==2,:)) nanmean(pracn0acc(split==3,:)); ...
    nanmean(pracn1acc(split==1,:)) nanmean(pracn1acc(split==2,:)) nanmean(pracn1acc(split==3,:)); ...
    nanmean(pracn2acc(split==1,:)) nanmean(pracn2acc(split==2,:)) nanmean(pracn2acc(split==3,:)); ...
    nanmean(prac3dacc(split==1,:)) nanmean(prac3dacc(split==2,:)) nanmean(prac3dacc(split==3,:))];
E = [nanstd(pracn0acc(split==1,:)) nanstd(pracn0acc(split==2,:)) nanstd(pracn0acc(split==3,:)); ...
    nanstd(pracn1acc(split==1,:)) nanstd(pracn1acc(split==2,:)) nanstd(pracn1acc(split==3,:)); ...
    nanstd(pracn2acc(split==1,:)) nanstd(pracn2acc(split==2,:)) nanstd(pracn2acc(split==3,:)); ...
    nanstd(prac3dacc(split==1,:)) nanstd(prac3dacc(split==2,:)) nanstd(prac3dacc(split==3,:))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
% bar(1:3:12,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
% hold on
% bar(2:3:12,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
% bar(3:3:12,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
% bars = reshape(bars',1,numel(bars));
% E = reshape(E',1,numel(E));
% errorbar(bars,E,'*k')
% xticks([2:3:12])
% xticklabels(tasklabels)
% legend({'Low SAPS','Mid SAPS','High SAPS'})
% title('Interaction of SAPS group with practice accuracy')
scatter(data.practiceacc(:,4),data.SAPS,[],SAPScolors(2,:),'Filled')
[r,p] = corr(data.practiceacc(~isnan(data.SAPS),4),data.SAPS(~isnan(data.SAPS)))
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

%% Tertile split NFC scores
% individual differences in NFC and cognition
split = tertileSplit(data.NFC);
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

n1effect = []; n2effect = []; n3effect = []; % keep track of numbers pulled out here for stats later
figure
for task = 1:(length(tasklabels)-1) %cycle through what is being displayed, display by what happened before
    subplot(2,2,task)
    for subj = 1:n %go subject by subject
        group = split(subj);
        displayed = data.task_displayed(subj,:);
        idx = find(displayed == tasknumbers(task+1)); 
        matches = data.nbackmatches(subj,:);
        if task == 3
            matches = data.ndetectmatches(subj,:);
        end
        perf = data.perf(subj,:);
        BDMs = data.values(subj,:);
        completed = find(data.task_progression(subj,:)==tasks(task+1));
        for trial = 1:length(idx)
            now = idx(trial);
            if sum(completed<now)>0 %they've done the last once before
                last = completed(completed<now); last = last(end);
                delay = now-last;
                if ~isnan(group)
                    %color = colors{group};
                    %scatter(matches(last),BDMs(now),color,'Filled')
                    eval(['n' num2str(task) 'effect = [n' num2str(task) 'effect; subj matches(last) delay BDMs(now) group];'])
                end
            end
        end
    end
    eval(['low = n' num2str(task) 'effect(n' num2str(task) 'effect(:,5)==1,:); mid = n' num2str(task) 'effect(n' num2str(task) 'effect(:,5)==2,:); high = n' num2str(task) 'effect(n' num2str(task) 'effect(:,5)==3,:);']) 
    toplot = countEntries([mid(:,2),mid(:,4)]);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1)*10,NFCcolors(2,:),'Filled');
    hold on
    toplot = countEntries([low(:,2),low(:,4)]);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1)*10,NFCcolors(1,:),'Filled');
    toplot = countEntries([high(:,2),high(:,4)]);
    scatter(toplot(:,1),toplot(:,2),(toplot(:,3)+1)*10,NFCcolors(3,:),'Filled');
    ylabel('BDM points')
    title(['BDM for ' tasklabels(task+1)])
    legend({'Mid NFC','Low NFC','High NFC'})
    xlabel(['n-back matches in last block'])
    fig = gcf;
    fig.Color = 'w';
    ylim([1 5.1])
end

low = n1effect(n1effect(:,5)==1,:); mid = n1effect(n1effect(:,5)==2,:); high = n1effect(n1effect(:,5)==3,:); 
low(isnan(low(:,4)),:) = []; mid(isnan(mid(:,4)),:) = []; high(isnan(high(:,4)),:) = []; 

[r,p] = corr(low(:,2),low(:,4));
% disp('low NFC spearman corr of matches and BDM, 1-back')

[r,p] = corr(mid(:,2),mid(:,4));
% disp('mid NFC spearman corr of matches and BDM, 1-back')

[r,p] = corr(high(:,2),high(:,4));
% disp('high NFC spearman corr of matches and BDM, 1-back')

low = n2effect(n2effect(:,5)==1,:); mid = n2effect(n2effect(:,5)==2,:); high = n2effect(n2effect(:,5)==3,:); 
low(isnan(low(:,4)),:) = []; mid(isnan(mid(:,4)),:) = []; high(isnan(high(:,4)),:) = []; 

[r,p] = corr(low(:,2),low(:,4));
% disp('low NFC spearman corr of matches and BDM, 2-back')

[r,p] = corr(mid(:,2),mid(:,4));
% disp('mid NFC spearman corr of matches and BDM, 2-back')

[r,p] = corr(high(:,2),high(:,4));
% disp('high NFC spearman corr of matches and BDM, 2-back')

subplot(2,2,4)
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

% are high NFC-ers more accurate than everyone else on the 2-back?
[h,p] = ttest2(tasks_overall(split==3,3),tasks_overall(split==2,3));
[h,p] = ttest2(tasks_overall(split==3,3),tasks_overall(split==1,3));


% how does NFC interact with 'baseline EF' i.e. practice accuracy?
figure
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

%is high NFC-er's baseline EF better/higher? maybe
[~,~,stats] = anova1(pracn2acc,split);

maxlength = min([sum(split==1) sum(split==2) sum(split==3)]);
toolong = [sum(split==1) sum(split==2) sum(split==3)]>maxlength;
idxes = find(split==find(toolong));
split(idxes(end)) = NaN; %trim it so it's 9, 9, and 9
matrix = [pracn0acc pracn1acc pracn2acc prac3dacc split];
matrix = sortrows(matrix,5); matrix(isnan(matrix(:,5)),:) = []; matrix(:,5) = [];
[~,~,stats] = anova2(matrix,length(unique(split(~isnan(split)))));

%2-way anova on overall accuracy by task and NFC group
tasklist = [repmat(1,n,1); repmat(2,n,1); repmat(3,n,1); repmat(4,n,1)];
matrix = [reshape(tasks_overall,numel(tasks_overall),1) repmat(split,4,1) tasklist];
[~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{'NFC','task'});

[h,p] = ttest2(data.overall(split==1),data.overall(split==3));
[h,p] = ttest2(data.overall(split==2),data.overall(split==3));
[h,p] = ttest2(data.overall(split==1),data.overall(split==2));

%% try to understand random effects by plotting one plot/subject, BDM vs.
%nswitches, NFC as title %%%%%%%
split = tertileSplit(data.NFC);

% figure
% for row = 1:n
%     subplot(6,7,row)
%     for task = 2:3
%         nback = find(data.task_progression(row,:)==tasks(task));
%         nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
%         next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
%         x = data.nbackmatches(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
%         y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
%         if task == 2
%             color = 'ob';
%         elseif task == 3
%             color = 'or';
%         end
%         scatter(x,y,color,'Filled')
%         hold on
%         title(num2str(data.NFC(row)))
%         ax = gca; fig = gcf;
%         ax.FontSize = 12;
%         fig.Color = 'w';
%     end
%     if row>(n-7)
%         xlabel('N back matches')
%     end
%     if mod(row,7)==1
%         ylabel('BDM value')
%     end
%     if row == n
%         legend({'1-back','2-back'},'Location','Best')
%     end
%     linkaxes
% end

n1 = data.task_displayed == 1;
n2 = data.task_displayed == 2;
ndetect = data.task_displayed==7;
figure
subplot(2,2,1)
for subj = 1:n
    diff(subj) = [nanmean(data.values(subj,n2(subj,:)))-nanmean(data.values(subj,n1(subj,:)))];
    scatter(diff(subj),data.NFC(subj))
    hold on
end
[r,p] = corr(diff(~isnan(data.NFC))',data.NFC(~isnan(data.NFC)),'Type','Spearman');
title([num2str(r) ' : relationship of NFC and mean diff between task values'])
ylabel('NFC')
xlabel('2-back minus 1-back fair wages')
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

subplot(2,2,2)
for subj = 1:n
    scatter(diff(subj),nanstd(data.values(subj,:)))
    hold on
end
title('Relationship of STD of BDMs and mean diff between task values')
ylabel('STD of BDM fair wages (all from 1 subj)')
xlabel('2-back minus 1-back fair wages')
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

% plot mean fair wage on each task by tertile split groups
measures = [data.NFC data.SAPS];
labels = {'NFC','SAPS'};
for measure = 1:2
    split = [];
    split = tertileSplit(measures(:,measure));
    if measure == 1
        colors = NFCcolors;
    else
        colors = SAPScolors;
    end
    subplot(2,2,2+measure)
    lowNFCvalues = data.values(split==1,:);
    midNFCvalues = data.values(split==2,:);
    highNFCvalues = data.values(split==3,:);
    errorbar([nanmean(lowNFCvalues(n1(split==1,:))) nanmean(lowNFCvalues(ndetect(split==1,:))) nanmean(lowNFCvalues(n2(split==1,:)))],[nanstd(lowNFCvalues(n1(split==1,:))) nanstd(lowNFCvalues(ndetect(split==1,:))) nanstd(lowNFCvalues(n2(split==1,:)))]/sqrt(sum(split==1)),'Color',colors(1,:),'Linewidth',1.5)
    hold on
    errorbar([nanmean(midNFCvalues(n1(split==2,:))) nanmean(midNFCvalues(ndetect(split==2,:))) nanmean(midNFCvalues(n2(split==2,:)))],[nanstd(midNFCvalues(n1(split==2,:))) nanstd(midNFCvalues(ndetect(split==2,:))) nanstd(midNFCvalues(n2(split==2,:)))]/sqrt(sum(split==2)),'Color',colors(2,:),'Linewidth',1.5)
    errorbar([nanmean(highNFCvalues(n1(split==3,:))) nanmean(highNFCvalues(ndetect(split==3,:))) nanmean(highNFCvalues(n2(split==3,:)))],[nanstd(highNFCvalues(n1(split==3,:))) nanstd(highNFCvalues(ndetect(split==3,:))) nanstd(highNFCvalues(n2(split==3,:)))]/sqrt(sum(split==3)),'Color',colors(3,:),'Linewidth',1.5)
    legend(['Low ' labels{measure}],['Mid ' labels{measure}],['High ' labels{measure}])
    title(['Mean Task Wage Requested by ' labels{measure} ' group'])
    ylabel('Wage')
    xlabel('Task')
    xticklabels(tasklabels(2:end))
    xticks([1:3])
    ax = gca; fig = gcf;
    ax.FontSize = 12; fig.Color = 'w';

    tasklist = [repmat(1,n,1); repmat(2,n,1); repmat(3,n,1)];
    ratings = [nanmean(n1subjvalue,2); nanmean(n3subjvalue,2); nanmean(n2subjvalue,2)];
    matrix = [ratings repmat(split,3,1) tasklist];
    %[~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{'NFC','task'});
end

%% modeling results x individual differences
% first, NFC by linear mixed effects
split = tertileSplit(data.NFC);

figure
subplot(1,2,1)
bar(1,nanmean(data.n2intercept(split==1)),'FaceColor',NFCcolors(1,:));
hold on
bar(2,nanmean(data.n2intercept(split==2)),'FaceColor',NFCcolors(2,:)); 
bar(3,nanmean(data.n2intercept(split==3)),'FaceColor',NFCcolors(3,:));
errorbar([nanmean(data.n2intercept(split==1)) nanmean(data.n2intercept(split==2)) nanmean(data.n2intercept(split==3))], ...
    [nanstd(data.n2intercept(split==1)) nanstd(data.n2intercept(split==2)) nanstd(data.n2intercept(split==3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))], ...
    'k','LineWidth',1.5)
scatter(ones(sum(split==1),1),data.n2intercept(split==1),[],'k','Filled')
scatter(2.*ones(sum(split==2),1),data.n2intercept(split==2),[],'k','Filled')
scatter(3.*ones(sum(split==3),1),data.n2intercept(split==3),[],'k','Filled')
ylabel({'Mean Intercepts'})
title('Mixed Effect Intercepts by NFC group')
xticks([1:3])
xticklabels({'Low NFC','Mid NFC','High NFC'})
ax = gca; fig = gcf; 
fig.Color = 'w'; ax.FontSize = 12;

subplot(1,2,2)
bar(1,nanmean(data.n2slope(split==1)),'FaceColor',NFCcolors(1,:));
hold on
bar(2,nanmean(data.n2slope(split==2)),'FaceColor',NFCcolors(2,:)); 
bar(3,nanmean(data.n2slope(split==3)),'FaceColor',NFCcolors(3,:));
errorbar([nanmean(data.n2slope(split==1)) nanmean(data.n2slope(split==2)) nanmean(data.n2slope(split==3))], ...
    [nanstd(data.n2slope(split==1)) nanstd(data.n2slope(split==2)) nanstd(data.n2slope(split==3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))], ...
    'k','LineWidth',1.5)
scatter(ones(sum(split==1),1),data.n2slope(split==1),[],NFCcolors(1,:),'Filled')
scatter(2.*ones(sum(split==2),1),data.n2slope(split==2),[],NFCcolors(2,:),'Filled')
scatter(3.*ones(sum(split==3),1),data.n2slope(split==3),[],NFCcolors(3,:),'Filled')
ylabel({'Mean Slopes'})
title('Mixed Effect Slopes by NFC group')
xticks([1:3])
xticklabels({'Low NFC','Mid NFC','High NFC'})
ax = gca; fig = gcf; 
fig.Color = 'w'; ax.FontSize = 12;

[h,p] = ttest2(data.n2intercept(split==1),data.n2intercept(split==3))
[h,p] = ttest2(data.n2slope(split==1),data.n2slope(split==3))

%is the high high NFC intercept a competence effect?
idx = ~isnan(data.n2intercept)&split==3;
[r,~] = corr(data.n2intercept(idx),data.practiceacc(idx,3));
practicereps = sum(~isnan(data.practiceacc),2);
[r,p] = corr(data.n2intercept(idx),practicereps(idx));

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
   if p<0.05
       disp(['relationship between ' tasklabels{task} ' overall acc and diff rating r = ' num2str(r) ' p = ' num2str(p)])
   end
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

prune = [data.n2intercept, data.diffrating];
idx = isnan(data.n2intercept)|isnan(data.diffrating);
prune(idx,:) = [];
[r,p] = corr(prune(:,1),prune(:,2));
figure
subplot(1,2,1)
scatter(prune(:,1),prune(:,2),'Filled')
lsline
title('Intercept on N2 BDM Points by Overall Difficulty Rating')
xlabel('Subject Intercept')
ylabel('Subject Difficulty Rating (1 easy, 5 hard)')

subplot(1,2,2)
iters = sum(~isnan(n2subjvalue),2);
y = data.n2intercept + (data.n2slope.*iters); 
scatter(data.diffrating,y,'Filled')
lsline
title('Final n2 Slope Effect by Difficulty Rating')
ylabel('Final n2 value')
xlabel('Difficulty Rating at End')
fig = gcf; fig.Color = 'w';
idx = isnan(y)|isnan(data.diffrating);
pruned = [y, data.diffrating];
pruned(idx,:) = [];
[r,p] = corr(pruned(:,1),pruned(:,2));

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
toplot = countEntries(NFCdata{group});
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

modelfits = load('data/modelfits_2_models_13-Jan-2021.mat'); %miss c, response c, lure c - most recoverable
%modelfits = load('data/modelfits_3_models_09-Jan-2021_b.mat'); %uc, responsec, lurec - best fitting model thus far
best_model = modelfits.modelStruct.best_model;
nparams = best_model.nparams; 
params = best_model.lowparams; paramnames = best_model.paramnames;
measures = [data.NFC data.SAPS];
names = {'NFC','SAPS'};
for measure = 1:size(measures,2)
split = tertileSplit(measures(:,measure));
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
    disp([paramnames{p} ': ' num2str(pval) 'difference between low and mid ' names{measure}])
    [h,pval] = ttest2(params(split==3,p),params(split==2,p));
    disp([paramnames{p} ': ' num2str(pval) 'difference between high and mid ' names{measure}])
    [h,pval] = ttest2(params(split==1,p),params(split==3,p));
    disp([paramnames{p} ': ' num2str(pval) 'difference between low and high ' names{measure}])
end
fig = gcf; fig.Color = 'w';
end

%% why high dropout rates?
figure
subplot(1,2,1)
% plot for experiment 1, no longer relevant here
% bar([sum(excluded.values(excluded.exp_version==1,1:4)) 0 n1])
% labels = [excluded.labels(1,1:4) ' ' 'finished'];
% xticklabels([labels])
% title('Dropout over course of experiment 1')
% ylabel('n')
% ax = gca; fig = gcf;
% ax.FontSize = 14;
% fig.Color = 'w';
for i = 1:height(excluded)
    curve = excluded.practice_accuracy(i,:);
    plot(1:length(curve),curve)
    hold on
end
title('Accuracy by Practice #, for excluded subjects')
ylabel('Accuracy')
xlabel('Block #')
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

subplot(1,2,2)
bar([sum(excluded.values(excluded.exp_version==4,1:4)) 0 n])
labels = [excluded.labels(1,1:4) ' ' 'finished'];
xticklabels([labels])
title('Dropout over course of experiment 3')
ylabel('n')
ax = gca; fig = gcf;
ax.FontSize = 14;
fig.Color = 'w';

%% Timing analyses for 1-detect task - no longer relevant but may be again
figure
subplot(2,2,1)
detecttiming = [];
for subj = 1:n
    thing = data.detect_time_breakdown{subj}; thing(:,1) = subj;
    detecttiming = [detecttiming; thing];
    scatter(subj,nanmean(thing(:,end)),'Filled')
    hold on
end
title('Mean TOT for each subject')

subplot(2,2,2)
for subj = 1:n
    thing = data.detect_time_breakdown{subj}; 
    scatter(subj,nanmean(thing(:,2)),'Filled')
    hold on
end
title('Mean time on trials for each subject')

subplot(2,2,3)
for subj = 1:n
    thing = data.detect_time_breakdown{subj};
    scatter(subj,nanmean(thing(:,3)),'Filled')
    hold on
end
title('Mean time on feedback for each subject')

subplot(2,2,4)
for subj = 1:n
    thing = data.detect_time_breakdown{subj};
    scatter(subj,nanmean(thing(:,4)),'Filled')
    hold on
end
title('Mean time on ITI for each subject')

figure
for example = 1:length(detecttiming)
    row = detecttiming(example,:);
    subj = row(1);
    color = subjcolors(subj,:);
    %scatter(example,row(:,2)/1000,[],color,'Filled') %how much is trial timing contributing to total timing?
    hold on
    point = (row(:,3)+row(:,2))/1000;
    scatter(example,point,[],color,'Filled') %how much is feedback timing contributing?
    point = (row(:,3)+row(:,2)+row(:,4))/1000;
    %scatter(example,point,[],color,'Filled') %how much is ITI timing contributing to outliers?
    scatter(example,row(:,5),'k','Filled')
end

figure
[r,p] = corr(detecttiming(:,6),detecttiming(:,5))
scatter(detecttiming(:,5),detecttiming(:,6),'Filled')
title('Does num responses have an effect on timing?')
xlabel('TOT')
ylabel('# responses that block')
if p < 0.05
    lsline
end

