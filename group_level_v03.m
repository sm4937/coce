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
files = {'24.03.2020.mat','26.03.2020.mat','31.03.2020.mat','07.04.2020.mat'}; %% %wrong # of variables here %'06.04.2020.mat'
[data,excluded] = loadcostdata(files); %load data in
n = height(data);

tasks = [categorical(cellstr('detection'));categorical(cellstr('n1')); categorical(cellstr('n2'))];
tasklabels = {'0-back','1-back','2-back'};
default_length = 32;

n1 = sum(version==1);
n2 = sum(version==2);
n3 = sum(version==3);%specific n's for graphing/useful to have
inattentive = data.perf<70;

% save important variables for later
for subj = 1:n
    for task = 1:length(tasks)
        list = find(data.task_progression(subj,:)==tasks(task));
        tasks_overall(subj,task) = nanmean(data.perf(subj,list));
    end
end
tasks_rts = [nanmean(data.detectrts,2) nanmean(data.n1rts,2) nanmean(data.n2rts,2)];

data = run_stats(data);

data.taskfreqs = [sum(data.task_progression==tasks(1),2) sum(data.task_progression==tasks(2),2) sum(data.task_progression==tasks(3),2)];
%% Check data quality - how long each task takes, how many times they've been completed
figure
subplot(1,2,1)
y = data.TOT(data.task_progression==tasks(1));
x = data.TOT(data.task_progression==tasks(2));
z = data.TOT(data.task_progression==tasks(3));
scatter(ones(length(y),1),y)
hold on
scatter(2*ones(length(x),1),x)
scatter(3*ones(length(z),1),z)
% Fix outliers - where did those come from?
%[h,p] = ttest2(z,y);
xlim([0.75 3.25])
ylabel('Time on task (seconds)')
xticks([1:3])
xticklabels({'0-back','1-back','2-back'})
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;

[h,p] = ttest2(data.TOT(data.task_progression==tasks(1)),data.TOT(data.task_progression==tasks(2)));
[h,p] = ttest2(data.TOT(data.task_progression==tasks(3)),data.TOT(data.task_progression==tasks(1)));

%prune outliers from task 1 - does that fix it?
dettotstd = nanstd(data.TOT(data.task_progression==tasks(1)));
dettotmean = nanmean(data.TOT(data.task_progression==tasks(1)));

outliers = abs(data.TOT-dettotmean)>(3.*dettotstd);
trimmed = data.TOT;
trimmed(outliers) = dettotmean;

[h,p] = ttest2(trimmed(data.task_progression==tasks(1)),trimmed(data.task_progression==tasks(2)));
[h,p] = ttest2(trimmed(data.task_progression==tasks(3)),trimmed(data.task_progression==tasks(1)));

%there is a significant difference in how long it takes to do task 1 - fix
%in version 4 - it is only 400 milliseconds on average but it cannot
%persist

%no difference in tasks 2 and 3 (n-backs)

errorbar(1:3,[mean(y) mean(x) mean(z)],[std(y)/sqrt(length(y)) std(x)/sqrt(length(x)) std(z)/sqrt(length(z))],'k')

subplot(1,2,2)
y = sum(data.taskfreqs(:,1));
x = sum(data.taskfreqs(:,2));
w = sum(data.taskfreqs(:,3));
%bar([y,x,w,z])
bar([y,x,w])
%xticklabels({'detection','combine','n-switch,p0.1','n-switchp0.9'})
xticklabels({'0-back','1-back','2-back'})
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title('# times each task completed across subjects')
ylabel('Count')

%Mean accuracy by task, mean RTs, make sure "harder" tasks are actually
%harder
figure
subplot(1,3,1)
fig = gcf;
fig.Color = 'w';
bar([nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3))])
hold on
errorbar(1:length(tasklabels),[nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3))],[nanstd(tasks_overall(:,1))/sqrt(n) nanstd(tasks_overall(:,2))/sqrt(n) nanstd(tasks_overall(:,3))/sqrt(n)],'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
ylim([50 100])
title(['Accuracy by task'])

subplot(1,3,2)
allsubjs = [];
for subj = 1:n
   line = [nanmean(data.detectacc(subj,:)) nanmean(data.n1acc(subj,:)) nanmean(data.n2acc(subj,:))];
   allsubjs = [allsubjs; line];
   plot(1:length(line),line,'--')
   hold on
end
errorbar(nanmean(allsubjs),nanstd(allsubjs)/sqrt(n),'k','LineWidth',1.5)
title('Accuracy subj by subj for each task')
xticklabels(tasklabels)
fig = gcf;
fig.Color = 'w';
xlim([0.75 3.25])
missingvals = sum(isnan(allsubjs),2)>0;
allsubjs(missingvals,:) = []; %prune for stats
[r,p] = corr(allsubjs(:,1),allsubjs(:,2));
[r,p] = corr(allsubjs(:,2),allsubjs(:,3));
[r,p] = corr(allsubjs(:,1),allsubjs(:,3));

subplot(1,3,3)
fig = gcf;
fig.Color = 'w';
bar([nanmean(tasks_rts(:,1)) nanmean(tasks_rts(:,2)) nanmean(tasks_rts(:,3))])
hold on
errorbar(1:length(tasklabels),[nanmean(tasks_rts(:,1)) nanmean(tasks_rts(:,2)) nanmean(tasks_rts(:,3))],[nanstd(tasks_rts(:,1))/sqrt(n) nanstd(tasks_rts(:,2))/sqrt(n) nanstd(tasks_rts(:,3))/sqrt(n)],'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
title('Mean RT by task')

%just check out subject BDM strategy
figure
subplot(2,2,1)
for i = 1:n
    scatter(1:default_length,data.values(i,:),'o','Filled')
    hold on
end
%errorbar(nanmean(data.values),nanstd(data.values)/sqrt(n),'k','LineWidth',1)
title('Mean fair wage per subject per block')
fig = gcf; ax = gca;
fig.Color = 'w';
ax.FontSize = 12;
legend({'Subj 1','Subj 2','...','Subj n'})
ylabel('BDM points requested')
xlabel('Block')

subplot(2,2,2)
n1 = data.task_displayed == 1;
n2 = data.task_displayed == 2;
bar([nanmean(data.values(n1)) nanmean(data.values(n2))],'BarWidth',0.5)
hold on
for subj = 1:n
    errorbar([nanmean(data.values(subj,n1(subj,:))) nanmean(data.values(subj,n2(subj,:)))],[nanstd(data.values(subj,n1(subj,:)))/sqrt(sum(n1(subj,:))) nanstd(data.values(subj,n2(subj,:)))/sqrt(sum(n1(subj,:)))],'LineWidth',1.25)
end
xticks([1:2])
xlim([0.75 2.25])
xticklabels({'1-back','2-back'})
title('Mean BDM points by task displayed')
legend({'Task Means','Subj 1','Subj 2','...','Subj n'})
ylabel('Mean BDM points requested')
xlabel('Task')

[h,p] = ttest(data.values(n1),data.values(n2))
disp ('t-test task 1 versus task 2 BDM values') %task 1 values versus task 2 values

subplot(2,2,3)
n1subjlearning = NaN(n,default_length+1); n2subjlearning = NaN(n,default_length+1);
n1subjvalue = NaN(n,default_length+1); n2subjvalue = NaN(n,default_length+1);
for row = 1:n %cycle through subjects
    for task = 1:2 %cycle through tasks
        for trial = 2:default_length
            iters = sum(data.task_progression(row,1:(trial-1)) == tasks(task+1))+1;
            if data.task_progression(row,trial) == tasks(task+1)
                eval(['n' num2str(task) 'subjlearning(row,iters) = data.perf(row,trial);'])
                eval(['n' num2str(task) 'subjvalue(row,iters) = data.values(row,trial);'])
            end
        end
    end
end
errorbar(nanmean(n1subjlearning),nanstd(n1subjlearning)/sqrt(n),'b','LineWidth',2)
hold on
errorbar(nanmean(n2subjlearning),nanstd(n2subjlearning)/sqrt(n),'r','LineWidth',2)
title('Task performance per task iteration')
legend({'1-back','2-back'})
fig = gcf; ax = gca;
fig.Color = 'w';
ax.FontSize = 12;
ylabel('Accuracy')
xlabel('Task iteration')
xticks([0:20])

subplot(2,2,4)
errorbar(nanmean(n1subjvalue),nanstd(n1subjvalue)/sqrt(n),'b','LineWidth',2)
hold on
errorbar(nanmean(n2subjvalue),nanstd(n2subjvalue)/sqrt(n),'r','LineWidth',2)
title('BDM by task iteration')
legend({'1-back','2-back'})
fig = gcf; ax = gca;
fig.Color = 'w';
ax.FontSize = 12;
ylabel('BDM points requested')
xlabel('Task iteration')
xticks([0:20])

figure
%what does this bdm x iter thing look like subj by subj?
for task = 1:2
    subplot(1,2,task)
    eval(['toplot = n' num2str(task) 'subjvalue;'])
    for subj = 1:n
        plot(1:default_length+1,toplot(subj,:))
        hold on
    end
    title(['BDM by iteration: ' num2str(task) '-back'])
end

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
subplot(1,3,1)
ax = gca; fig = gcf;
for i = 1:n
    plot(data.detectacc(i,:),'bo')
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
xlabel('Block #')
ylabel('Accuracy')

subplot(1,3,2)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.n1acc(i,:),'ro')
end
errorbar(nanmean(data.n1acc),nanstd(data.n1acc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (1-back)')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylabel('Accuracy')

subplot(1,3,3)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.n2acc(i,:),'ro')
end
errorbar(nanmean(data.n2acc),nanstd(data.n2acc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (2-back)')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylabel('Accuracy')
linkaxes

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
end

%% Understand n-back matches' effect on costs

%first, plot general nbackmatches effect for all subjects
figure
matrix = [];
for i=1:n
    task_list = data.task_progression(i,:);
    if ismember(categorical({'n1'}),task_list) %version 2
        nback = find(task_list==tasks(2)|task_list==tasks(3)&~inattentive(i,:));
        next = nback+1;
        l = length(task_list);
        y = data.values(i,next(next~=(l+1)));
        x = data.nbackmatches{i}(nback(nback~=l)); %exclude last trial for sizing reasons
        %costs = NFC_sigmoid(x,data.NFC(i));
        matrix = [matrix; y' x'];
        scatter(x,y,'o','Filled');
        hold on
    end
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
disp(['Nswitches vs. BDM request. \ Spearman rho: ' num2str(r) ', p value: ' num2str(p)])

same_long = []; diff_long = [];
figure
for task = 1:2 %cycle through what is being displayed, display by what happened before
    subplot(3,2,task)
    for subj = 1:n %go subject by subject
        displayed = data.task_displayed(subj,:)+1;
        idx = displayed == (task+1);
        matches = [NaN data.nbackmatches{subj,:}]; matches(end) = []; %last trial not influential
        BDMs = data.values(subj,:);
        completed = NaN(1,default_length);completed(data.task_progression(subj,:)==tasks(2)) = 2;completed(data.task_progression(subj,:)==tasks(3)) = 3;
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
    title(['BDM for ' num2str(task) '-back task vs. matches in prev round'])
    xlabel('# of n-back matches in last block')
end
fig = gcf;
fig.Color = 'w';
ylim([1 5.1])

n1effect = []; n2effect = []; % keep track of numbers pulled out here for stats later
for task = 1:2 %cycle through what is being displayed, display by what happened before
    subplot(3,2,task+2)
    for subj = 1:n %go subject by subject
        displayed = data.task_displayed(subj,:)+1;
        idx = find(displayed == (task+1)); 
        matches = data.nbackmatches{subj,:};
        misses = data.nbackmisses(subj,:);
        perf = data.perf(subj,:);
        BDMs = data.values(subj,:);
        completed = find(data.task_progression(subj,:)==tasks(task+1));
        for trial = 1:length(idx)
            now = idx(trial);
            if sum(completed<now)>0 %they've done the last once before
                last = completed(completed<now); last = last(end);
                delay = now-last;
                scatter(matches(last),BDMs(now),'Filled')
                eval(['n' num2str(task) 'effect = [n' num2str(task) 'effect; subj matches(last) misses(now) delay BDMs(now) perf(now)];'])
            end
            hold on
        end
    end
    ylabel('BDM points')
    title(['BDM for ' num2str(task) '-back task'])
    xlabel(['n-back matches in last block of ' num2str(task) '-back'])
end
fig = gcf;
fig.Color = 'w';
ylim([1 5.1])

[r,p] = corr(n1effect(:,2),n1effect(:,5),'Type','Spearman');
[r,p] = corr(n2effect(:,2),n2effect(:,5),'Type','Spearman');

for task = 1:2
    subplot(3,2,task+4)
    for subj = 1:n
        eval(['y = data.n' num2str(task) 'matcheffect(subj,:);'])
        plot(1:3,y,'LineWidth',1.5)
        hold on
    end
    title(['Match "effect" on BDM from last ' num2str(task) '-back'])
    eval(['y = data.n' num2str(task) 'matcheffect;'])
    errorbar(nanmean(y),nanstd(y)/sqrt(n),'k','LineWidth',1.5)
    xticklabels({'3','4','5'})
    xlabel('Matches in Last Round of Task')
    legend({'Subj 1','Subj 2','...','Subj n'})
end

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
histogram(data.NFC,'BinWidth',0.5)
title('Hist of NFC scores (normalized)')
subplot(1,3,2)
histogram(data.SAPS,'BinWidth',0.25)
title('Hist of SAPS scores (normalized)')
meanBDM = nanmean(data.values(:,2:end),2);
subplot(1,3,3)
scatter(data.NFC,data.SAPS,'Filled')
title('NFC by SAPS score')
xlabel('NFC score')
ylabel('SAPS score')
[r,p] = corr(data.NFC(~isnan(data.NFC)),data.SAPS(~isnan(data.NFC)),'Type','Spearman');
if p<0.05
    hold on
    lsline
end

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
    nanmean(tasks_overall(split==1,3)) nanmean(tasks_overall(split==2,3)) nanmean(tasks_overall(split==3,3))];
E = [nanstd(tasks_overall(split==1,1)) nanstd(tasks_overall(split==2,1)) nanstd(tasks_overall(split==3,1)); ...
    nanstd(tasks_overall(split==1,2)) nanstd(tasks_overall(split==2,2)) nanstd(tasks_overall(split==3,2)); ...
    nanstd(tasks_overall(split==1,3)) nanstd(tasks_overall(split==2,3)) nanstd(tasks_overall(split==3,3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:9,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:9,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:9,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
title('Accuracy by SAPS group')
xticks([2:3:9])
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
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
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
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
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

figure
for task = 1:2
    y = data.nbackfurthererrors;
    relevant = y(data.task_progression==tasks(task+1));
    subplot(2,2,task)
    histogram(relevant,'BinWidth',0.25)
    title(['Distribution of ER after first error - ' num2str(task) '-back'])
    xlabel('Post-Error Error Rate')
    for group = 1:3
        subplot(2,2,task+2)
        mini = y(split==group,:);
        relevant = mini(data.task_progression(split==group,:)==tasks(task+1));
        %means(group) = nanmean(relevant); E(group) = nanstd(relevant);
        %bar(group,nanmean(relevant),'FaceColor',SAPScolors(group,:))
        histogram(relevant,'BinWidth',0.10,'FaceColor',SAPScolors(group,:))
        hold on
        title(['ER after first error - ' num2str(task) '-back'])
        ylabel('Post-Error Error Rate')
        legend({'Low SAPS','Mid SAPS','High SAPS'})
    end
    %errorbar(means,E,'k','LineWidth',1.5)
end %of task loop
fig = gcf; fig.Color = 'w';

figure
%slower performance, higher perfectionism?
bars = [nanmean(tasks_rts(split==1,1)) nanmean(tasks_rts(split==2,1)) nanmean(tasks_rts(split==3,1)); ...
    nanmean(tasks_rts(split==1,2)) nanmean(tasks_rts(split==2,2)) nanmean(tasks_rts(split==3,2)); ...
    nanmean(tasks_rts(split==1,3)) nanmean(tasks_rts(split==2,3)) nanmean(tasks_rts(split==3,3))];
E = [nanstd(tasks_rts(split==1,1)) nanstd(tasks_rts(split==2,1)) nanstd(tasks_rts(split==3,1)); ...
    nanstd(tasks_rts(split==1,2)) nanstd(tasks_rts(split==2,2)) nanstd(tasks_rts(split==3,2)); ...
    nanstd(tasks_rts(split==1,3)) nanstd(tasks_rts(split==2,3)) nanstd(tasks_rts(split==3,3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:9,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:9,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:9,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticks([2:3:9])
xticklabels(tasklabels)
title('Mean RTs by SAPS group')
legend({'Low SAPS','Mid SAPS','High SAPS'})
xticklabels(tasklabels)

figure
subplot(1,2,1)
pracn0acc = data.practiceacc(:,1);
pracn1acc = data.practiceacc(:,2);
pracn2acc = data.practiceacc(:,3); %first practice round for each task
bars = [nanmean(pracn0acc(split==1,:)) nanmean(pracn0acc(split==2,:)) nanmean(pracn0acc(split==3,:)); ...
    nanmean(pracn1acc(split==1,:)) nanmean(pracn1acc(split==2,:)) nanmean(pracn1acc(split==3,:)); ...
    nanmean(pracn2acc(split==1,:)) nanmean(pracn2acc(split==2,:)) nanmean(pracn2acc(split==3,:))];
E = [nanstd(pracn0acc(split==1,:)) nanstd(pracn0acc(split==2,:)) nanstd(pracn0acc(split==3,:)); ...
    nanstd(pracn1acc(split==1,:)) nanstd(pracn1acc(split==2,:)) nanstd(pracn1acc(split==3,:)); ...
    nanstd(pracn2acc(split==1,:)) nanstd(pracn2acc(split==2,:)) nanstd(pracn2acc(split==3,:))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:9,bars(:,1),'FaceColor',SAPScolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:9,bars(:,2),'FaceColor',SAPScolors(2,:),'BarWidth',0.2)
bar(3:3:9,bars(:,3),'FaceColor',SAPScolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticks([2:3:9])
xticklabels(tasklabels)
legend({'Low SAPS','Mid SAPS','High SAPS'})
title('Interaction of SAPS group with practice accuracy')
ax = gca; fig = gcf; 
fig.Color = 'w'; ax.FontSize = 12;

% look into self-rated competence i.e. how "easy" or "hard" they rated the
% task to be
subplot(1,2,2)
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

%% Tertile split NFC scores
%style variables
NFCcolors = [.5 0 .5; .75 0 .75; .95 0 .95];
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
    x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
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
    x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
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
    x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
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

n1effect = []; n2effect = []; % keep track of numbers pulled out here for stats later
figure
for task = 1:2 %cycle through what is being displayed, display by what happened before
    subplot(1,3,task)
    for subj = 1:n %go subject by subject
        group = split(subj);
        displayed = data.task_displayed(subj,:)+1;
        idx = find(displayed == (task+1)); 
        matches = data.nbackmatches{subj,:};
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
    title(['BDM for ' num2str(task) '-back task'])
    legend({'Mid NFC','Low NFC','High NFC'})
    xlabel(['n-back matches in last block of ' num2str(task) '-back'])
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

subplot(1,3,3)
% accuracy by group
bars = [nanmean(tasks_overall(split==1,1)) nanmean(tasks_overall(split==2,1)) nanmean(tasks_overall(split==3,1)); ...
    nanmean(tasks_overall(split==1,2)) nanmean(tasks_overall(split==2,2)) nanmean(tasks_overall(split==3,2)); ...
    nanmean(tasks_overall(split==1,3)) nanmean(tasks_overall(split==2,3)) nanmean(tasks_overall(split==3,3))];
E = [nanstd(tasks_overall(split==1,1)) nanstd(tasks_overall(split==2,1)) nanstd(tasks_overall(split==3,1)); ...
    nanstd(tasks_overall(split==1,2)) nanstd(tasks_overall(split==2,2)) nanstd(tasks_overall(split==3,2)); ...
    nanstd(tasks_overall(split==1,3)) nanstd(tasks_overall(split==2,3)) nanstd(tasks_overall(split==3,3))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:9,bars(:,1),'FaceColor',NFCcolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:9,bars(:,2),'FaceColor',NFCcolors(2,:),'BarWidth',0.2)
bar(3:3:9,bars(:,3),'FaceColor',NFCcolors(3,:),'BarWidth',0.2)
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticks([2:3:9])
xticklabels(tasklabels)
title('Accuracy by NFC group')
legend({'Low NFC','Mid NFC','High NFC'})
xticklabels(tasklabels)

% how does NFC interact with 'baseline EF' i.e. practice accuracy?
figure
pracn0acc = data.practiceacc(:,1);
pracn1acc = data.practiceacc(:,2);
pracn2acc = data.practiceacc(:,3); %first practice round for each task
bars = [nanmean(pracn0acc(split==1,:)) nanmean(pracn0acc(split==2,:)) nanmean(pracn0acc(split==3,:)); ...
    nanmean(pracn1acc(split==1,:)) nanmean(pracn1acc(split==2,:)) nanmean(pracn1acc(split==3,:)); ...
    nanmean(pracn2acc(split==1,:)) nanmean(pracn2acc(split==2,:)) nanmean(pracn2acc(split==3,:))];
E = [nanstd(pracn0acc(split==1,:)) nanstd(pracn0acc(split==2,:)) nanstd(pracn0acc(split==3,:)); ...
    nanstd(pracn1acc(split==1,:)) nanstd(pracn1acc(split==2,:)) nanstd(pracn1acc(split==3,:)); ...
    nanstd(pracn2acc(split==1,:)) nanstd(pracn2acc(split==2,:)) nanstd(pracn2acc(split==3,:))]./[sqrt(sum(split==1)) sqrt(sum(split==2)) sqrt(sum(split==3))];
bar(1:3:9,bars(:,1),'FaceColor',NFCcolors(1,:),'BarWidth',0.2)
hold on
bar(2:3:9,bars(:,2),'FaceColor',NFCcolors(2,:),'BarWidth',0.2)
bar(3:3:9,bars(:,3),'FaceColor',NFCcolors(3,:),'BarWidth',0.2)
xticks([2:3:9])
bars = reshape(bars',1,numel(bars));
E = reshape(E',1,numel(E));
errorbar(bars,E,'*k')
xticklabels(tasklabels)
legend({'Low NFC','Mid NFC','High NFC'})
title('Interaction of NFC group with practice accuracy')
ax = gca; fig = gcf; 
fig.Color = 'w'; ax.FontSize = 12;

%% try to understand random effects by plotting one plot/subject, BDM vs.
%nswitches, NFC as title %%%%%%%
split = tertileSplit(data.NFC);

figure
for row = 1:n
    subplot(6,7,row)
    for task = 2:3
        nback = find(data.task_progression(row,:)==tasks(task));
        nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
        next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
        x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
        y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
        if task == 2
            color = 'ob';
        elseif task == 3
            color = 'or';
        end
        scatter(x,y,color,'Filled')
        hold on
        title(num2str(data.NFC(row)))
        ax = gca; fig = gcf;
        ax.FontSize = 12;
        fig.Color = 'w';
    end
    if row>(n-7)
        xlabel('N back matches')
    end
    if mod(row,7)==1
        ylabel('BDM value')
    end
    if row == n
        legend({'1-back','2-back'},'Location','Best')
    end
    linkaxes
end

n1 = data.task_displayed == 1;
n2 = data.task_displayed == 2;
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

subplot(2,2,3)
lowNFCvalues = data.values(split==1,:);
midNFCvalues = data.values(split==2,:);
highNFCvalues = data.values(split==3,:);
errorbar([nanmean(lowNFCvalues(n1(split==1,:))) nanmean(lowNFCvalues(n2(split==1,:)))],[nanstd(lowNFCvalues(n1(split==1,:))) nanstd(lowNFCvalues(n2(split==1,:)))]/sqrt(sum(split==1)),'c','Linewidth',1.5)
hold on
errorbar([nanmean(midNFCvalues(n1(split==2,:))) nanmean(midNFCvalues(n2(split==2,:)))],[nanstd(midNFCvalues(n1(split==2,:))) nanstd(midNFCvalues(n2(split==2,:)))]/sqrt(sum(split==2)),'b','Linewidth',1.5)
errorbar([nanmean(highNFCvalues(n1(split==3,:))) nanmean(highNFCvalues(n2(split==3,:)))],[nanstd(highNFCvalues(n1(split==3,:))) nanstd(highNFCvalues(n2(split==3,:)))]/sqrt(sum(split==3)),'r','Linewidth',1.5)
legend('Low NFC','Mid NFC','High NFC')
title('Mean Task Wage Requested by NFC group')
ylabel('Wage')
xlabel('Task')
xticklabels({'1-back','2-back'})
xticks([1 2])
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

[h,p] = ttest2(lowNFCvalues(n1(split==1,:)),lowNFCvalues(n2(split==1,:)))
disp('t-test for low NFC people')
[h,p] = ttest2(midNFCvalues(n1(split==2,:)),midNFCvalues(n2(split==2,:)))
disp('t-test for mid NFC people')
[h,p] = ttest2(highNFCvalues(n1(split==3,:)),highNFCvalues(n2(split==3,:)))
disp('t-test for high NFC people')

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
[r,p] = corr(data.n2intercept(idx),data.practiceacc(idx,3));
practicereps = sum(~isnan(data.practiceacc),2);
[r,p] = corr(data.n2intercept(idx),practicereps(idx));

%survivor analysis by group as well
figure
for task = 1:2
    E = []; means = [];
    y = data.nbackfurthererrors;
    relevant = y(data.task_progression==tasks(task+1));
    subplot(2,2,task)
    histogram(relevant,'BinWidth',0.25)
    title(['Distribution of ER after first error - ' num2str(task) '-back'])
    xlabel('Post-Error Error Rate')
    for group = 1:3
        subplot(2,2,task+2)
        mini = y(split==group,:);
        relevant = mini(data.task_progression(split==group,:)==tasks(task+1));
        means(group) = nanmean(relevant); E(group) = nanstd(relevant)/sqrt(sum(split==group));
        bar(group,nanmean(relevant),'FaceColor',NFCcolors(group,:))
        hold on
        title(['Mean ER after first error - ' num2str(task) '-back'])
        ylabel('Post-Error Error Rate')
        xticks([1:3])
        xticklabels({'Low NFC','Mid NFC','High NFC'})
    end
    errorbar(means,E,'k','LineWidth',1.5)
end %of task loop
fig = gcf; fig.Color = 'w';

% look into self-rated competence i.e. how "easy" or "hard" they rated the
% task to be
figure
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

subplot(3,2,5)
scatter(data.age,meanBDM)
xlabel('Age')
ylabel('Mean BDM')
title('Age by Mean BDM')

subplot(3,2,6)
scatter(data.age,data.overall)
xlabel('Age')
ylabel('Overall Acc')
title('Age by Overall Accuracy')

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
bar([sum(excluded.values(excluded.exp_version==3,1:4)) 0 n])
labels = [excluded.labels(1,1:4) ' ' 'finished'];
xticklabels([labels])
title('Dropout over course of experiment 3')
ylabel('n')
ax = gca; fig = gcf;
ax.FontSize = 14;
fig.Color = 'w';


