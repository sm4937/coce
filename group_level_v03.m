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
files = {'24.03.2020.mat','26.03.2020.mat'};
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
errorbar(1:3,[mean(y) mean(x) mean(z)],[std(y)/sqrt(length(y)) std(x)/sqrt(length(x)) std(z)/sqrt(length(z))],'ko')
xlim([0.75 3.25])
ylabel('Time on task (seconds)')
xticks([1:3])
xticklabels({'0-back','1-back','2-back'})
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;

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
subplot(1,2,1)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
bar([nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3))])
hold on
errorbar(1:length(tasklabels),[nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3))],[nanstd(tasks_overall(:,1))/sqrt(n) nanstd(tasks_overall(:,2))/sqrt(n) nanstd(tasks_overall(:,3))/sqrt(n)],'k*','LineWidth',2)
xticklabels(tasklabels)
xtickangle(45)
ylim([50 100])
title(['Accuracy by task'])

subplot(1,2,2)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
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
for row = 1:n
    for x = 1:length(data.task_progression(row,:))
        task = data.task_displayed(row,x);
        if task==1
            color = 'b';
        elseif task==2
            color = 'r';
        end
        scatter(x,data.values(row,x),['o' color],'Filled')
        hold on
    end
end
title('Mean fair wage per task per block')
legend({'1-back','2-back'})
fig = gcf; ax = gca;
fig.Color = 'w';
ax.FontSize = 12;
ylabel('BDM points requested')
xlabel('Block')

subplot(2,2,3)
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

[h,p] = ttest(data.values(n1),data.values(n2)); %task 1 values versus task 2 values

subplot(2,2,4)
%plot task avoidance, well, sort of
for i = 1:n
    plot(1:3,data.taskfreqs(i,:),'LineWidth',1)
    hold on
end
ylabel('# of times task completed')
title('Task avoidance by subject')
xlabel('Task')
xticks([1:3])
xticklabels(tasklabels)
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

% plot task performance by block, any learning or decay?
figure
subplot(1,2,1)
ax = gca; fig = gcf;
hold on
for i = 1:n
    plot(data.n2acc(i,:),'ro')
end
errorbar(nanmean(data.n2acc),nanstd(data.n2acc)/sqrt(n),'r','LineWidth',1)
title('Learning Curves (2-back)')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylabel('Accuracy')

subplot(1,2,2)
ax = gca; fig = gcf;
% hold on
% for i = 1:n
%     plot(data.n2rts(i,:),'ro')
% end
% errorbar(nanmean(data.n2rts),nanstd(data.n2rts)/sqrt(n),'k','LineWidth',1)
% title('Learning Curves (2-back)')
% fig.Color = 'w';
% ax.FontSize = 12;
% xlabel('Block #')
% ylabel('Mean log(RT) (ms)')
% removing RT "learning" curve for now
for i = 1:n
    plot(data.detectacc(i,:),'bo')
    hold on
end
errorbar(nanmean(data.detectacc),nanstd(data.detectacc)/sqrt(n),'k','LineWidth',1)
title('Learning Curves (0-back)')
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylabel('Accuracy')
linkaxes

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

%% deal with individual differences in perfectionism, NFC
figure
subplot(1,2,1)
histogram(data.NFC)
title('Hist of NFC scores (normalized)')
subplot(1,2,2)
histogram(data.SAPS)
title('Hist of SAPS scores (normalized)')
meanBDM = nanmean(data.values(:,2:end),2);

figure
subplot(1,2,1)
%plot center BDM requested by NFC score
scores = data.NFC;
scatter(scores,meanBDM,'o','Filled')
xlabel('Need for Cognition Score')
ylabel('Center Mean BDM score')

subplot(1,2,2)
%median split SAPS scores
split = NaN(length(data.SAPS),1);
medi = median(data.SAPS(~isnan(data.SAPS)));
split(data.SAPS>medi) = 2;
split(data.SAPS<=medi) = 1;
%plot accuracy by SAPS score

for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2));
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.perf(row,nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==1
        dots = 'or';
    elseif split(row)==2
        dots = 'ob';
    else
        dots = 'ok';
    end
    scatter(x,y,dots,'Filled')
    hold on
end
legend({'High Perfectionism',' ','No Data','Low Perfectionism'},'Location','Best')
legend('boxoff')
xlabel('Accuracy')
ylabel('BDM value')
xlim([50 100])

%median split NFC scores
split = NaN(length(data.NFC),1);
medi = median(data.NFC(~isnan(data.NFC)));
split(data.NFC>medi) = 2;
split(data.NFC<=medi) = 1;
%plot nswitches vs. BDM as function of NFC score
low = []; high = [];
figure
subplot(1,2,1)
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

subplot(1,2,2)
for row = 1:n
    nback = find(data.task_progression(row,:)==tasks(3)|data.task_progression(row,:)==tasks(2)); 
    nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
    next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
    x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
    y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
    if split(row)==2
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

[r,p] = corr(high(~isnan(high(:,2)),1),high(~isnan(high(:,2)),2),'Type','Spearman'); %high NFC has no relationship of n switches vs. value
[r,p] = corr(low(~isnan(low(:,2)),1),low(~isnan(low(:,2)),2),'Type','Spearman'); %low NFC has relationship p = 0.02 with n = 22

%try to understand random effects by plotting one plot/subject, BDM vs.
%nswitches, NFC as title %%%%%%%

figure
nplot = 0;
for row = 1:n
    for task = 2:3
        nplot = nplot + 1;
        subplot(n,2,nplot)
        nback = find(data.task_progression(row,:)==tasks(task));
        nback(ismember(nback,find(inattentive(row,:)))) = []; %prune blocks where subjects obviously weren't paying attention
        next = nback+1; nback(next==default_length+1) = []; next(next==default_length+1) = [];
        x = data.nbackmatches{row}(nback); %x(24) = []; %last value doesn't have a corresponding BDM value
        y = data.values(row,next); %y(1) = []; %what do they ask for next time? first value doesn't have a corresponding perf level
        if split(row)==2
            color = 'oc';
        elseif split(row) == 1
            color = 'om';
        else
            color = 'ok';
        end
        scatter(x,y,color,'Filled')
        title(num2str(data.NFC(row)))
        ax = gca; fig = gcf;
        ax.FontSize = 12;
        fig.Color = 'w';
        if task == 2
            legend({'1-back'})
        else
            legend({'2-back'})
        end
    end
    if row>(n-5)
        xlabel('N back matches')
    end
    if mod(row,5)==1
        ylabel('BDM value')
    end
end

n1 = data.task_displayed == 1;
n2 = data.task_displayed == 2;
figure
subplot(1,2,1)
for subj = 1:n
    diff(subj) = [nanmean(data.values(subj,n2(subj,:)))-nanmean(data.values(subj,n1(subj,:)))];
    scatter(diff(subj),data.NFC(subj))
    hold on
end
[r,p] = corr(diff(~isnan(data.NFC))',data.NFC(~isnan(data.NFC)));
title([num2str(r) ' : relationship of NFC and mean diff between task values'])
ylabel('NFC')
xlabel('2-back minus 1-back fair wages')
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';

subplot(1,2,2)
lowNFCvalues = data.values(split==1,:);
highNFCvalues = data.values(split==2,:);
errorbar([nanmean(lowNFCvalues(n1(split==1,:))) nanmean(lowNFCvalues(n2(split==1,:)))],[nanstd(lowNFCvalues(n1(split==1,:))) nanstd(lowNFCvalues(n2(split==1,:)))]/sqrt(sum(split==1)),'b','Linewidth',1.5)
hold on
errorbar([nanmean(highNFCvalues(n1(split==2,:))) nanmean(highNFCvalues(n2(split==2,:)))],[nanstd(highNFCvalues(n1(split==2,:))) nanstd(highNFCvalues(n2(split==2,:)))]/sqrt(sum(split==2)),'r','Linewidth',1.5)
legend('Low NFC','High NFC')
title('Mean Task Wage Requested by NFC group')
ylabel('Wage')
xlabel('Task')
xticklabels({'1-back','2-back'})
xticks([1 2])
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';


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


