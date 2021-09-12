%% A script to print stats and paper stuff, with more upfront formatting and 
% fewer exploratory analyses
disp(['Sample of ' num2str(n) ' subjects (' num2str(sum(data.sex==2)) ' female). Mean(std) age: ' num2str(nanmean(data.age)) '(' num2str(nanstd(data.age)) ')'])
disp(['Unspecified sex n: ' num2str(sum(isnan(data.sex))) '; unspecified age: ' num2str(sum(isnan(data.age)))])
disp(['Mean total TOT: ' num2str(nanmean(data.totalTOT)) ', median total TOT: ' num2str(median(data.totalTOT))])
% REMINDER OF TASK ORDERS in VARIABLES %
%tasks = [detection,n1,3detection,n2];

% % MEAN ACCURACY BY TASK % % 
figure
subplot(2,2,1)
fig = gcf; fig.Color = 'w';
errorbar(1:length(tasklabels),[nanmean(tasks_overall(:,1)) nanmean(tasks_overall(:,2)) nanmean(tasks_overall(:,3)) nanmean(tasks_overall(:,4))],[nanstd(tasks_overall(:,1)) nanstd(tasks_overall(:,2)) nanstd(tasks_overall(:,3)) nanstd(tasks_overall(:,4))]/sqrt(n),'k','LineWidth',1.25)
xticks(1:length(tasklabels))
xticklabels(tasklabels)
%xtickangle(45)
xlabel('Task')
ylim([70 100]); xlim([0.5 4.5])
ylabel('Mean accuracy')

% % STATS ON ACCURACY, MEAN RT, and DIFFICULTY RATINGS
% COMPARE MEAN ACCURACY BY TASK %
clear vector; vector = tasks_overall(:); taskidentity = [ones(n,1); 2*ones(n,1); 3*ones(n,1); 4*ones(n,1)];
[~,~,stats] = anovan(vector,taskidentity);
%tasks = [detection,n1,3detection,n2];
Table1 = table; Table1.onedetect = nanmean(vector(taskidentity==1)); Table1.oneback = nanmean(vector(taskidentity==2));
Table1.threedetect = nanmean(vector(taskidentity==3)); Table1.twoback = nanmean(vector(taskidentity==4)); 
Table1.Properties.RowNames{1} = 'Accuracy';
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==4),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==3));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==4));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==4));

% COMPARE MEAN RTs ON TASK %
clear vector; vector = tasks_rts(:); 
[~,~,stats] = anovan(vector,taskidentity,'display','off');
%tasks = [detection,n1,3detection,n2];
temp = table; temp.onedetect = nanmean(vector(taskidentity==1)); temp.oneback = nanmean(vector(taskidentity==2));
temp.threedetect = nanmean(vector(taskidentity==3)); temp.twoback = nanmean(vector(taskidentity==4)); 
Table1 = [Table1; temp]; Table1.Properties.RowNames{2} = 'RT (msec)';
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==4),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==3));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==4));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==4));

% COMPARE DIFFICULTY RATINGS FROM END OF TASK %
clear vector; vector = data.taskratings(:); 
[~,~,stats] = anovan(vector,taskidentity,'display','off');
%tasks = [detection,n1,n2,3detection];
temp = table; temp.onedetect = nanmean(vector(taskidentity==1)); temp.oneback = nanmean(vector(taskidentity==2));
temp.threedetect = nanmean(vector(taskidentity==3)); temp.twoback = nanmean(vector(taskidentity==4)); 
Table1 = [Table1; temp]; Table1.Properties.RowNames{3} = 'Difficulty';
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==4),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==3));
[h,p] = ttest(vector(taskidentity==1),vector(taskidentity==4));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==2));
[h,p] = ttest(vector(taskidentity==3),vector(taskidentity==4));

% COMPARE MEAN BDM RATINGS %
n1 = data.task_displayed==tasknumbers(2); %1
n2 = data.task_displayed==tasknumbers(4); %2
ndetect = data.task_displayed==tasknumbers(3); %7

temp = table; temp.onedetect = NaN; temp.oneback = nanmean(data.values(n1));
temp.threedetect = nanmean(data.values(ndetect)); temp.twoback = nanmean(data.values(n2)); 
Table1 = [Table1; temp]; Table1.Properties.RowNames{4} = 'Fair wage';
[h,p] = ttest2(data.values(n1),data.values(n2));
disp (['t-test 1-back versus 2-back BDM values p = ' num2str(p)]); 
[h,p] = ttest2(data.values(ndetect),data.values(n2));
disp (['t-test 2-back versus 3-detect BDM values p = ' num2str(p)]); %task 2 values versus task 3 values
[h,p] = ttest2(data.values(n1),data.values(ndetect));
disp (['t-test 1-back versus 3-detect BDM values p = ' num2str(p)]); %task 1 values versus task 3 values
% basically it's task 2 versus the world

n0subjlearning = NaN(n,default_length+1); n0subjvalue = NaN(n,default_length+1); n0subjrt = NaN(n,default_length+1);
n1subjlearning = NaN(n,default_length+1); n2subjlearning = NaN(n,default_length+1); n3subjlearning = NaN(n,default_length+1);
n1subjvalue = NaN(n,default_length+1); n2subjvalue = NaN(n,default_length+1); n3subjvalue = NaN(n,default_length+1);
n1subjrt = NaN(n,default_length+1); n2subjrt = NaN(n,default_length+1); n3subjrt = NaN(n,default_length+1);
for row = 1:n %cycle through subjects
    for task = 1:4 %cycle through tasks
        for trial = 2:default_length
            iters = sum(data.task_progression(row,1:(trial-1)) == tasks(task))+1;
            if data.task_progression(row,trial) == tasks(task)
                eval(['n' num2str(task-1) 'subjlearning(row,iters) = data.perf(row,trial);'])
                eval(['n' num2str(task-1) 'subjvalue(row,iters) = data.values(row,trial);'])
                eval(['n' num2str(task-1) 'subjrt(row,iters) = data.BDMrt(row,trial);'])
            end
        end
    end
end

% % MEAN FAIR WAGE BY TASK % %
subplot(2,2,2)
dotsize = 5;
scatter(ones(n,1),nanmean(n1subjvalue,2),dotsize*ones(n,1),taskcolors(2,:),'filled','DisplayName','1-back')
hold on
scatter(2*ones(n,1),nanmean(n2subjvalue,2),dotsize*ones(n,1),taskcolors(3,:),'filled','DisplayName','3-detect')
scatter(3*ones(n,1),nanmean(n3subjvalue,2),dotsize*ones(n,1),taskcolors(4,:),'filled','DisplayName','2-back')

%overlay group means
errorbar([nanmean(data.values(n1)) nanmean(data.values(ndetect)) nanmean(data.values(n2))],[nanstd(data.values(n1)) nanstd(data.values(ndetect)) nanstd(data.values(n2))]./sqrt(n),'k','LineWidth',1.25)

xticks(1:(length(tasklabels)-1))
xlim([0.75 3.25])
ylim([1 5])
xticklabels(tasklabels(2:end))
ylabel('Mean fair wage')
xlabel('Task')

% % ACCURACY BY TASK ITERATION % %
subplot(2,2,3)
errorbar(nanmean(n0subjlearning(:,1:11)),nanstd(n0subjlearning(:,1:11))/sqrt(n),'Color',taskcolors(1,:),'LineWidth',1.25,'DisplayName',tasklabels{1})
hold on
errorbar(nanmean(n1subjlearning),nanstd(n1subjlearning)/sqrt(n),'Color',taskcolors(2,:),'LineWidth',1.25,'DisplayName',tasklabels{2})
errorbar(nanmean(n2subjlearning),nanstd(n2subjlearning)/sqrt(n),'Color',taskcolors(3,:),'LineWidth',1.25,'DisplayName',tasklabels{3})
errorbar(nanmean(n3subjlearning),nanstd(n3subjlearning)/sqrt(n),'Color',taskcolors(4,:),'LineWidth',1.25,'DisplayName',tasklabels{4}) %task 3 is actually ndetect
legend('boxoff')
fig = gcf; ax = gca;
fig.Color = 'w';
ylabel('Accuracy')
xlabel('Iteration #')

% % FAIR WAGE BY TASK ITERATION % % 
subplot(2,2,4)
errorbar(nanmean(n1subjvalue),nanstd(n1subjvalue)./sqrt(n),'Color',taskcolors(2,:),'LineWidth',1.25,'DisplayName',tasklabels{2})
hold on
errorbar(nanmean(n2subjvalue),nanstd(n2subjvalue)./sqrt(n),'Color',taskcolors(3,:),'LineWidth',1.25,'DisplayName',tasklabels{3})
errorbar(nanmean(n3subjvalue),nanstd(n3subjvalue)./sqrt(n),'Color',taskcolors(4,:),'LineWidth',1.25,'DisplayName',tasklabels{4})
ylabel('Fair wage')
legend('boxoff')
xlabel('Iteration #')

% What's the relationship of fair wage and accuracy?
meanBDM = [nanmean(n1subjvalue,2) nanmean(n2subjvalue,2) nanmean(n3subjvalue,2)];
meanAcc = tasks_overall(:,2:4);

[r,p] = corr(meanBDM(:,1),meanAcc(:,1)) %1-back
[r,p] = corr(meanBDM(:,2),meanAcc(:,2)) %3-detect
[r,p] = corr(meanBDM(:,3),meanAcc(:,3)) %2-back

% % % subjects are RELATIVELY STABLE IN THEIR BDMs % % 
% subplot(2,2,4)
% for task = 4 %grab 2-back specifically
%     rating_idx = data.task_displayed==tasknumbers(task);
%     eval(['ratings = n' num2str(task-1) 'subjvalue;'])
%     hold on
%     eval(['completions = sum(~isnan(n' num2str(task-1) 'subjlearning),2);'])
%     ratings = [completions ratings]; ratings = sortrows(ratings,1);
%     n_rounds = unique(completions);
%     for i = 1:length(n_rounds)
%         tomean = ratings(ratings(:,1)==n_rounds(i),2:end);
%         toplot = nanmean(tomean,1);
%         errorbar(toplot,nanstd(tomean,[],1)./sqrt(size(tomean,1)),'Color',taskcolors(task,:))
%         hold on
%     end
%     xlabel('# task iters completed')
%     ylabel('Mean BDM rating')   
%     title([tasklabels{task} ' ratings by completed iterations'])
% end
% fig = gcf; fig.Color = 'w';
% xlim([0 11])

% % BDM RT STUFF % %
% are people getting more decisive on how many BDM points they want?
big_matrix = [];
for task = 1:(length(unique(data.task_displayed(~isnan(data.task_displayed)))))
    matrix = [];
    for subj = 1:n
        eval(['rts = n' num2str(task) 'subjrt(subj,:);'])
        display = data.task_displayed(subj,:);
        init = NaN(1,default_length);
        init(~isnan(rts)) = rts(~isnan(rts));
        matrix = [matrix; init];
    end
    big_matrix = [big_matrix; matrix]; % task agnostic measure of rt by task iteration
    %errorbar(nanmean(matrix),nanstd(matrix)/sqrt(n),'Color',taskcolors(task+1,:),'LineWidth',1.5)
    %hold on
end
%title('BDM RT by iteration by task displayed')

[h,p] = ttest2(big_matrix(:,1),big_matrix(:,5)); %choosing iter 1 versus iter 5
disp(['t-test iter 1 rts vs iter 5 rts p = ' num2str(p)])
[h,p] = ttest2(data.BDMrt(:,1),data.BDMrt(:,31));
disp(['t-test block 1 rts vs block 31 rts p = ' num2str(p)])

% correlate block # (out of 32) with perf measures
blockbyperf_mat = [reshape(data.perf',n*32,1) repmat([1:32]',n,1)];
blockbyperf_mat(isnan(blockbyperf_mat),:) = [];
[r,p] = corr(blockbyperf_mat(:,1),blockbyperf_mat(:,2));
disp(['Relationship round number/accuracy p = ' num2str(p)])
clear blockbyperf_mat;
blockbyperf_mat = [reshape(data.meanRTs',n*32,1) repmat([1:32]',n,1)];
blockbyperf_mat(isnan(blockbyperf_mat),:) = [];
[r,p] = corr(blockbyperf_mat(:,1),blockbyperf_mat(:,2));
disp(['Relationship round number/mean RTs p = ' num2str(p)])

% test effect of iteration on BDM rating, run 2-way anova on iteration and
% task
n1subjvalue = n1subjvalue(:,1:10); n2subjvalue = n2subjvalue(:,1:10); n3subjvalue = n3subjvalue(:,1:10);
vector = [n1subjvalue(:);n2subjvalue(:);n3subjvalue(:)];
taskidentity = [ones(length(n1subjvalue(:)),1);2*ones(length(n2subjvalue(:)),1);3*ones(length(n3subjvalue(:)),1)];
iternum = repmat([1 2 3 4 5 6 7 8 9 10],n*3,1);
[~,~,stats] = anovan(vector,[taskidentity iternum(:)],'display','off');

% run the same ANOVA on task accuracy
n1subjlearning = n1subjlearning(:,1:10); n2subjlearning = n2subjlearning(:,1:10); n3subjlearning = n3subjlearning(:,1:10);
vector = [n1subjlearning(:);n2subjlearning(:);n3subjlearning(:)];
taskidentity = [ones(length(n1subjlearning(:)),1);2*ones(length(n2subjlearning(:)),1);3*ones(length(n3subjlearning(:)),1)];
iternum = repmat([1 2 3 4 5 6 7 8 9 10],n*3,1);
[~,~,stats] = anovan(vector,[taskidentity iternum(:)],'display','off');

%% NFC and SAPS stuff
measures = [data.NFC data.SAPS];
measures(sum(isnan(measures),2)>1,:) = [];
agemeasures = [data.NFC data.SAPS data.age];
ages = agemeasures(sum(isnan(agemeasures),2)==0,:);
disp(['Mean(std) NFC = ' num2str(nanmean(data.NFC)) '(' num2str(nanstd(data.NFC)) '); mean(std) SAPS = ' num2str(nanmean(data.SAPS)) '(' num2str(nanstd(data.SAPS)) ')'])
disp(['Missing ' num2str(sum(isnan(data.NFC))) ' NFC; missing ' num2str(sum(isnan(data.SAPS))) ' SAPS.'])
[r,p] = corr(measures(:,1),measures(:,2));
disp(['Corr NFC/SAPS r = ' num2str(r) '; p = ' num2str(p)])
[r,p] = corr(ages(:,1),ages(:,3));
disp(['Corr NFC/age r = ' num2str(r) '; p = ' num2str(p)])
[r,p] = corr(ages(:,2),ages(:,3));
disp(['Corr SAPS/age r = ' num2str(r) '; p = ' num2str(p)])

% PLOT DISTRIBUTIONS
figure
subplot(2,2,1)
histogram(data.NFC)
xlabel('Score'); ylabel('# subjects')
title('Distribution of NFC')
ax = gca; ax.FontSize = 12; 
subplot(2,2,2)
histogram(data.SAPS)
xlabel('Score'); ylabel('# subjects')
title('Distribution of SAPS')
fig = gcf; fig.Color = 'w';
ax = gca; ax.FontSize = 12; 

% % LINEAR RELATIONSHIPS OF NFC/SAPS and ACCURACY/MEAN RTs/DIFF
% RATINGS/BDMs
names = {'NFC','SAPS'};
measures = [data.NFC data.SAPS];
trim = [measures tasks_overall tasks_rts data.taskratings];
trim(sum(isnan(trim),2)>0,:) = [];
BDMs = data.values(sum(isnan(trim),2)==0,:); taskBDMs = [];
task_by_round{1} = n1(sum(isnan(trim),2)==0,:); task_by_round{2} = ndetect(sum(isnan(trim),2)==0,:); task_by_round{3} = n2(sum(isnan(trim),2)==0,:);
for subj = 1:size(BDMs,1)
    taskBDMs = [taskBDMs; nanmean(BDMs(subj,task_by_round{1}(subj,:))) nanmean(BDMs(subj,task_by_round{2}(subj,:))) nanmean(BDMs(subj,task_by_round{3}(subj,:)))];
end

for m = 1:2   
    % ACCURACY CORRs
    [r,p] = corr(trim(:,m),nanmean(trim(:,3:6),2));
    disp(['Relationship of ' names{m} ' and overall acc: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,3));
    disp(['Relationship of ' names{m} ' and 1-detect acc: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,4));
    disp(['Relationship of ' names{m} ' and 1-back acc: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,5));
    disp(['Relationship of ' names{m} ' and 3-detect acc: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,6));
    disp(['Relationship of ' names{m} ' and 2-back acc: r = ' num2str(r) '; p = ' num2str(p)])
    disp('. . .')
    % RT CORRs
    [r,p] = corr(trim(:,m),nanmean(trim(:,7:10),2));
    disp(['Relationship of ' names{m} ' and overall RT: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,7));
    disp(['Relationship of ' names{m} ' and 1-detect RT: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,8));
    disp(['Relationship of ' names{m} ' and 1-back RT: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,9));
    disp(['Relationship of ' names{m} ' and 3-detect RT: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,10));
    disp(['Relationship of ' names{m} ' and 2-back RT: r = ' num2str(r) '; p = ' num2str(p)])
    disp('. . .')
    % DIFF RATINGS
    [r,p] = corr(trim(:,m),nanmean(trim(:,11:14),2));
    disp(['Relationship of ' names{m} ' and overall rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,11));
    disp(['Relationship of ' names{m} ' and 1-detect rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,12));
    disp(['Relationship of ' names{m} ' and 1-back rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,13));
    disp(['Relationship of ' names{m} ' and 3-detect rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),trim(:,14));
    disp(['Relationship of ' names{m} ' and 2-back rating: r = ' num2str(r) '; p = ' num2str(p)])
    disp('. . .')
    % MEAN BDMs
    [r,p] = corr(trim(:,m),nanmean(BDMs,2));
    disp(['Relationship of ' names{m} ' and mean BDM rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),taskBDMs(:,1));
    disp(['Relationship of ' names{m} ' and 1-back BDM rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),taskBDMs(:,2));
    disp(['Relationship of ' names{m} ' and 3-detect BDM rating: r = ' num2str(r) '; p = ' num2str(p)])
    [r,p] = corr(trim(:,m),taskBDMs(:,3));
    disp(['Relationship of ' names{m} ' and 2-back BDM rating: r = ' num2str(r) '; p = ' num2str(p)])
    disp('. . .')
end

%% % GROUP ANALYSIS OF BDMS by NFC AND SAPS GROUPS % % 
measures = [data.NFC data.SAPS];
labels = {'NFC','SAPS'};
for measure = 1:2
    split = [];
    split = tertileSplit(measures(:,measure));
    disp([labels{measure} ' group Ns: ' num2str([sum(split==1) sum(split==2) sum(split==3)])])
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
    title(['Mean Fair Wage by ' labels{measure} ' group'])
    ylabel('Wage')
    xlabel('Task')
    xticklabels(tasklabels(2:end))
    xticks([1:3])
    ax = gca; ax.FontSize = 12;

    tasklist = [repmat(1,n,1); repmat(2,n,1); repmat(3,n,1)];
    ratings = [nanmean(n1subjvalue,2); nanmean(n3subjvalue,2); nanmean(n2subjvalue,2)];
    matrix = [ratings repmat(split,3,1) tasklist];
    [~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{labels{measure},'task'},'display','off');
    
end

% % regress mean BDM values by SAPS and NFC, and quadratic terms on them,
% too % %
Y = nanmean(data.values,2);
X = [ones(n,1) data.NFC data.SAPS data.NFC.^2 data.SAPS.^2];
invalid = sum(isnan(X),2)>0; Y(invalid,:) = []; X(invalid,:) = [];

[betas,BINT,~,~,stats] = regress(Y,X); 
% get betas for quadratic term
predicted = X*betas;
distance = predicted-Y;
MSE = distance'*distance; %squared distance
%disp(['betas on intercept, NFC, SAPS, NFC^2, and SAPS^2: ' num2str(betas') ', p = ' num2str(stats(3)), ', MSE = ' num2str(MSE)])
% The only real term here is the intercept

%Examine individual NFC group differences post-ANOVA main effect
split = tertileSplit(data.NFC);
lowNFCvalues = data.values(split==1,:);midNFCvalues = data.values(split==2,:);highNFCvalues = data.values(split==3,:);
low = nanmean(lowNFCvalues,2); mid = nanmean(midNFCvalues,2); high = nanmean(highNFCvalues,2);
[h,p] = ttest2(low,mid);
disp(['Mid versus low NFC fair wage ratings: p = ' num2str(p)])
[h,p] = ttest2(high,low);
disp(['High versus low NFC fair wage ratings: p = ' num2str(p)])
[h,p] = ttest2(mid,high);
disp(['Mid versus high NFC fair wage ratings: p = ' num2str(p)])

%Do individual differences relate to baseline executive function?
%0-detect, 3-detect, 1-back, 2-back (in that order in practice)
meanpracaccs = nanmean(data.practiceacc,2);
[r,p] = corr(meanpracaccs(~isnan(data.NFC)),data.NFC(~isnan(data.NFC)));
disp(['Relationship of overall practice accuracy & NFC r = ' num2str(r) ', p = ' num2str(p)])
[r,p] = corr(meanpracaccs(~isnan(data.SAPS)),data.SAPS(~isnan(data.SAPS)));
disp(['Relationship of overall practice accuracy & SAPS r = ' num2str(r) ', p = ' num2str(p)])
practiceaccs = data.practiceacc(:,1:4); 
matrix = [practiceaccs(:) [ones(n,1); repmat(2,n,1); repmat(3,n,1); repmat(4,n,1)] repmat(split,4,1)];
[~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{'task','NFC'},'display','off');

%From that main effect of NFC on practice accuracy, run some post-hoc tests
%of NFC groups & practice accuracy
[h,p] = ttest2(practiceaccs(split==1,:),practiceaccs(split==2,:))
disp('t-test low vs mid NFC in practice rounds 1-4')
[h,p] = ttest2(practiceaccs(split==2,:),practiceaccs(split==3,:))
disp('t-test mid vs high NFC in practice rounds 1-4')
[h,p] = ttest2(practiceaccs(split==1,:),practiceaccs(split==3,:))
disp('t-test low vs high NFC in practice rounds 1-4')

%Does this difference persist into the actual experiment?
matrix = [tasks_overall(:) [ones(n,1); repmat(2,n,1); repmat(3,n,1); repmat(4,n,1)] repmat(split,4,1)];
[~,~,stats] = anovan(matrix(:,1),{matrix(:,2),matrix(:,3)},'model','interaction','varnames',{'task','NFC'},'display','off');

%Does this difference account for NFC effect on fair wage?
matrix = [nanmean(n1subjvalue,2) nanmean(n2subjvalue,2) nanmean(n3subjvalue,2) data.NFC];
matrix(sum(isnan(matrix),2)>0,:) = [];
[r,p] = corr(matrix)

%% modeling results (parameters) versus behavioral results (NFC and SAPS)
addpath('HBI/')
%modelfits = load('data/modelfits_2_models_13-Jan-2021.mat'); %miss c, response c, lure c - most recoverable
modelfits = load('HBI/HBI_modelStruct.mat');
no_fits = load('simdata/toanalyze.mat','trim'); no_fits = no_fits.trim;
% remove subjects who weren't fit to models
% grab best model
best_model = modelfits.modelStruct.best_model;
responsibility = modelfits.modelStruct.best_model.overallfit.responsibility;
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
    disp([names{measure} ' vs ' paramnames(p) ' r = ' num2str(r) ', p = ' num2str(pval)])
end
fig = gcf; fig.Color = 'w';


figure
for p = 1:nparams
    
    % run group-level t-tests
    [h,pval1] = ttest2(params(split==1,p),params(split==2,p));
    [h,pval2] = ttest2(params(split==3,p),params(split==2,p));
    [h,pval3] = ttest2(params(split==1,p),params(split==3,p));
    ps = [NaN pval1 pval3; pval1 NaN pval2; pval3 pval2 NaN];
    
    subplot(3,2,p)
    superbar([nanmean(params(split==1,p)) nanmean(params(split==2,p)) nanmean(params(split==3,p))], ...
                'E',[nanstd(params(split==1,p)) nanstd(params(split==2,p)) nanstd(params(split==3,p))]./sqrt(n), ...
                'P',ps,'BarFaceColor',colors,'PStarShowNS',false);
    title(paramnames(p))
    xticks([1:3])
    xticklabels({['Low ' names{measure}],['Mid ' names{measure}],['High ' names{measure}]})
    
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
[h,p] = ttest2(measures([assignments==models_at_play(1)]'&sum(~isnan(measures),2)==2,1),measures([assignments==models_at_play(2)]'&sum(~isnan(measures),2)==2,1))
disp('NFC across models 1 and 2')
[h,p] = ttest2(measures([assignments==models_at_play(1)]'&sum(~isnan(measures),2)==2,1),measures([assignments==models_at_play(1)]'&sum(~isnan(measures),2)==2,2))
disp('SAPS across models 1 and 2')
[h,p] = ttest2(measures([assignments==models_at_play(2)]'&sum(~isnan(measures),2)==2,1),measures([assignments==models_at_play(3)]'&sum(~isnan(measures),2)==2,1))
disp('NFC across models 2 and 3')
[h,p] = ttest2(measures([assignments==models_at_play(2)]'&sum(~isnan(measures),2)==2,1),measures([assignments==models_at_play(3)]'&sum(~isnan(measures),2)==2,2))
disp('SAPS across models 2 and 3')

% Assign subjects to their model, plot parameters that way
all_params = best_model.overallfit.parameters; 
for measure = 1:2
    if measure == 1
        split = tertileSplit(data.NFC); colors = NFCcolors;
        % Also, run multiple linear regressions
        X = [ones(n,1) data.NFC data.NFC.^2];
    elseif measure == 2
        split = tertileSplit(data.SAPS); colors = SAPScolors;
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
        paramnames = strrep(paramnames,'u ','update '); strrep(paramnames,'main ','maintenance ');
        costs = find(contains(paramnames,'cost'));
        values = all_params{modelnum}(:,costs);
        for c = 1:size(costs,2)
            count = count+1;
            
            % run group-based & continuous stats
            [h,pval1] = ttest2(values(split==1,c),values(split==2,c));
            [h,pval2] = ttest2(values(split==3,c),values(split==2,c));
            [h,pval3] = ttest2(values(split==1,c),values(split==3,c));
            ps = [NaN pval1 pval3; pval1 NaN pval2; pval3 pval2 NaN];
            
            Y = nanmean(values,2); Y(invalid,:) = [];
            [betas,BINV,~,~,stats] = regress(Y,X);
            % get betas for quadratic term
            predicted = X*betas;
            distance = predicted-Y; MSE = distance'*distance; %squared distance
            disp([names{measure} ' betas on ' paramnames{costs(c)} ', linear & quadratic: ' num2str(betas')])
            
            % plot effect of measure group on parameter values
            subplot(4,3,count)
            superbar([nanmean(values(split==1,c)) nanmean(values(split==2,c)) nanmean(values(split==3,c))], ...
                'E',[nanstd(values(split==1,c)) nanstd(values(split==2,c)) nanstd(values(split==3,c))]./sqrt(n), ...
                'P',ps,'BarFaceColor',colors,'PStarShowNS',false);
            title(paramnames{costs(c)})
            xticks([1:3])
            xticklabels({['Low ' names{measure}],['Mid ' names{measure}],['High ' names{measure}]})
            
        end %of cycling over each cost
        
        plotted = [plotted paramnames(costs)];
    end %of cycling over each model
    
end %of cycling over SAPS & NFC

