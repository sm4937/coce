%% Load in subject data for use in simulations
% load in actual subject trials for this
clear all 
preloadflag = false;
if preloadflag
    files = {'29.04.2020.mat','07.05.2020.mat','12.05.2020.mat','13.05.2020.mat','18.05.2020.mat'};
    load('simdata\fullsubjnumbers.mat') %grab usable subjects
    tosim = [];
    for file = 1:length(files)
        load(['data\' files{file}])
        for row = 1:length(list)
            subj = list(row);
            data = Untitled(Untitled.subjnum==subj,:);
            if height(data)>0 %that subj might not be in that file
            data = sortrows(data,find(ismember(data.Properties.VariableNames,'trial_index'))); %a little jumbled up still
            oneperson = table;
            tasknum = double(string(data.tasknum));
            idx = ~isnan(tasknum); %clean this up a little before sim stuff
            oneperson.block = tasknum(idx);
            oneperson.subj = repmat(row,height(oneperson),1);
            oneperson.n = data.n(idx); oneperson.n(1:10) = 0; %the first 10 trials are 0-back practice
            oneperson.detect = data.detect(idx)==categorical(cellstr('1'));
            oneperson.nback = data.nback(idx)==categorical(cellstr('1'));
            values = double(string(data.value));
            display = double(string(data.task_displayed)); %the last four values here are for difficulty ratings
            display(end-15:end) = []; %so trim them
            oneperson.BDM = NaN(height(oneperson),1); oneperson.display = NaN(height(oneperson),1);
            d = sum(isnan(oneperson.n))-sum(~isnan(values));
            oneperson.BDM(isnan(oneperson.n)) = [values(~isnan(values)); NaN(d,1)];
            d = sum(isnan(oneperson.n))-sum(~isnan(display));
            oneperson.display(isnan(oneperson.n)) = [display(~isnan(display)); NaN(d,1)];
            oneperson.display(oneperson.display==7) = 3;
            oneperson.BDM = (oneperson.BDM-1)*25; %from 0 to 100, to clarify math stuff
            oneperson.task = NaN(height(oneperson),1);
            tasklist = data.task(idx);
            oneperson.task(tasklist==categorical(cellstr('detection'))) = 0;
            oneperson.task(tasklist==categorical(cellstr('n-back'))&oneperson.n==1) = 1;
            oneperson.task(tasklist==categorical(cellstr('n-back'))&oneperson.n==2) = 2;
            oneperson.task(tasklist==categorical(cellstr('ndetection'))) = 3;
            oneperson.correct = false(height(oneperson),1);
            oneperson.correct = data.correct(idx)==categorical(cellstr('1'));
            oneperson.stimnum = data.stimnum(idx);
            oneperson.nmatches = double(string(data.nmatches(idx)));
            oneperson.nmatches(isnan(oneperson.nmatches)) = 0; %replace NaN's with 0
            oneperson.keypress = data.key_press(idx);
            tosim = [tosim; oneperson];
            end
        end
    end %end of loading/formatting
    save('simdata\simdata_n30.mat','tosim')
else
    load('simdata\simdata_n30.mat')
end

n_tasks = length(unique(tosim.task(~isnan(tosim.task))));
tasklabels = {'1-detect','1-back','2-back','3-detect'};
%% Pull possible costs from subjects' behavior and task presented to them
% Goes through the tasks round by round to tally the costs incurred by each
% subject on each round. Does it slightly according to task rules and slightly
% according to an internal process which is common across tasks.

practiceblocks = tosim(tosim.block==-1,:);
% Will become relevant when I get practice stimnums out of new subjects.
% But for now, irrelevant.

components = []; %bookkeeping 
noise = 0; %not really doing this in a biologically plausible way, just setting noise to 0
for subj = 1:length(unique(tosim.subj))
    oneperson = tosim(tosim.subj==subj,:);
    for block = 0:max(oneperson.block) 
        blockdata = oneperson(oneperson.block==block,:);
        nupdates = 0; nlures = 0; %initialize cost variables
        maintained = NaN(1,2); %keep track of storage over time, to get a sense of maintained info over time
        task = blockdata.task(2)+1;
        n = blockdata.n(2); matches = blockdata.nmatches(2);
        storage = NaN(1,n); %empty storage for n-back
        for trial = 2:height(blockdata)
            stim = blockdata.stimnum(trial);
            if n == 0 %do detection, pretty simple
                maintained(trial,:) = 0;
            else %n-backs and n-detects
                if task==3 %only in 2-back can there be WM lures inside buffer
                    if storage(2)==stim
                        nlures = nlures + 1; %not a real 2-back, but a lure trial
                    end
                end
                if sum(stim == storage)==n %full storage remaining the same
                    % do nothing, no cost, no update
                else %update storage with some noise
                    storage(1) = []; 
                    die = rand(); %get a random number
                    if die < noise %noise param wins
                        storage(end+1) = NaN;
                    else %noise param doesn't win
                        storage(end+1) = stim; %clear out old info, add in new info, shift stuff over
                    end
                    nupdates = nupdates + 1; %cost of update/gating in/out
                end
                maintained(trial,:) = storage;
            end
        end
        % grab noisiness for whole block instead of trial
        noisiness = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==false)/max([sum(blockdata.detect|blockdata.nback) 1]);
        nmatches = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==true);
        nmisses = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==false);
        nresponses = sum(~isnan(blockdata.keypress(2:end)));
        maintained(1,:) = []; maintained = nanmean(sum(maintained>0,2));
        
        components = [components; subj nmatches maintained nupdates nmisses blockdata.display(1)+1 task blockdata.BDM(1) noisiness nresponses nlures];
    end %end of block by block loop
end

toanalyze = table;
toanalyze.subj = components(:,1); toanalyze.nmatches = components(:,2);
toanalyze.maintained = components(:,3); toanalyze.nupdates = components(:,4);
toanalyze.nmisses = components(:,5); toanalyze.display = components(:,6); 
toanalyze.task = components(:,7); toanalyze.BDM = components(:,8); 
toanalyze.noisiness = components(:,9); toanalyze.nresponses = components(:,10);
toanalyze.nlures = components(:,11);
%save('simdata/toanalyze.mat','toanalyze')

%% Check out individual differences in measures
nsubjs = length(unique(tosim.subj));

figure; subplot(2,2,2)
for misses = 0:4
    for subj = 1:nsubjs
        nmisses(subj) = max(toanalyze.nmisses(toanalyze.subj==subj));
        BDMs = toanalyze.BDM(toanalyze.subj==subj,:);
        nmissesall(:,subj) = [toanalyze.nmisses(toanalyze.subj==subj); nan(32-sum(toanalyze.subj==subj),1)];
        misseffect(misses+1,:,subj) = NaN(1,32); %initialize empty
        idx = find(nmissesall(:,subj)==misses)+1; idx(idx>length(BDMs)) = []; %trim out edges
        misseffect(misses+1,idx,subj) = BDMs(idx); %catalog BDM effects from misses
    end
    scatter(1:nsubjs,sum(nmissesall==misses),'Filled')
    hold on
end
legend({'0','1','2','3','4'})
title('All Misses per Subject')
ylabel('Frequency'); xlabel('Subject')
subplot(2,2,3)
histogram(nmisses)
title('Max number of misses for each subject')
subplot(2,2,1)
histogram(toanalyze.nmisses)
title('All misses values all subjects')
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
title('Mean BDM following each # of misses')
%mean wage for each subject following n misses
totest = [squeeze(nanmean(misseffect(1,:,:),2)) zeros(nsubjs,1); squeeze(nanmean(misseffect(2,:,:),2)) ones(nsubjs,1); ...
    squeeze(nanmean(misseffect(3,:,:),2)) 2*ones(nsubjs,1); squeeze(nanmean(misseffect(4,:,:),2)) 3*ones(nsubjs,1); ...
    squeeze(nanmean(misseffect(5,:,:),2)) 4*ones(nsubjs,1)];
[~,~,stats] = anova1(totest(:,1),totest(:,2));

%% Simulate some stuff to see if behavior can be captured by these costs
nsubjs = length(unique(toanalyze.subj))-1;
subjcolors = rand(nsubjs*2);subjcolors(:,4:end) = []; %delete unnecessary columns
taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
%plot one subj's learning block by block, to see how the dynamics of this
%go
realvssim = [toanalyze.subj toanalyze.display toanalyze.BDM components(:,8)];
figure
subplot(3,2,1)
for subj = 1:2
    toplot = realvssim(realvssim(:,1)==subj,:);
    scatter(1:32,toplot(:,3),[],subjcolors(subj,:),'Filled') 
    hold on
    scatter(1:32,toplot(:,4),[],subjcolors(subj,:)) 
end
legend({'Real','Sim'})
xlabel('Block')
title('Simulated vs real BDMs')
fig = gcf; fig.Color = 'w';

subplot(3,2,2)
for task = 1:length(unique(toanalyze.display(~isnan(toanalyze.display)))) %should be three rated tasks
    for subj = 1:2
        toplot = realvssim(realvssim(:,1)==subj&realvssim(:,2)==task+1,:);
        scatter(find(toplot(:,2)==task+1),toplot(:,3),[],taskcolors(task,:),'Filled') 
        hold on
        scatter(find(toplot(:,2)==task+1),toplot(:,4),[],taskcolors(task,:)) 
    end
end
legend({'Real','Sim'})
xlabel('Block')
title('Simulated vs real BDMs - task specific')
fig = gcf; fig.Color = 'w';

subplot(3,2,3)
%plot group means instead of individual subjects
errorbar(nanmean(BDMrows),nanstd(BDMrows)/sqrt(nsubjs),'LineWidth',1.5)
hold on
errorbar(nanmean(simBDMrows),nanstd(simBDMrows)/sqrt(nsubjs),'LineWidth',1.5)
legend({'Real','Sim'})
title('Group mean BDMs by block')

ofinterest = toanalyze;
ofinterest(ofinterest.subj==8,:) = []; %delete subj 8 because they're problematic
nsubjs = length(unique(ofinterest.subj));
tasks = reshape(ofinterest.display,nsubjs,32);

blankrow = NaN(1,32); rts = []; sts = [];
for task = 2:4
    subplot(3,2,2+task)
    subjs = unique(ofinterest.subj);
    for row = 1:nsubjs
        subj = subjs(row,:);
        realtaskspecific = blankrow;
        realtaskspecific(tasks(row,:)==task) = BDMrows(subj,tasks(row,:)==task);
        rts = [rts; realtaskspecific];
        simtaskspecific = blankrow;
        simtaskspecific(tasks(row,:)==task) = simBDMrows(subj,tasks(row,:)==task);
        sts = [sts; simtaskspecific];
    end
    errorbar(nanmean(rts),nanstd(rts)/sqrt(nsubjs),'Color',taskcolors(task-1,:),'LineWidth',1.5)
    hold on
    errorbar(nanmean(sts),nanstd(sts)/sqrt(nsubjs),'*-','Color',taskcolors(task-1,:),'LineWidth',1.5)
    title(['Sim vs. real bdm in ' tasklabels(task)])
    xlabel('block')
    legend({'Real','Sim'})
    ylabel('average BDM value')
end

%% Run linear mixed effect models on simulations
%mean center the appropriate scores
% for task = 2:4
tbl = toanalyze;
tbl.nmatches = tbl.nmatches - nanmean(tbl.nmatches);
tbl.nmisses = tbl.nmisses - nanmean(tbl.nmisses); 
tbl.maintained = tbl.maintained - nanmean(tbl.maintained);
tbl.nupdates = tbl.nupdates - nanmean(tbl.nupdates);
tbl.newBDM = 1+(tbl.newBDM/25);

random_two = fitlme(tbl,'newBDM ~ 1 + nmatches + nmisses + maintained + nupdates + ( 1 + nmatches | subj) + ( 1 + nmisses | subj)');
random_three = fitlme(tbl,'newBDM ~ 1 + nmatches + nmisses + maintained + nupdates + ( 1 + nmatches | subj) + ( 1 + nmisses | subj) + ( 1 + nupdates | subj)');
random_new = fitlme(tbl,'newBDM ~ 1 + nmatches + nmisses + maintained + nupdates + ( 1 + nmatches | subj) + ( 1 + nupdates | subj)');

AICs = [-random_two.LogLikelihood.*2+(2.*numel(random_two.covarianceParameters{1})), ...
    -random_new.LogLikelihood.*2+(2.*numel(random_new.covarianceParameters{1})), ...
    -random_three.LogLikelihood.*2+(2.*numel(random_three.covarianceParameters{1}))]; %, ...
%     -no_random_test.LogLikelihood.*2+(2.*0), ...
%     -random_int_test.LogLikelihood.*2+(2.*numel(random_int_test.covarianceParameters{1})), ...
%     -random_both_test.LogLikelihood.*2+(2.*numel(random_both_test.covarianceParameters{1}))];\

model_struct.m1 = random_two;
model_struct.m2 = random_new;
model_struct.m3 = random_three;

n = length(unique(toanalyze.subj))-1;

figure
bar(AICs)
hold on
[score,which] = min(AICs);
%which = 3; %the second and third models are equivalent now...
plot(which,score,'*k','LineWidth',1.5)
title(['AIC by model, fit to all task data'])
ylabel('AIC')
xlabel('Model')
xticklabels({'Nmatches, Nmisses','Nmatches, Nupdates','Nmatches, Nupdates, Nmisses'})
xtickangle(45)

% look into values for minimum model

eval(['winner = model_struct.m' num2str(which)])

[~,~,rEffects] = randomEffects(winner);
estimates_inter = fixedEffects(winner);
estimates_inter = [estimates_inter(1) estimates_inter(2) estimates_inter(3) estimates_inter(5)];

figure
subplot(6,6,6)
plot(0:0.05:4,normpdf(0:0.05:4,estimates_inter(1),std(rEffects.Estimate(1:2:n*2))))
hold on
%scatter(rEffects.Estimate(1:2:n*2)+estimates_inter(1),normpdf(rEffects.Estimate(1:2:n*2)+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate(1:2:n*2))))
%the first line is plotting the individual values on the exact line
yvalue = max(normpdf(rEffects.Estimate(1:2:n*2)+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate(1:2:n*2))))+1;
scatter(rEffects.Estimate(1:2:n*2)+estimates_inter(1),yvalue*ones(n,1))
title('Fit Match Intercepts')

subplot(6,6,11)
plot(-1:0.05:1,normpdf(-1:0.05:1,estimates_inter(2),std(rEffects.Estimate(2:2:n*2))))
hold on
%scatter(rEffects.Estimate(2:2:n*2)+estimates_inter(2),normpdf(rEffects.Estimate(2:2:n*2)+estimates_inter(2),estimates_inter(2),std(rEffects.Estimate(2:2:n*2))))
yvalue = max(normpdf(rEffects.Estimate(2:2:n*2)+estimates_inter(2),estimates_inter(2),std(rEffects.Estimate(2:2:n*2))))+2;
scatter(rEffects.Estimate(2:2:n*2)+estimates_inter(2),yvalue*ones(n,1))
title('Fit Match Slopes')

%match intercepts versus slopes
subplot(6,6,12)
scatter(rEffects.Estimate(2:2:n*2)+estimates_inter(2),rEffects.Estimate(1:2:n*2)+estimates_inter(1),'Filled')
[r,p] = corr(rEffects.Estimate(2:2:n*2)+estimates_inter(2),rEffects.Estimate(1:2:n*2)+estimates_inter(1));
if p < 0.05
    lsline
end
%xlabel('Match slopes')
%ylabel('Match intercepts')
%match slopes, intercepts, correlation

subplot(6,6,16)
plot(0:0.05:4,normpdf(0:0.05:4,estimates_inter(1),std(rEffects.Estimate((n*2)+1:2:(n*4)))))
hold on
%scatter(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),normpdf(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate((n*2)+1:2:(n*4)))))
yvalue = max(normpdf(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate((n*2)+1:2:(n*4)))))+1;
scatter(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),yvalue*ones(n,1))
title('Fit Miss Intercepts')

subplot(6,6,17)
scatter(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),rEffects.Estimate(2:2:n*2)+estimates_inter(2),'Filled')
%xlabel('Miss intercepts')
%ylabel('Match slopes')
[r,p] = corr(rEffects.Estimate(2:2:n*2)+estimates_inter(2),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,18)
scatter(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),rEffects.Estimate(1:2:n*2)+estimates_inter(1),'Filled')
%xlabel('Miss intercepts')
%ylabel('Match intercepts')
[r,p] = corr(rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),rEffects.Estimate(1:2:n*2)+estimates_inter(1));
if p < 0.05
    lsline
end


subplot(6,6,21)
plot(-1:0.05:1,normpdf(-1:0.05:1,estimates_inter(3),std(rEffects.Estimate((n*2)+2:2:(n*4)))))
hold on
%scatter(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),normpdf(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),estimates_inter(3),std(rEffects.Estimate((n*2)+2:2:(n*4)))))
yvalue = max(normpdf(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),estimates_inter(3),std(rEffects.Estimate((n*2)+2:2:(n*4)))))+2;
scatter(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),yvalue*ones(n,1))
title('Fit Miss Slopes')

subplot(6,6,22)
scatter(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),'Filled')
%xlabel('Miss slopes')
%ylabel('Miss intercepts')
[r,p] = corr(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1));
if p < 0.05
    lsline
end
%miss slope, miss intercepts, correlation of both

subplot(6,6,23)
scatter(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate(2:2:n*2)+estimates_inter(2),'Filled')
%xlabel('Miss slopes')
%ylabel('Match slopes')
[r,p] = corr(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate(2:2:n*2)+estimates_inter(2));
if p < 0.05
    lsline
end

subplot(6,6,24)
scatter(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate(1:2:n*2)+estimates_inter(1),'Filled')
%xlabel('Miss slopes')
%ylabel('Match intercepts')
[r,p] = corr(rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),rEffects.Estimate(1:2:n*2)+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,26)
plot(0:0.05:4,normpdf(0:0.05:4,estimates_inter(1),std(rEffects.Estimate((n*4)+1:2:(n*6)))))
hold on
%scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),normpdf(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate((n*4)+1:2:(n*6)))))
yvalue = max(normpdf(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),estimates_inter(1),std(rEffects.Estimate((n*4)+1:2:(n*6)))))+1;
scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),yvalue*ones(n,1))
title('Fit Update Intercepts')

subplot(6,6,27)
scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),'Filled')
%xlabel('Update Intercepts')
%ylabel('Miss slopes')
[r,p] = corr(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3));
if p < 0.05
    lsline
end

subplot(6,6,28)
scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),'Filled')
%xlabel('Update intercepts')
%ylabel('Miss intercepts')
[r,p] = corr(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,29)
scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate(2:2:n*2)+estimates_inter(2),'Filled')
%xlabel('Update intercepts')
%ylabel('Match slopes')
[r,p] = corr(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate(2:2:n*2)+estimates_inter(2));
if p < 0.05
    lsline
end

subplot(6,6,30)
scatter(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate(1:2:n*2)+estimates_inter(1),'Filled')
%xlabel('Update intercepts')
%ylabel('Match intercepts')
[r,p] = corr(rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),rEffects.Estimate(1:2:n*2)+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,31)
plot(-1:0.05:1,normpdf(-1:0.05:1,estimates_inter(4),std(rEffects.Estimate((n*4)+2:2:(n*6)))))
hold on
%scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),normpdf(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),estimates_inter(4),std(rEffects.Estimate((n*4)+2:2:(n*6)))))
yvalue = max(normpdf(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),estimates_inter(4),std(rEffects.Estimate((n*4)+2:2:(n*6)))))+2;
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),yvalue*ones(n,1))
title('Fit Update Slopes')
fig = gcf; fig.Color = 'w';

subplot(6,6,32)
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1),'Filled')
%xlabel('Update Slopes')
%ylabel('Update Intercepts')
[r,p] = corr(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*4)+1:2:(n*6))+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,33)
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3),'Filled')
%xlabel('Update slopes')
%ylabel('Miss slopes')
[r,p] = corr(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*2)+2:2:(n*4))+estimates_inter(3));
if p < 0.05
    lsline
end

subplot(6,6,34)
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1),'Filled')
%xlabel('Update slopes')
%ylabel('Miss intercepts')
[r,p] = corr(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate((n*2)+1:2:(n*4))+estimates_inter(1));
if p < 0.05
    lsline
end

subplot(6,6,35)
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate(2:2:n*2)+estimates_inter(2),'Filled')
%xlabel('Update slopes')
%ylabel('Match slopes')
[r,p] = corr(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate(2:2:n*2)+estimates_inter(2));
if p < 0.05
    lsline
end

subplot(6,6,36)
scatter(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate(1:2:n*2)+estimates_inter(1),'Filled')
%xlabel('Update slopes')
%ylabel('Match intercepts')
[r,p] = corr(rEffects.Estimate((n*4)+2:2:(n*6))+estimates_inter(4),rEffects.Estimate(1:2:n*2)+estimates_inter(1));
if p < 0.05
    lsline
end

% end

%% Run simple simulations for task characteristics, individual characteristics, simple models

%% Simulate slopes based on NFC

alpha = 0.50; %baseline slope
NFC = 1; %fake, as if NFC were between 0 and 1
B = -10:10;
S = alpha*(1./(1+exp(-B.*NFC)));
figure
scatter(B,S)
ylabel('Slope')
title(['Formula alpha*sigmoid(-B*NFC), NFC: ' num2str(NFC)])
xlabel('Beta Value')
xticks(B)
xlim([B(1) B(end)])
xticklabels(B)
ylim([-1 1])

NFCs = data.NFC;
maxNFC = 5;
NFCs = NFCs./5;

Beta = 3;
alpha = 3;
slopes = alpha*(1-(1./(1+exp(-Beta.*NFCs))));
%slopes = alpha + (Beta.*NFCs);
nswitches = 6:15;
Bs = rand(n,1)*3;
y = slopes.*nswitches + Bs;

figure
subplot(2,1,2)
toplot = 5:10;
for subj = toplot
    scatter(nswitches,y(subj,:),'o','Filled')
    hold on
end
xlabel('N Switches')
ylabel('BDM value')
title('Simulated effect of nswitches according to NFC')
labels = [repmat('NFC: ',length(toplot),1) num2str(NFCs(toplot))];
legend({labels})
ax = gca; 
ax.FontSize = 12;

subplot(2,1,1)
scatter(slopes,NFCs)
xlabel('Slope')
ylabel('NFC score normalized')
title('Slope of subcomponent effect by NFC score')
ax = gca; fig = gcf;
fig.Color = 'w'; ax.FontSize = 12;

%% simulate delay plots by fixed BDM values

BDMs = [0:4:100];
BDMs = 1+(BDMs./25);
nblocks = 32;
sim_list = [];
for p = 1:length(BDMs) %cycle through fixed values of BDMs - mirror consistency of our subjects
    value = BDMs(p);
    tasks = ones(1,nblocks/2); tasks = [tasks 2.*ones(1,nblocks/2)];
    tasks = tasks(randperm(nblocks));
    for block = 1:32
        offer = BDMs(ceil(rand*length(BDMs)));
        if value>offer
            tasks(block) = 0;
        end
    end
    sim_list = [sim_list; tasks value];
end

figure
for subj = 1:length(BDMs)
    value = BDMs(subj);
    row = sim_list(:,end)==value;
    subplot(5,6,subj)
    histogram(sim_list(row,1:end-1))
    xticklabels({'0','1','2'})
    xlabel('Task')
    ylabel('Freq.')
    title(['value = ' num2str(value)])
end

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




