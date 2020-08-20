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
nsubjs = 2; %how many simulated subjects? fixed for now
subjcolors = rand(nsubjs*2);subjcolors(:,4:end) = []; %delete unnecessary columns

%  % SET PARAMETERS %  %
% AKA
% model parameters!
alpha = 0.1; %same for everyone
init = 0.4; %same for everyone
noise = 0.2;

params = struct;
params.alpha = alpha; %learning rate
params.main = 0.2; %cost of maintenance
params.inter = 0.2; %cost of interfering stimuli
params.match = 0.2; %cost of matches
params.miss = 0.2; %cost of misses
params.forget = 0.5;
params.init = init;
params.noise = noise;

parameters{1} = params;

%overwrite the ones that are different for subject 2
params.main = 0.1; %cost of maintenance
params.inter = 0.1; %cost of interfering stimuli
params.match = 0.1; %cost of matches
params.miss = 0.1; %cost of misses
params.forget = 0.6;

parameters{2} = params; %second simulated subject
%% A first-pass cost model, based on task characteristics
% %features: 
% %adapting, changing as task gets learned
% %maintenance has a cost, interfering stimuli have a cost
% 
% %initialize costs of each task
% cost = init*ones(1,n_tasks); %initialize to 25, start low, allow learning during practice round
% real_cost = NaN(1,n_tasks);
% %real_cost(1) = 1;
% 
% %simulate some stuff!!
% figure
% for subj = 1:nsubjs
%     params = parameters{subj};
%     main = params.main; inter = params.inter; match = params.match; miss = params.miss; forget = params.forget;
%     for block = -1:max(tosim.block) %need to separate out practice blocks later
%         blockdata = tosim(tosim.block==block,:);
%         round = block+2;
%         task = blockdata.task(2)+1;
%         real_cost(task) = blockdata.BDM(1); %update actual BDM value on each round
%         n = blockdata.n(2);
%         blockcost = cost(task);
%         blockcost = blockcost + main*n; %add maintenance cost
%         nmatches = nansum(blockdata.nback(blockdata.correct)) + nansum(blockdata.detect(blockdata.correct));
%         nmisses = nansum(blockdata.nback(~blockdata.correct)) + nansum(blockdata.detect(~blockdata.correct));
%         blockcost = blockcost + match*nmatches;
%         blockcost = blockcost - miss*nmisses; %apologetic
%         idxes = find(blockdata.nback|blockdata.detect);
%         distractors = 0;
%         if n>0
%         for i = 1:nmatches
%             idx = idxes(i);
%             within = blockdata(idx-n:idx,:);
%             other = nansum(within.stimnum(1:n-1)~=within.stimnum(n)); %count stimuli that are not the same
%             % for 3detect this should always be 0
%             distractors = distractors + other;
%         end
%         end
%         blockcost = blockcost+(distractors*inter); %inter parameter * distractor
%         cost(task) = cost(task) + alpha*blockcost;
%         not = [1:n_tasks];
%         not(task) = [];
%         cost(not) = cost(not) - forget*(cost(not)-init); %forget tasks that weren't just completed
%     end
%     cost = (cost/25)+1; %go from 1 to 5
%     real_cost = (real_cost/25)+1; %go from 1 to 5
%     
%     scatter(1:length(tasklabels),cost,[],subjcolors(subj,:),'Filled')
%     hold on
%     scatter(1:length(real_cost),real_cost,[],subjcolors(subj+1,:),'Filled')
%     real_cost
%     legend({'Simulated Costs','Real Costs','Simulated Costs','Real Costs'})
%     xticks([1:4])
%     xticklabels(tasklabels)
%     title('Final Costs of Tasks - 1 Sim Subj')
% end

%% Second-pass model, where Matlab does the task itself, and costs come from the operations needed for that
% no threshold yet 
% work through the n-back round by round

practiceblocks = tosim(tosim.block==-1,:);
% how to incorporate practice when I don't have... stimnums... hmm...
% plus the nmatches is always 3 for everything (except ndetect, which has
% 2) 
% leave this for now

components = []; %bookkeeping

%let costs of tasks evolve over time
for subj = 1:length(unique(tosim.subj))
    params = parameters{1};
    %initialize costs of each task
    init = params.init*100;
    cost = init*ones(1,n_tasks); %initialize to 25, start low, allow learning during practice round
    real_cost = NaN(1,n_tasks);
    main = params.main*10; inter = params.inter*10; match = params.match*10; miss = params.miss*10; forget = params.forget; noise = params.noise;
    oneperson = tosim(tosim.subj==subj,:);
    for block = 0:max(oneperson.block) %need to separate out practice blocks later
        blockdata = oneperson(oneperson.block==block,:);
        nmatches = 0; nmisses = 0; nupdates = 0; noutputs = 0; %initialize cost variables
        maintained = NaN(1,2); %keep track of storage over time, to get a sense of maintained info over time
        task = blockdata.task(2)+1;
        n = blockdata.n(2); matches = sum(blockdata.nback==true|blockdata.detect==true); %number of intended matches
        storage = NaN(1,n); %empty storage for n-back
        blockcost = cost(task); %pull remembered cost from lookup table
        for trial = 2:height(blockdata)
            stim = blockdata.stimnum(trial);
            if n == 0 %do detection
                if stim == 3
                    nmatches = nmatches + (blockdata.detect(trial) == true & blockdata.correct(trial)==true);
                end
                maintained(trial,:) = 0;
            else %n-backs and n-detects
                if stim == storage(1) % a match!
                    nmatches = nmatches + 1;
                else % no match!
                    if blockdata.nback(trial)==true
                        nmisses = nmisses + 1; %count misses
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
        maintained(1,:) = []; maintained = nanmean(sum(maintained>0,2));
        blockcost = blockcost + maintained*main; %average maintenance by cost of that
        blockcost = blockcost + nmatches*match;
        blockcost = blockcost + nupdates*inter; %interfering cost attributed to updating cost
        blockcost = blockcost - nmisses*miss;
        cost(task) = cost(task) + alpha*blockcost;
        real_cost(task) = blockdata.BDM(1); %update actual BDM value on each round
        not = [1:n_tasks];
        not(task) = [];
        cost(not) = cost(not) - forget*(cost(not)-init); %forget tasks that weren't just completed
        
        components = [components; subj nmatches maintained nupdates nmisses blockdata.display(1)+1 task cost(task) real_cost(task)];
    end %end of block by block loop
    cost = (cost/25)+1; %go from 1 to 5
    real_cost = (real_cost/25)+1; %go from 1 to 5

%     scatter(1:length(tasklabels),cost,[],subjcolors(subj,:),'Filled')
%     hold on
%     scatter(1:length(real_cost),real_cost,[],subjcolors(subj+1,:),'Filled')
%     legend({'Simulated Costs','Real Costs','Simulated Costs','Real Costs'})
%     xticks([1:4])
%     xticklabels(tasklabels)
%     title('Final Costs of Tasks - 1 Sim Subj')

end

toanalyze = table;
toanalyze.subj = components(:,1); toanalyze.nmatches = components(:,2);
toanalyze.maintained = components(:,3); toanalyze.nupdates = components(:,4);
toanalyze.nmisses = components(:,5); toanalyze.display = components(:,6); 
toanalyze.task = components(:,7); toanalyze.BDM = components(:,9); 
% pop things labelled into a table for clarity
% need to alter BDM column to be reflective of post-task costs (i.e.
% stagger by one round)
toanalyze.newBDM = NaN(height(toanalyze),1);

BDMrows = NaN(length(unique(toanalyze.subj)),32); simBDMrows = NaN(length(unique(toanalyze.subj)),32); %also, make meanable variable for BDMs/block
for subj = 1:length(unique(toanalyze.subj))
    newcolumn = NaN(32,1);
    if subj ~= 8 %trim this out for now because for some reason it's broken
    oldcolumn = toanalyze.BDM(toanalyze.subj==subj);
    tasklist = toanalyze.task(toanalyze.subj==subj);
    displayed = toanalyze.display(toanalyze.subj==subj);
    for task = 2:4
        realvalues = [];
        completed = find(tasklist==task);
        score = find(displayed==task);
        for i = 1:length(completed)
            next = score(score>completed(i)); 
            if length(next)>0
                next = next(1); %immediately following round
                realvalues = [realvalues; oldcolumn(next)];
            else
                realvalues = [realvalues; NaN];
            end
        end
        newcolumn(completed) = realvalues; %shift everything over by one, last value is now NaN since it influences nothing
    end
    toanalyze.newBDM(toanalyze.subj==subj) = newcolumn;
    BDMrows(subj,:) = newcolumn';
    simBDMrows(subj,:) = components(components(:,1)==subj,8)';
    end
end

%save('simdata/toanalyze.mat','toanalyze')

nsubjs = length(unique(toanalyze.subj))-1;
subjcolors = rand(nsubjs*2);subjcolors(:,4:end) = []; %delete unnecessary columns
taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
%plot one subj's learning block by block, to see how the dynamics of this
%go
realvssim = [toanalyze.subj toanalyze.display toanalyze.newBDM components(:,8)];
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




