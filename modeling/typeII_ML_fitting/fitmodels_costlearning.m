%% Hand-written expectation maximization function
% Written by Sarah Master, 06-08/2020
clear all %start from scratch
fitflag = false;
%generate simulated data on WM storage/updates/maintenance model
load('simdata\toanalyze.mat')

global modeltofit onesim fit_epsilon_opt mu sigma realsubjectsflag
realsubjectsflag = true;
%actions = stimuli;
%rows = actions
%columns = stimuli

%which model or models to recover?
model_names = {'respc_epsilon_init_alpha_lurec','respc_epsilon_init_alpha_lurec_mc'};
% mc lurec respc model BEST individual recovery so far
modelStruct = {};

% recover subject parameters, and distribution of those parameters
%pick the model to recover
if fitflag == true
for rec = 1:length(model_names)
    thrownout = 500; 
    starting = [];
    if rec>1
        nestedflag = checkNested(model_names{rec},model_names{rec-1});
        if nestedflag %start with guesses from previously fit, more general model
            eval(['starting = modelStruct.' model_names{rec-1} '.lowparams;'])
        end
    end
    modeltofit = coc_createModels(model_names{rec});
    output = EMfit_sm(toanalyze,modeltofit,starting); %RUN FITTING HERE!
    eval(['modelStruct.' model_names{rec} ' = output;'])
end %of generate loop

nparamlist = []; model_names = fieldnames(modelStruct);
for label = 1:length(model_names)
    modelname = model_names{label};
    eval(['nparamlist(label) = modelStruct.' modelname '.nparams;'])
    eval(['fitsbysubj = modelStruct.' modelname '.fitbysubj;'])
    eval(['modelStruct.' modelname '.AICsbysubj = (fitsbysubj*2) + (2*nparamlist(label));'])
    AICsbysubj(label,:) = (fitsbysubj*2)+(2*nparamlist(label));
end
AICs = mean(AICsbysubj,2); modelStruct.AICs = AICs;
[bestscore,which] = min(AICs);
best_model = eval(['modelStruct.' model_names{which}]);
modelStruct.best_model = best_model;
save(['data/modelfits_' num2str(length(model_names)) '_models_' date '.mat'],'modelStruct')
end %of fitting if

taskcolors = [0.75 0.75 0.75;0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75];
tasklabels = {'1-back','3-detect','2-back'};
    %display variable real quick ^ 
%% display parameters, model scores, etc.
%load('data/modelfits_2_models_23-Dec-2020.mat') %uc, matchc, one w delta
%load('data/modelfits_3_models_30-Dec-2020.mat') %uc, one w mainc, one w delta
%load('data/modelfits_3_models_02-Jan-2020.mat') %comparisons of uc, mainc,
%matchc models. doesn't appear to be any effect of mainc on costs
% uc is low cost, still
%load('data/modelfits_2_models_06-Jan-2021.mat') %uc, response c
%load('data/modelfits_2_models_13-Jan-2021.mat') %m c, response c, lure c
load('data/modelfits_3_models_09-Jan-2021_b.mat') %uc, responsec, lurec

% same as above to compare fits
load('simdata/toanalyze.mat')
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];

model_names = fieldnames(modelStruct); %post-loading, pull stuff out
model_names(contains(model_names,'AICs')) = []; model_names(contains(model_names,'best_model')) = [];
for name = 1:length(model_names)
    model_labels{name} = strrep(model_names{name},'_','-');
    eval(['model = modelStruct.' model_names{name} ';'])
    %extract AICs by subj
    AICsbysubj(name,:) = model.AICsbysubj(name,:);
end
AICs = modelStruct.AICs; bestscore = min(AICs);
best_model = modelStruct.best_model;
nparams = best_model.nparams; 
nsubjs = size(best_model.lowparams,1); lowparams = best_model.lowparams;
mu = best_model.highparams(1:nparams); sigma = best_model.highparams((nparams+1):end);
for p = 1:length(best_model.paramnames) %reformat param names for better plotting
    name = best_model.paramnames{p};
    paramnames{p} = strrep(name,'c',' costs');
    paramnames{p} = strrep(paramnames{p},'m','miss');
end


figure
subplot(1,2,1)
%plot group means
bar(AICs-bestscore)
title('Relative mean AIC scores by group')
xticklabels(model_labels)
xtickangle(45)

%plot AICs by subj
[~,tops] = min(AICsbysubj,[],1);
subplot(1,2,2)
measures = [];
for m = 1:length(model_names)
    measures(m) = sum(tops==m);
end
bar(measures)
xticklabels(model_labels)
xtickangle(45)
fig = gcf; fig.Color = 'w';
title('# of subjs best fit by each model')


%how do the fits look?
figure
for p = 1:nparams
    subplot(4,2,p)
    scatter(1:nsubjs,lowparams(:,p),[],paramcolors(p,:),'Filled');
    name = paramnames{p};
    title(['Fit ' name])
    xlabel('Subject')
    ylabel('Param value')
end %of nparam loop
subplot(4,2,nparams+1)
for p = 1:nparams
bar(p,mean(lowparams(:,p)),'FaceColor',paramcolors(p,:))
hold on
end
errorbar(nanmean(lowparams),nanstd(lowparams)/sqrt(nsubjs),'*k')
ylabel('mean value')
xticks([1:nparams])
xticklabels(paramnames)
fig = gcf; fig.Color = 'w';

figure; costs = [];
for p = 1:length(best_model.paramnames)
    costs = [costs; strcmp('c',best_model.paramnames{p}(end))];
end
costs = find(costs); 
trim = lowparams(:,1)>=(2.*std(lowparams(:,1))+mean(lowparams(:,1))); %remove outliers?
lowparams(trim,:) = [];
bar(mean(lowparams(:,costs)),'FaceColor',[0 0.7 0])
hold on
errorbar(nanmean(lowparams(:,costs)),nanstd(lowparams(:,costs))./sqrt(nsubjs),'*k')
plot(1:length(costs),lowparams(:,costs),'--k')
fig = gcf; fig.Color = 'w';
xticklabels(paramnames(costs))
xtickangle(30)
ylabel('Mean parameter value')
xlabel('Cost parameter')

%look at distribution of each parameter and their tradeoffs
figure
count = 0;
for p = 1:nparams
    for z = 1:nparams
        count = count+1;
        subplot(nparams,nparams,count)
        scatter(lowparams(:,p),lowparams(:,z),[],paramcolors(p,:),'Filled')
        lsline
        if z == p
            histogram(lowparams(:,p),'FaceColor',paramcolors(p,:))
        end
        xlabel(paramnames(p))
        ylabel(paramnames(z))
        [r,pval] = corr(lowparams(:,p),lowparams(:,z));
    end
end
fig = gcf; fig.Color = 'w';

%% model simulations and validation figures

simdata = simulate_cost_model(best_model,lowparams,toanalyze);
ntrials = sum(simdata(:,1)==1); nsubjs = length(unique(simdata(:,1)));
%scale things appropriately
simdata(:,3) = (simdata(:,3)./25)+1;

figure
subplot(2,1,1)
bar([nanmean(simdata(simdata(:,2)==2,3)) nanmean(simdata(simdata(:,2)==3,3)) nanmean(simdata(simdata(:,2)==4,3))])
hold on 
for subj = 1:nsubjs
    onesubj = simdata(simdata(:,1)==subj,:);
    errorbar([nanmean(onesubj(onesubj(:,2)==2,3)) nanmean(onesubj(onesubj(:,2)==3,3)) nanmean(onesubj(onesubj(:,2)==4,3))],[nanstd(onesubj(onesubj(:,2)==2,3)) nanstd(onesubj(onesubj(:,2)==3,3)) nanstd(onesubj(onesubj(:,2)==4,3))]/sqrt(length(onesubj)),'LineWidth',1);
end
title('Mean rating of tasks (sim)')
xticks([1:3])
ylim([1 5])
xticklabels({'1-back','2-back','3-back'})
xlabel('Task')
ylabel('Mean rating (sim)')

subplot(2,1,2)
task1 = []; task2 = []; task3 = [];
for task = 1:(ntasks-1)
    for subj = 1:nsubjs
        onesubj = simdata(simdata(:,1)==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onesubj(:,4)==task+1)) = onesubj(onesubj(:,4)==task+1,3)';
        eval(['task' num2str(task) ' = [task' num2str(task) '; curve];'])
    end
end
errorbar(nanmean(task1),nanstd(task1)/sqrt(nsubjs),'Color',taskcolors(1,:),'LineWidth',2)
hold on
errorbar(nanmean(task2),nanstd(task2)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
errorbar(nanmean(task3),nanstd(task3)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
legend({'1-back','2-back','3-detect'})
title('Group learning curves for each task')
fig = gcf; fig.Color = 'w';

figure
subjs = randperm(nsubjs,6);
subjs = [20 22 13 7 23 18];
tasksymbols = {'s','o','d'};
trialcolors = linspace(0.1,1,ntrials); trialcolors = repmat(trialcolors',1,3);
taskcount = false(1,3);
for i = 1:length(subjs)
    subj = subjs(i);
    subplot(3,2,i)
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3); 
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; real = (real./25)+1;
    for trial = 1:sum(~isnan(real))
        task = onereal.display(trial);
        if ~isnan(task)
            tasksymbol = tasksymbols{task-1};
            if ~taskcount(task-1)
                scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled','DisplayName',tasklabels{task-1})
                taskcount(task-1) = true; %have you already labeled this task?
            else
                scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled')
            end
            hold on
            line
            %re-formatting for Peter
%             scatter(trial,real(trial),[],taskcolors(task,:),tasksymbol,'Filled')
%             hold on
%             scatter(trial,simulated(trial),[],taskcolors(task,:),'o','Filled')
        end
    end
    title(['Subj ' num2str(subj)])
    %legend('location','best')
    %legend({'Real','Simulated'})
    xlabel('real values'); ylabel('simulated');
end
fig = gcf; fig.Color = 'w';

%rank subjects by how they're fit by the model
for subj = 1:nsubjs
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; real = (real./25)+1;
    idx = true(length(simulated),1); idx(isnan(real(idx))) = false; 
    MSE(subj) = sqrt(nanmean((real(idx)-simulated(idx)).^2));
end
subjrankings = sortrows([MSE' (1:nsubjs)'],1);

figure
for p = 1:nparams 
    subplot(4,2,p)
    for i = 1:length(subjrankings)
        subj = subjrankings(i,2);
        scatter(i,lowparams(subj,p),[],paramcolors(p,:),'Filled')
        hold on
        name = paramnames{p};
        title(['Fit ' name])
        xlabel('Fit rank')
        ylabel('Param value')
    end
end
fig = gcf; fig.Color = 'w';