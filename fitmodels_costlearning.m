%% Hand-written expectation maximization function
% Written by Sarah Master, 06-08/2020
clear all %start from scratch
fitflag = false;
%generate simulated data on WM storage/updates/maintenance model
load('simdata\toanalyze.mat')
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));

global modeltofit onesim fit_epsilon_opt mu sigma realsubjectsflag bounds
realsubjectsflag = true;
bounds = [0 10];
%actions = stimuli;
%rows = actions
%columns = stimuli

model_names = {'uc','uc_epsilon','uc_epsilon_init','uc_init','uc_epsilon_init_mc', ...
    'uc_init_mc','uc_epsilon_init_mainc','uc_init_mainc','uc_epsilon_init_mc_mainc',... 
    'uc_epsilon_init_mc_mainc_alpha'};
for name = 1:length(model_names)
    model_labels{name} = strrep(model_names{name},'_','-');
end

modelstorec = [10]; 
%which model or models to recover?
modelStruct = {};

% recover subject parameters, and distribution of those parameters

%pick the model to recover
if fitflag == true
for rec = 1:length(modelstorec)
    thrownout = 500; 
    modeltofit = coc_createModels(model_names{modelstorec(rec)});
    output = EMfit_sm(toanalyze,modeltofit);
    eval(['modelStruct.' model_names{modelstorec(rec)} ' = output;'])
end %of generate loop
save(['simdata/fittosubjs_' model_names{modelstorec}])


nparamlist = [];
for label = 1:length(modelstorec)
    nparamlist(label) = length(strsplit(model_labels{label},'-'));
    modelname = model_names{modelstorec(label)};
    eval(['model_scores(label) = modelStruct.' modelname '.map;'])
end

%AICs = confusion_matrix*2 + 2*nparamlist(1,1:size(confusion_matrix,2));

figure
bar(model_scores)
title('Model scores')
xticklabels(model_labels(modelstorec))

[score,which] = min(model_scores);
best_name = model_names{modelstorec(which)};
best_model = eval(['modelStruct.' best_name]);
nparams = best_model.nparams;
save(['data/modelfits_' num2str(nparams) 'param_' date '.mat'],'best_model')
end %of fitting if

taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
    %display variable real quick ^ 
%% display, just like, parameters
load('data/modelfits_6param_14-Aug-2020.mat')
load('simdata/toanalyze.mat')
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));
nparams = best_model.nparams;
nsubjs = size(best_model.lowparams,1); lowparams = best_model.lowparams;
mu = best_model.highparams(1:nparams); sigma = best_model.highparams((nparams+1):end);
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0];

%how do the fits look?
figure
for p = 1:nparams
    subplot(4,2,p)
    scatter(1:nsubjs,lowparams(:,p),[],paramcolors(p,:),'Filled');
    name = best_model.paramnames{p};
    title(['Fit ' name 's'])
    xlabel('Subject')
    ylabel('Param value')
end %of nparam loop
subplot(4,2,nparams+1)
bar(1,mean(lowparams(:,1)),'FaceColor',paramcolors(1,:))
hold on
bar(2,mean(lowparams(:,2)),'FaceColor',paramcolors(2,:))
bar(3,mean(lowparams(:,3)),'FaceColor',paramcolors(3,:))
bar(4,mean(lowparams(:,4)),'FaceColor',paramcolors(4,:))
bar(5,mean(lowparams(:,5)),'FaceColor',paramcolors(5,:))
bar(6,mean(lowparams(:,6)),'FaceColor',paramcolors(6,:))
errorbar(nanmean(lowparams),nanstd(lowparams)/sqrt(nsubjs),'*k')
ylabel('mean value')
xticks([1:nparams])
xticklabels(best_model.paramnames)
fig = gcf; fig.Color = 'w';

figure
count = 0;
for p = 1:nparams
    for z = 1:nparams
        count = count+1;
        subplot(6,6,count)
        scatter(lowparams(:,p),lowparams(:,z),[],paramcolors(p,:),'Filled')
        lsline
        if z == p
            histogram(lowparams(:,p),'FaceColor',paramcolors(p,:))
        end
        xlabel(best_model.paramnames(p))
        ylabel(best_model.paramnames(z))
        [r,pval] = corr(lowparams(:,p),lowparams(:,z));
    end
end
fig = gcf; fig.Color = 'w';

%% model simulations and validation figures

simdata = simulate_cost_model(best_model,lowparams,toanalyze);
ntrials = sum(simdata(:,1)==1); nsubjs = length(unique(simdata(:,1)));
%scale things appropriately
simdata(:,3) = ((simdata(:,3).*bounds(2))./25)+1;

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
        curve(1:sum(onesubj(:,2)==task+1)) = onesubj(onesubj(:,2)==task+1,3)';
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
for i = 1:length(subjs)
    subj = subjs(i);
    subplot(3,2,i)
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.newBDM; real = (real./25)+1;
    for trial = 1:length(real)
        task = onereal.task(trial);
        scatter(trial,real(trial),[],taskcolors(task,:),'*')
        hold on
        scatter(trial,simulated(trial),[],taskcolors(task,:),'o','Filled')
    end
    title(['Subj ' num2str(subj)])
    legend({'Real','Simulated'})
end
fig = gcf; fig.Color = 'w';

%rank subjects by how they're fit by the model
for subj = 1:nsubjs
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.newBDM; real = (real./25)+1;
    MSE(subj) = sqrt(nanmean((real-simulated).^2));
end
subjrankings = sortrows([MSE' (1:nsubjs)'],1);

figure
for p = 1:nparams 
    subplot(3,2,p)
    for i = 1:length(subjrankings)
        subj = subjrankings(i,2);
        scatter(i,lowparams(subj,p),[],paramcolors(p,:),'Filled')
        hold on
        name = best_model.paramnames{p};
        title(['Fit ' name 's'])
        xlabel('Fit rank')
        ylabel('Param value')
    end
end
fig = gcf; fig.Color = 'w';

