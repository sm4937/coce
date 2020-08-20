%% Hand-written expectation maximization script
% Written by Sarah Master, 06/2020
clear all %start from scratch
%generate simulated data on WM storage/updates/maintenance model
load('simdata\toanalyze.mat')
nsubjs = 30;
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));

global modeltofit onesim fit_epsilon_opt noiselessri mu sigma bounds realsubjectsflag
realsubjectsflag = false;
bounds = [0 10]; %bounds on parameters
%actions = stimuli;
%rows = actions
%columns = stimuli

model_names = {'uc','uc_epsilon','uc_epsilon_init','uc_init','uc_epsilon_init_mc','uc_init_mc','uc_epsilon_init_mainc','uc_init_mainc','uc_epsilon_init_mc_mainc','uc_epsilon_init_mc_mainc_alpha'};
for name = 1:length(model_names)
    model_labels{name} = strrep(model_names{name},'_','-');
end
tic
genrec = [10]; %1:length(model_names);
%which model or models to recover?
confusion_matrix = [];
for gen = 1:length(genrec)
    %pull model to simulate
    %ground truth high-level parameters (distributions of parameters)
    mu = 0.5; sigma = 0.2; 
    %lower level fixed parameters
    modeltosim = coc_createModels(model_names{genrec(gen)});
    nparams = modeltosim.nparams;

    %global simdata %initialize empty
    simdata = []; realparamlist = []; 
    toexecute = 'normrnd(mu,sigma) ';
    for subj = 1:nsubjs  %simulate the model for one subject at a time
    eval(['params = [' repmat(toexecute,1,nparams) '];'])
    params(params<bounds(1)) = bounds(1); params(params>1) = 1; %center between 0 & 1/10/100
    realparamlist(subj,:) = params; real_epsilon = [];
    end
    realmu = mean(realparamlist);%repmat(mu,1,nparams); 
    realsigma = std(realparamlist);
    
    simdata = simulate_cost_model(modeltosim,realparamlist,toanalyze);

    %% recover subject parameters, and distribution of those parameters
    %pick the model to recover
    for rec = 1:length(genrec)
    modelscores = []; thrownout = []; %initialize for each model you're going to try to recover
    modeltofit = coc_createModels(model_names{genrec(rec)});
    nparams = modeltofit.nparams;
    
    %recover!
    modeltofit = EMfit_sm(simdata,modeltofit);
    lowparams = modeltofit.lowparams.*bounds(2); bestmu = modeltofit.highparams(1:nparams); bestsigma = modeltofit.highparams(nparams+1:end);
    MAP = modeltofit.map;
    modelscores = [modelscores; MAP];
    % % MODEL COMPARISON SCORES % %
    %calculate relative model fit by getting mean LLH with parameters sampled
    %from best fitting distribution, with fit parameters
    k = nsubjs; llhfromsample = []; %holding modeltofit constant
    for sample = 1:k
        params = [];
        for i = 1:nparams
            params(i) = normrnd(bestmu(i),bestsigma(i));
        end
        llhfromsample(sample) = llh_costlearning(params);
    end
    confusion_matrix(gen,rec) = mean(llhfromsample);
    
    if rec==gen

    figure
    subplot(4,2,1)
    bar([bestmu - realmu,bestsigma-realsigma])
    ylabel('Distance from truth')
    title(['High-level fits'])
    
    for p = 1:nparams
        subplot(4,2,p+1)
        scatter(lowparams(:,p),realparamlist(:,p),'Filled');
        name = modeltofit.paramnames{p};
        title(['Fit ' name 's versus real ones - final'])
        xlabel('Fit')
        ylabel('Real')
        legend(num2str(mu(p)))
        fig = gcf; fig.Color = 'w';
    end %of nparam loop
    for subj = 1:nsubjs
        mu = bestmu; sigma = bestsigma;
        onesim = simdata(simdata(:,1)==subj,:);
        bestfit(subj) = map_costlearning(realparamlist(subj,:));
    end
    if ~modeltofit.epsilon
        subplot(4,2,p+2)
        scatter(fit_epsilon_opt,real_epsilon_opt,'Filled')
        title(['Fit epsilons versus real ones - final'])
        xlabel('Fit')
        ylabel('Real')
        legend(num2str(mu(p)))
    end
    end %of plotting if statement

    if min(modelscores)>min(thrownout)
        oop = true
    end
    confusion_matrix(gen,rec) = min(modelscores);
    end % of recover loop
    
end %of generate loop
save(['simdata/genrectomodel' num2str(gen)])
time_elapsed = [num2str(toc/60) ' minutes elapsed.']
nparamlist = [];
for label = 1:length(model_labels)
    nparamlist(label) = length(strsplit(model_labels{label},'-'));
end

%AICs = confusion_matrix*2 + 2*nparamlist(1,1:size(confusion_matrix,2));

figure
heatmap(model_labels(genrec),model_labels(genrec),confusion_matrix)
title('Model confusion matrix')
xlabel('recovering model')
ylabel('generating model')
taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
    %display variable real quick ^ 

%% model simulations and validation figures

if length(genrec) == 1
figure
if ntasks < 2
    subplot(4,2,1)
    finalratings = simdata(ntrials:ntrials:end,3);
    scatter(realparamlist(:,1),finalratings)
    title('Final cost for task 1 by UC value')
    xlabel('UC')
    ylabel('Final cost of task 1')
    fig = gcf; fig.Color = 'w';

    subplot(4,2,2)
    %learning curve for a few subjects
    subjects = 1:nsubjs;
    subjects = subjects(randperm(nsubjs));
    for sim = 1:5
        toplot = simdata(simdata(:,1) == subjects(sim),3);
        scatter(1:ntrials,toplot,'Filled')
        alphalabels(sim) = cellstr(num2str(realparamlist(subjects(sim),1)));
        hold on
    end
    legend(alphalabels)
    xlabel('Block')
    ylabel('Rating')
    title('Cost learning curves by UC value')
else %if ntasks > 2
    subplot(4,2,1)
    bar([nanmean(simdata(simdata(:,2)==2,3)) nanmean(simdata(simdata(:,2)==3,3)) nanmean(simdata(simdata(:,2)==4,3))])
    hold on 
    for subj = 1:nsubjs
        onesubj = simdata(simdata(:,1)==subj,:);
        errorbar([nanmean(onesubj(onesubj(:,2)==2,3)) nanmean(onesubj(onesubj(:,2)==3,3)) nanmean(onesubj(onesubj(:,2)==4,3))],[nanstd(onesubj(onesubj(:,2)==1,3)) nanstd(onesubj(onesubj(:,2)==2,3)) nanstd(onesubj(onesubj(:,2)==3,3))]/sqrt(length(onesubj)),'LineWidth',1);
    end
    title('mean rating of tasks')
    xlabel('task')
    ylabel('mean rating')

    subplot(4,2,2)
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
    title('Group learning curves for each task')
end
%plot likelihood surface for each parameter
onesim = simdata(simdata(:,1) == 1,:);
realparams = realparamlist(1,:);
mu = realmu; sigma = realsigma;

samples = bounds(1):1:bounds(2);
for p = 1:nparams
    subplot(4,2,p+2)
    constant = true(1,nparams); constant(p) = false;
    for sample = 1:length(samples)
        moving = samples(sample);
        values = realparams; values(~constant) = moving;
        llh = llh_costlearning(values);
        scatter(moving,llh)
        hold on
    end
    title(['Likelihood surface of simulation over ' modeltofit.paramnames{p}])
    plot(realparams(p),1000,'*k')
    plot(0,0,0,ylim,'*k')
    plot(1,0,1,ylim,'*k')
end
end % of eliminating plots when full generate/recover being run 
