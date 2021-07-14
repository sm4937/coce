%% Hand-written expectation maximization script
% Written by Sarah Master, 06/2020
clear all %start from scratch
%generate simulated data on WM storage/updates/maintenance model
load('simdata\toanalyze.mat')
nsubjs = 30;
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));

global modeltofit onesim fit_epsilon_opt real_epsilon_opt mu sigma realsubjectsflag realparamlist
realsubjectsflag = false;

modelStruct = struct;
model_names = {'mc_epsilon_init_alpha_respc_lurec'};
notes = 'same code fit to subj';
for name = 1:length(model_names)
    model_labels{name} = strrep(model_names{name},'_','-');
end
tic
genrec = 1:length(model_names);
%which model or models to recover?
confusion_matrix = [];
for gen = 1:length(genrec)
    %pull model to simulate
    modeltosim = coc_createModels(model_names{genrec(gen)});
    nparams = modeltosim.nparams;

    %global simdata %initialize empty
    simdata = []; realparamlist = []; 
    %ground truth high-level parameters (distributions of parameters)
    mu = 0.1.*ones(nparams,1); sigma = 0.2.*ones(nparams,1);
    if modeltosim.init; idx = find(contains(modeltosim.paramnames,'init')); mu(idx) = 0.6; sigma(idx) = 0.1; end %set init to real value
    if modeltosim.alpha; idx = find(contains(modeltosim.paramnames,'alpha')); mu(idx) = 0.05; sigma(idx) = 0.03; end %set alpha specially
    params = [];
    for p = 1:nparams
        params(:,p) = normrnd(mu(p),sigma(p),nsubjs,1); %one column of parameters at a time
    end
    if modeltosim.alpha
        idx = find(contains(modeltosim.paramnames,'alpha'));
        params(params(:,idx) < 0,idx) = 0;
    end
    if modeltosim.epsilon
        idx = find(contains(modeltosim.paramnames,'epsilon'));
        params(params(:,idx) <= 0,idx) = 0.05;
    end
    realparamlist = params; 
    realmu = mean(realparamlist);
    realsigma = std(realparamlist);
    simdata = simulate_cost_model(modeltosim,realparamlist,toanalyze);
    
    % What's the correlatedness of each of my components?
    [r,p] = corr(simdata(:,5:11))
    
    %% recover subject parameters, and distribution of those parameters
    %pick the model to recover
    for rec = 1:length(genrec)
    modelscores = []; thrownout = []; %initialize for each model you're going to try to recover
    modeltofit = coc_createModels(model_names{genrec(rec)});
    nparams = modeltofit.nparams;
    
    %recover!
    starting = [];
    if rec>1
        nestedflag = checkNested(model_names{genrec(rec)},model_names{genrec(rec-1)});
        if nestedflag %start with guesses from previously fit, more general model
            eval(['starting = modelStruct.' model_names{genrec(rec-1)} '.lowparams;'])
        end
    end
    modeltofit = EMfit_sm(simdata,modeltofit,starting); %RUN FITTING HERE!
    lowparams = modeltofit.lowparams; bestmu = modeltofit.highparams(1:nparams); bestsigma = modeltofit.highparams(nparams+1:end);
    MAP = modeltofit.map;
    modelscores = [modelscores; MAP];
    % % MODEL COMPARISON SCORES % %
    %calculate relative model fit by getting mean LLH with parameters sampled
    %from best fitting distribution, with fit parameters
%     k = nsubjs; llhfromsample = []; %holding modeltofit constant
%     for sample = 1:k
%         params = [];
%         for i = 1:nparams
%             params(i) = normrnd(bestmu(i),bestsigma(i));
%             params(params<bounds(1)) = 0; params(params>1) = 1;
%         end
%         [~,llhfromsample(sample)] = getprobs_costlearning(params);
%     end
%     confusion_matrix(gen,rec) = mean(llhfromsample);
    
    if rec==gen
    figure
    subplot(4,2,1)
    bar([bestmu - realmu])
    xticklabels(modeltofit.paramnames)
    ylabel('Mus - distance from truth')
    title(['High-level fits'])
    
    for p = 1:nparams
        subplot(4,2,p+1)
        scatter(realparamlist(:,p),lowparams(:,p),'Filled');
        hold on
        name = modeltofit.paramnames{p};
        title([name 's: ' notes])
        xlabel('Real')
        ylabel('Fit')
        %legend(num2str(bestmu(p)))
        fig = gcf; fig.Color = 'w';
        plot([0 max(lowparams(:,p))], [0 max(lowparams(:,p))],'--')
    end %of nparam loop
    for subj = 1:nsubjs
        mu = bestmu; sigma = bestsigma;
        onesim = simdata(simdata(:,1)==subj,:);
        bestfit(subj) = getprobs_costlearning(lowparams(subj,:));
    end
    if ~modeltofit.epsilon
        subplot(4,2,p+2)
        scatter(real_epsilon_opt,fit_epsilon_opt,'Filled')
        hold on
        plot([0 xlim], [0 xlim],'--')
        title(['Fit epsilons versus real ones - final'])
        ylabel('Fit')
        xlabel('Real')
        legend(num2str(mu(p)))
    end
    end %of plotting if statement

    confusion_matrix(gen,rec) = min(modelscores);
    eval(['modelStruct.' model_names{genrec(rec)} ' = modeltofit;'])
    end % of recover loop
end %of generate loop
save(['simdata/genrectomodel_model' num2str(genrec(gen))],'')
time_elapsed = [num2str(toc/60) ' minutes elapsed.']
nparamlist = [];
for label = 1:length(model_labels) %label models on figures in a clearer way
    nparamlist(label) = length(strsplit(model_labels{label},'-'));
end

%AICs = confusion_matrix*2 + 2*nparamlist(1,1:size(confusion_matrix,2));

figure
heatmap(model_labels(genrec),model_labels(genrec),confusion_matrix);
title('Model confusion matrix')
xlabel('recovering model')
ylabel('generating model')

%% model simulations and validation figures
taskcolors = [0.75 0.75 0.75;0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75];
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];

figure
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
        curve = NaN(1,100);
        curve(1:sum(onesubj(:,2)==task+1)) = onesubj(onesubj(:,2)==task+1,3)';
        eval(['task' num2str(task) ' = [task' num2str(task) '; curve];'])
    end
end
errorbar(nanmean(task1),nanstd(task1)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
hold on
errorbar(nanmean(task2),nanstd(task2)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
errorbar(nanmean(task3),nanstd(task3)/sqrt(nsubjs),'Color',taskcolors(4,:),'LineWidth',2)
title('Group learning curves for each task')

%plot likelihood surface for each parameter
subj = 9;
onesim = simdata(simdata(:,1) == subj,:);
[r,~] = corr(onesim(:,5:9))
realparams = realparamlist(subj,:);
mu = realmu; sigma = realsigma;

samples = 0:0.01:1;
for p = 1:nparams
    subplot(5,2,p+2)
    constant = true(1,nparams); constant(p) = false;
    for sample = 1:length(samples)
        moving = samples(sample);%.*scalingvector(p);
        values = realparams; values(~constant) = moving;
        [~,llh] = map_costlearning(values);
        scatter(moving,llh)
        hold on
    end
    title(['Likelihood surface of simulation over ' modeltofit.paramnames{p}])
    plot(realparams(p),100,'*k')
    values = realparams;
    lillow = values; lillow(p) = values(p)-0.05; 
    lilhi = values; lilhi(p) = values(p)+0.05; 
    flatness = [getprobs_costlearning(lillow) getprobs_costlearning(values) getprobs_costlearning(lilhi)];
    disp(['flatness of ' modeltofit.paramnames{p} ' : ' num2str(flatness)])
end
[comparison(1),comparison(2)] = getprobs_costlearning(values);
disp(['neg LLH: ' num2str(comparison(1)) ', neg MAP: ' num2str(comparison(2))])

subplot(5,2,p+3)
bar(sqrt(nanmean((lowparams-realparamlist).^2)))
hold on
errorbar(sqrt(nanmean((lowparams-realparamlist).^2)),bestsigma)
%or errorbar with modeltofit.secondmoments?
title('MSE for each parameter')
xticklabels(modeltofit.paramnames)
fig = gcf; fig.Color = 'w';

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
        xlabel(modeltofit.paramnames(p))
        ylabel(modeltofit.paramnames(z))
        [r,pval] = corr(lowparams(:,p),lowparams(:,z));
    end
end
fig = gcf; fig.Color = 'w';

modeltosim = modeltofit;
newsimdata = simulate_cost_model(modeltosim,lowparams,toanalyze);
for subj = 1:nsubjs
    onesim = newsimdata(newsimdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = simdata(simdata(:,1)==subj,:);
    real = onereal(:,3);
    MSE(subj) = sqrt(nanmean((real-simulated).^2));
end
subjrankings = sortrows([MSE' (1:nsubjs)'],1);