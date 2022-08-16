%% Hand-written expectation maximization script
% Written by Sarah Master, 06/2020
clear all %start from scratch
%generate simulated data on WM storage/updates/maintenance model
load('simdata\toanalyze.mat')
nsubjs = 30;
tasks = toanalyze.task;
ntasks = length(unique(toanalyze.task));

global modeltofit onesim fit_epsilon_opt noiselessri mu sigma realsubjectsflag
realsubjectsflag = false;
%actions = stimuli;
%rows = actions
%columns = stimuli

model_names = {'uc','uc_epsilon','uc_epsilon_init','uc_init','uc_epsilon_init_mc','uc_init_mc','uc_epsilon_init_mainc','uc_init_mainc','uc_epsilon_init_mc_mainc'};
for name = 1:length(model_names)
    model_labels{name} = strrep(model_names{name},'_','-');
end
tic
genrec = [3]; %1:length(model_names);
%which model or models to recover?
confusion_matrix = [];
for gen = 1:length(genrec)
    %pull model to simulate
    %ground truth high-level parameters (distributions of parameters)
    mu = 0.5; sigma = 0.2; 
    %lower level fixed parameters
    init = 1; epsilon = 0.7; mc = 0; mainc = 0; alpha = 0.35;
    modeltosim = coc_createModels(model_names{genrec(gen)});
    nparams = modeltosim.nparams;

    %global simdata %initialize empty
    simdata = []; realparamlist = []; 
    toexecute = 'normrnd(mu,sigma) ';
    for subj = 1:nsubjs  %simulate the model for one subject at a time
    eval(['params = [' repmat(toexecute,1,nparams) '];'])
    params(params<0) = 0; params(params>1) = 1; %center between 0 and 1
    realparamlist(subj,:) = params; real_epsilon = [];
    uc = params(1); 
    if modeltosim.epsilon
        idx = find(contains(modeltosim.paramnames,'epsilon'));
        epsilon = params(idx);
        %epsilon = 0.2; %note to self: make sure hist of params is normally distrib
    end
    if modeltosim.init
        idx = find(contains(modeltosim.paramnames,'init'));
        init = params(idx)*5;
    end
    if modeltosim.mc
        idx = find(contains(modeltosim.paramnames,'mc'));
        mc = (params(idx)-0.5)*2;
    end
    if modeltosim.mainc
        idx = find(contains(modeltosim.paramnames,'mainc'));
        mainc = params(idx);
    end
    ratings = init*(ones(1,length(unique(tasks)))); %init for each subject
    onesubj = toanalyze(toanalyze.subj==subj,:);
    ntrials = height(onesubj);
    for trial = 1:ntrials
        task = onesubj.task(trial);
        updates = onesubj.nupdates(trial);
        misses = onesubj.nmisses(trial);
        misses = ceil(rand()*4);
        mains = onesubj.maintained(trial);
        mains = ceil(rand()*2);
        cost = uc*updates;
        if modeltosim.mc
            cost = cost + mc*misses;
        end
        if modeltosim.mainc
            cost = cost+ mainc*mains;
        end
        ratings(task) = ratings(task) + alpha*(cost); %not a learning rate of one
        noiselessri(trial) = ratings(task); 
        rating = ratings(task) + normrnd(0,epsilon); %mean 0, std epsilon
        real_epsilon = [real_epsilon; subj rating ratings(task) rating-ratings(task)];
        simdata = [simdata; subj task rating updates uc misses mc mains mainc];
    end %of one subj run-through
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
    end %of simulating all subjects

    realmu = repmat(mu,1,nparams); realsigma = repmat(sigma,1,nparams);

    %% recover subject parameters, and distribution of those parameters
    % set fmincon options, settings for model fitting
    options = optimoptions('fmincon','Display','off');
    toexecute = 'normrnd(mu,sigma) ';
    %pick the model to recover
    for rec = 1:length(genrec)
    modelscores = []; thrownout = []; %initialize for each model you're going to try to recover
    modeltofit = coc_createModels(model_names{genrec(rec)});
    nparams = modeltofit.nparams;
    
    %modeltofit = EMfit_sm(simdata,modeltofit);
    %lowparams = modeltofit.lowparams; bestmu = modeltofit.highparams(1:nparams); bestsigma = modeltofit.highparams(nparams+1:end);
    
    %put plots here
    %recover!
    pmin=zeros(1,nparams);pmax=ones(1,nparams);
    niters = 10; 
    %starting with just alpha
    bestmu = rand(1,nparams); bestsigma = rand(1,nparams);  
    highlconvergence = [];
    count = 0; thrownout = 500; %count number of outliers being fit 
    for iters = 1:niters %whole EM algorithm
        mus = rand(10,nparams); sigmas = rand(10,nparams); sigmas(sigmas<0.05)=0.05;
        mus(1,:) = bestmu; sigmas(1,:) = bestsigma; 
        bestdist = []; paramsbytry = []; smbytry = []; lowlconvergence = []; bestmus = []; bestsigmas = [];
        for tries = 1:10 %of different mus/sigmas
            bestforsubj = [];
            mu = mus(tries,:); sigma = sigmas(tries,:);
            for subj = 1:nsubjs
                onesim = simdata(simdata(:,1) == subj,:);
                fitparams = [];
                for z = 1 %z random starts for each subject
                    params = [];
                    for i = 1:nparams
                        params(i) = normrnd(mu(i),sigma(i));
                    end
                    params(params<0) = 0; params(params>1) = 1;
                    %params = realparamlist(subj,:);
                    
                    [bestparam,map,~,~,~,~,hess]=fmincon(@map_costlearning,params,[],[],[],[],pmin,pmax,[],options);
                    secondmoment = diag(inv(hess))';
                    %dbstop in MLEM_sm.m at 210 if sum(secondmoment < 0)>0
                    fitparams = [fitparams; subj z bestparam secondmoment map]; %for the iters over subject, which iters produced what
                end
                [bestmap,which] = min(fitparams(:,end));
                bestforsubj = [bestforsubj; subj fitparams(which,3:3+(nparams-1)) fitparams(which,3+nparams:end-1) bestmap];
                %track subj num, alpha, llh
            end
            if bestforsubj(:,3) == 0 & modeltofit.epsilon
                count = count + 1; %count, per model, how many times epsilon goes to 0
            end
            %re-infer mu and sigma?
            currentparams = bestforsubj(:,2:2+(nparams-1));
            currentsecondmoments = bestforsubj(:,(2+nparams):end-1);
            paramsbytry{tries} = currentparams;
            smbytry{tries} = currentsecondmoments;
            %newmu = mean(currentparams);
            %newsigma = std(currentparams);
            %meanLLH = mean(bestforsubj(:,end));
            mapsubjs = mean(bestforsubj(:,end)); %multiply all subj maps or add all llh's
            bestdist = [bestdist; mus(tries,:) sigmas(tries,:) mapsubjs];
            if length(genrec)==1; lowlconvergence(tries,:) = mean(currentparams-realparamlist);end
        end %of cycling through random/best mu and sigma starting points
        if length(genrec) ==1
            figure
%             if modeltofit.epsilon
%                 toplot = [];
%                 toplot = lowlconvergence;
%                 toplot(:,2) = 1./(1+exp(-toplot(:,2)));
%             end
            for p = 1:nparams
                plot(1:tries,lowlconvergence(:,p))
                hold on
            end
            legend(modeltofit.paramnames)
            ylabel('Mean distance from truth')
        end
        [MAP,which] = min(bestdist(:,end));
        bestmu = mean(paramsbytry{which}); 
        bestsigma = sqrt(mean(paramsbytry{which}.^2 + smbytry{which})-(bestmu.^2)); %variance
        %dbstop in MLEM_sm.m at 161 if sum(~isreal(bestsigma))>0
        while sum(~isreal(bestsigma))>0 %ill-fitting?
            thrownout = [thrownout; MAP];
            bestdist(which,:) = []; %erase last minimum, poorly fitting 
            idxes = 1:tries; idxes(which) = [];
            [MAP,which] = min(bestdist(:,end)); %replace with second-best
            bestmu = mean(paramsbytry{idxes(which)});
            bestsigma = sqrt(mean(paramsbytry{idxes(which)}.^2 + smbytry{idxes(which)})-(bestmu.^2)); %variance
        end
        if iters==niters & (sum(contains(modeltofit.paramnames,modeltosim.paramnames))==nparams) %show final subject-level fits
        figure
        for p = 1:nparams
            subplot(3,2,p)
            lowparams = bestforsubj(:,2:2+(nparams-1));
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
            subplot(2,2,p+1)
            scatter(fit_epsilon_opt,real_epsilon_opt,'Filled')
            title(['Fit epsilons versus real ones - final'])
            xlabel('Fit')
            ylabel('Real')
            legend(num2str(mu(p)))
        end
        end %of plotting if statement
        modelscores = [modelscores; MAP];
        highlconvergence(iters,:) = [bestmu bestsigma];
    end %of trying to find the best mu and sigma
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

    if length(genrec) == 1
        figure
        subplot(1,2,1)
        bar([bestmu - realmu,bestsigma-realsigma])
        ylabel('Distance from truth')
        title(['Final answer ' num2str(niters) ' iters'])
        
        subplot(1,2,2)
        plot(1:niters,highlconvergence(:,1:nparams)-realmu)
        hold on
        plot(1:niters,highlconvergence(:,nparams+1:end)-realsigma)
        title(['Convergence of high level parameters'])
        legg = repmat({'Mu'},1,nparams); legg(nparams+1:nparams*2)= repmat({'Sigma'},1,nparams);
        legend(legg)
        xlabel('Iteration of EM')
        ylabel('Distance from truth')
        fig = gcf; fig.Color = 'w';
    end

    [modeltofit.map,which] = min(modelscores);
    modeltofit.lowparams = bestforsubj(:,2:2+(nparams-1));
    modeltofit.highparams = [bestmu bestsigma];
    %save all of this somewhere useful
    genrec_highlevelparams{gen,rec} = [bestmu bestsigma];
    if min(modelscores)>min(thrownout)
        oop = true
    end
    %confusion_matrix(gen,rec) = min(modelscores);
    end % of recover loop
    
end %of generate loop
save(['simdata/genrectomodel' num2str(gen) '_' num2str(niters) 'iters'])
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

samples = 0:0.01:1;
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
    plot(pmin(p),0,pmin(p),ylim,'*k')
    plot(pmax(p),0,pmax(p),ylim,'*k')
end
end % of eliminating plots when full generate/recover being run 
