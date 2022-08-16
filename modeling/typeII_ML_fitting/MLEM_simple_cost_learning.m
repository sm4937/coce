%% Hand-written expectation maximization function
% Written by Sarah Master, 06/2020
clear all %start from scratch
%generate data to recover later
%start with one parameter and 32 trials, simulated
ntrials = 33;
nsubjs = 30;
stimuli = [ones(ntrials/3,1) 2*ones(ntrials/3,1) 3*ones(ntrials/3,1)];
stimuli = stimuli(randperm(ntrials));
ntasks = length(unique(stimuli));
costs = 1+stimuli; %change to learning about just one task
%actions = stimuli;
%rows = actions
%columns = stimuli
%ground truth high-level parameters (distributions of parameters)
mu = 0.5; sigma = 0.2; 
%lower level fixed parameters
init = 1; epsilon = 0.7;

%pull model to simulate
modeltosim = coc_createModels('alpha_epsilon_init');
nparams = modeltosim.nparams;

%global simdata %initialize empty
simdata = [];
%simulate the model for one subject at a time
toexecute = 'normrnd(mu,sigma) ';
for subj = 1:nsubjs
eval(['params = [' repmat(toexecute,1,nparams) '];'])
params(params<0) = 0; params(params>1) = 1; %center between 0 and 1
realparamlist(subj,:) = params;
alpha = params(1); 
if modeltosim.epsilon
    epsilon = params(2);
end
if modeltosim.init
    init = params(3)*5;
end
ratings = init*(ones(1,length(unique(stimuli)))); %init for each subject
    for trial = 1:ntrials
        stim = stimuli(trial);
        cost = costs(trial);
        ratings(stim) = ratings(stim) + 0.3*(alpha*cost); %not a learning rate of one
        rating = ratings(stim) + normrnd(0,epsilon); %mean 0, std epsilon
        simdata = [simdata; subj stim rating alpha];
    end
end

realmu = repmat(mu,1,nparams); realsigma = repmat(sigma,1,nparams);
taskcolors = [0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75; 0.75 0.75 0.75];
%display variable real quick ^ 

figure
if ntasks < 2
    subplot(2,2,1)
    finalratings = simdata(ntrials:ntrials:end,3);
    scatter(realparamlist(:,1),finalratings)
    title('Final cost for task 1 by alpha value')
    xlabel('Alpha')
    ylabel('Final cost of task 1')
    fig = gcf; fig.Color = 'w';

    subplot(2,2,2)
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
    title('Cost learning curves by alpha value')
else %if ntasks > 2
    subplot(2,2,1)
    bar([nanmean(simdata(simdata(:,2)==1,3)) nanmean(simdata(simdata(:,2)==2,3)) nanmean(simdata(simdata(:,2)==3,3))])
    hold on 
    for subj = 1:nsubjs
        onesubj = simdata(simdata(:,1)==subj,:);
        errorbar([nanmean(onesubj(onesubj(:,2)==1,3)) nanmean(onesubj(onesubj(:,2)==2,3)) nanmean(onesubj(onesubj(:,2)==3,3))],[nanstd(onesubj(onesubj(:,2)==1,3)) nanstd(onesubj(onesubj(:,2)==2,3)) nanstd(onesubj(onesubj(:,2)==3,3))]/sqrt(length(onesubj)),'LineWidth',1);
    end
    title('mean rating of tasks')
    xlabel('task')
    ylabel('mean rating')
    
    subplot(2,2,2)
    task1 = []; task2 = []; task3 = [];
    for task = 1:ntasks
        for subj = 1:nsubjs
            onesubj = simdata(simdata(:,1)==subj,:);
            curve = onesubj(onesubj(:,2)==task,3)';
            eval(['task' num2str(task) ' = [task' num2str(task) '; curve];'])
        end
    end
    errorbar(nanmean(task1),nanstd(task1)/sqrt(nsubjs),'Color',taskcolors(1,:),'LineWidth',2)
    hold on
    errorbar(nanmean(task2),nanstd(task2)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
    errorbar(nanmean(task3),nanstd(task3)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
    title('Group learning curves for each task')
end

global onesim modeltofit
subplot(2,2,3)
onesim = simdata(simdata(:,1) == 1,:);
modeltofit = modeltosim;
alphas = [mu-(sigma*5):0.01:mu+(sigma*5)];
epsilon = 0.7; %hold epsilon constant
for sample = 1:length(alphas)
    alpha = alphas(sample);
    llh = llhRL([alpha epsilon init]);
    scatter(alpha,llh)
    hold on
end
plot(realmu,1000,'*k')
title('Likelihood surface of simulation over alphas')

subplot(2,2,4)
epsilons = [mu-(sigma*5):0.01:mu+(sigma*5)];
alpha = 0.3; %now hold alpha constant
for sample = 1:length(epsilons)
    epsilon = epsilons(sample);
    llh = llhRL([alpha epsilon init]);
    scatter(epsilon,llh)
    hold on
end
plot(realmu,1000,'*k')
title('Likelihood surface of simulation over epsilon')

%% recover subject parameters, and distribution of those parameters
% set fmincon options, settings for model fitting
options = optimoptions(@fmincon,'Display','off');
toexecute = 'normrnd(mu,sigma) ';

%pick the model to recover
global modeltofit
modeltofit = coc_createModels('alpha_epsilon_init');
nparams = modeltofit.nparams;

%recover!
pmin=zeros(1,nparams);pmax=ones(1,nparams);
niters = 15;
%starting with just alpha
bestmu = rand(1,nparams); bestsigma = rand(1,nparams);
for iters = 1:niters %whole EM algorithm
    mus = rand(10,nparams); sigmas = rand(10,nparams);
    mus(1,:) = bestmu; sigmas(1,:) = bestsigma;
    bestdist = [];
    for tries = 1:10 %of different mus/sigmas
        bestforsubj = [];
        mu = mus(tries,:); sigma = sigmas(tries,:);
        for subj = 1:nsubjs
            onesim = simdata(simdata(:,1) == subj,:);
            fitparams = [];
            for z = 1 %z random starts for each subject
                for i = 1:nparams
                    params(i) = normrnd(mu(i),sigma(i));
                end
                params(params<0) = 0; params(params>1) = 1; %center between 0 and 1
                %init = pmin + rand(1,length(pmin)).*(pmax-pmin);
                [bestparam,localmin]=fmincon(@llhRL,params,[],[],[],[],pmin,pmax);
                %negllh = llhRL(alpha);
                %fitparams = [fitparams; subj z mu sigma alpha negllh];
                %dbstop in MLEM_sm.m at 82 if ~isreal(negllh)
                fitparams = [fitparams; subj z bestparam localmin]; %for the iters over subject, which iters produced what
            end
            [bestllh,which] = min(fitparams(:,end));
            bestforsubj = [bestforsubj; subj fitparams(which,3:3+(nparams-1)) bestllh];
            %track subj num, alpha, llh
        end
        if rand()<0.1 %show 10% of plots
        figure
        subplot(2,2,1)
        scatter(bestforsubj(:,2),realparamlist(:,1),'Filled');
        title('Fit alphas versus real ones')
        xlabel('Fit alphas')
        ylabel('Real alphas')
        legend(num2str(mu(1)))
        fig = gcf; fig.Color = 'w';
        
        subplot(2,2,2)
        scatter(bestforsubj(:,3),realparamlist(:,2),'Filled');
        title('Fit epsilons versus real ones')
        xlabel('Fit epsilons')
        ylabel('Real epsilons')
        legend(num2str(mu(2)))
        fig = gcf; fig.Color = 'w';
        
        subplot(2,2,3)
        scatter(bestforsubj(:,4),realparamlist(:,3),'Filled');
        title('Fit inits versus real ones')
        xlabel('Fit inits')
        ylabel('Real inits')
        legend(num2str(mu(2)))
        fig = gcf; fig.Color = 'w';

        subplot(2,2,4)
        scatter(bestforsubj(:,4),bestforsubj(:,2),'Filled');
        title('Fit inits versus fit alphas')
        xlabel('Fit inits')
        ylabel('Fit alphas')
        [r,p] = corr(bestforsubj(:,2),bestforsubj(:,3));
        legend(['corr = ' num2str(r)])
        fig = gcf; fig.Color = 'w';
        end
        
        %re-infer mu and sigma?
        currentparams = bestforsubj(:,nparams:end-1);
%       newmu = mean(currentparams);
%       newsigma = std(currentparams);
        meanLLH = mean(bestforsubj(:,end));
        bestdist = [bestdist; mus(tries,:) sigmas(tries,:) meanLLH];
    end %of cycling through random/best mu and sigma starting points
    [llh,which] = min(bestdist(:,3));
    bestmu = bestdist(which,1:nparams); bestsigma = bestdist(which,nparams+1:2*nparams); modelscore = llh;
end %of trying to find the best mu and sigma

figure
bar([bestmu - realmu,bestsigma-realsigma])
title([num2str(niters) ' iterations'])
