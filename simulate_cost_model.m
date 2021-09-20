function [simdata] = simulate_cost_model(modeltosim,allparams,toanalyze)
%simulate_cost_model Simulate specified cost model with specified
%parameters
%Get simdata out, which is plottable and analyzable
    %lower level parameters
    global real_epsilon_opt  
    nsubjs = size(allparams,1);
    %global simdata %initialize empty
    simdata = [];
    subjnums = unique(toanalyze.subj);
    for subj = 1:nsubjs  %simulate the model for one subject at a time
    real_epsilon = []; params = allparams(subj,:);
    subjnum = subjnums(subj);
    %this is important anymore or a dinosaur of a previous implementation
    %of something
    [uc,epsilon,init,mc,mainc,matchc,noisec,respc,lurec,alpha,delta] = setParamValues(params,modeltosim);
    costs = [uc mc mainc matchc noisec respc lurec];

    ratings = [1 init].*(ones(1,max(toanalyze.display))); %init for each subject
    trialScalar = 1;
    onesubj = repmat(toanalyze(toanalyze.subj==subjnum,:),trialScalar,1); %if
    %using only real subjects
    if modeltosim.alpha
        nupdates = zeros(length(onesubj.nupdates),1); nupdates(onesubj.nupdates>0,:) = zscore(onesubj.nupdates(onesubj.nupdates>0,:)); % need to edit nupdates because it has so many zeros from irrelevant task 1
        nmisses = zscore(onesubj.nmisses); nmaintained = zscore(onesubj.maintained); nmatches = zscore(onesubj.nmatches);
        noisiness = zscore(onesubj.noisiness); responses = zscore(onesubj.nresponses); nlures = zscore(onesubj.nlures);
    elseif (modeltosim.delta || modeltosim.deltai)
        nupdates = onesubj.nupdates; % need to edit nupdates because it has so many zeros from irrelevant task 1
        nmisses = onesubj.nmisses; nmaintained = onesubj.maintained; nmatches = onesubj.nmatches;
        noisiness = onesubj.noisiness; responses = onesubj.nresponses; nlures = onesubj.nlures;        
    end
    ntrials = sum(~isnan(onesubj.BDM)&~isnan(onesubj.display)); %height(onesubj);
    for trial = 1:ntrials
        torate = onesubj.display(trial);
        if ~isnan(torate) %skip those trials
            rating = ratings(torate) + normrnd(0,epsilon); %mean 0, std epsilon
            rating(rating<0)=0; rating(rating>100) = 100;
            real_epsilon = [real_epsilon; subj rating ratings(torate) rating-ratings(torate)];
        else
            rating = NaN;
        end
        if (modeltosim.delta || modeltosim.deltai) & trial > 1
            costs = setNewCosts(costs,delta,trial);
            %figure(10);scatter(trial*ones(1,sum(costs~=0)),costs(costs~=0));hold on
        end
        task = onesubj.task(trial);
        updates = nupdates(trial);misses = nmisses(trial); mains = nmaintained(trial); matches = nmatches(trial); noise = noisiness(trial); nresp = responses(trial); lures = nlures(trial);
        cost = costs*[nupdates(trial);nmisses(trial);nmaintained(trial);nmatches(trial);noisiness(trial);responses(trial);nlures(trial)]; %add all the costs together
        if modeltosim.alpha
            ratings(task) = ratings(task) + alpha*(cost-ratings(task)); %delta rule
        else %alpha fixed or no alpha
            %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
            %ratings(task) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
            ratings(task) = ratings(task)+ 1*(cost-ratings(task));
        end
        simdata = [simdata; subjnum task rating torate updates misses mains matches noise nresp lures];
    end %of one subj run-through
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
    end %of simulating all subjects
end

