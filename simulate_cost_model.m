function [simdata] = simulate_cost_model(modeltosim,allparams,toanalyze)
%simulate_cost_model Simulate specified cost model with specified
%parameters
%Get simdata out, which is plottable and analyzable
    %lower level parameters
    global real_epsilon_opt realparamlist model
    nparams = modeltosim.nparams; nsubjs = size(allparams,1);
    %global simdata %initialize empty
    simdata = [];
    tasks = unique(toanalyze.task);
    for subj = 1:nsubjs  %simulate the model for one subject at a time
    real_epsilon = []; params = allparams(subj,:);
    
    model = modeltosim;
    [uc,epsilon,init,mc,mainc,matchc,noisec,respc,lurec,alpha,delta] = setParamValues(params);
    costs = [uc mc mainc matchc noisec respc lurec];

    ratings = init.*(ones(1,max(toanalyze.display))); %init for each subject
    trialScalar = 1;
    %onesubj = repmat(toanalyze(toanalyze.subj==ceil(rand().*30),:),trialScalar,1); %if simulating some fake people
    onesubj = repmat(toanalyze(toanalyze.subj==subj,:),trialScalar,1); %if
    %using only real subjects
    nupdates = zeros(length(onesubj.nupdates),1); nupdates(onesubj.nupdates>0,:) = zscore(onesubj.nupdates(onesubj.nupdates>0,:)); % need to edit nupdates because it has so many zeros from irrelevant task 1
    nmisses = zscore(onesubj.nmisses); nmaintained = zscore(onesubj.maintained); nmatches = zscore(onesubj.nmatches);
    noisiness = zscore(onesubj.noisiness); responses = zscore(onesubj.nresponses); nlures = zscore(onesubj.nlures);
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
        if modeltosim.delta & trial > 1
            costs = setNewCosts(costs,delta,trial);
        end
        task = onesubj.task(trial);
        updates = nupdates(trial);misses = nmisses(trial); mains = nmaintained(trial); matches = nmatches(trial); noise = noisiness(trial); nresp = responses(trial); lures = nlures(trial);
        cost = costs*[nupdates(trial);nmisses(trial);nmaintained(trial);nmatches(trial);noisiness(trial);responses(trial);nlures(trial)]; %add all the costs together
        if modeltosim.alpha
            ratings(task) = ratings(task) + alpha*(cost-ratings(task)); %delta rule
        else %alpha fixed or no alpha
            %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
            ratings(task) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
        end
        simdata = [simdata; subj task rating torate updates misses mains matches noise nresp lures];
    end %of one subj run-through
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
    end %of simulating all subjects
end

