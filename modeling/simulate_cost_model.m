function [simdata] = simulate_cost_model(modeltosim,allparams,toanalyze)
%simulate_cost_model Simulate fake data with pre-set parameters, pre-set
%model, pre-set parameter values. Primarily for model testing and for model
%validation post-fitting on real subject data.

%simdata, the output variable, is formatted the same way as the real
%subject data, for ease of comparison across dataframes.

global real_epsilon_opt
% this is updated at the very bottom of the script

nsubjs = size(allparams,1);
%global simdata %initialize empty
simdata = [];
subjnums = unique(toanalyze.subj);
for subj = 1:nsubjs  %simulate the model for one subject at a time
    real_epsilon = []; params = allparams(subj,:);
    subjnum = subjnums(subj);
    %this is important anymore or a dinosaur of a previous implementation
    %of something
    [uc,epsilon,init,missc,mainc,matchc,noisec,respc,lurec,errorc,fac,alpha,delta] = set_param_values(params,modeltosim);
    costs = [uc missc mainc matchc noisec respc lurec errorc fac];
    % put all costs in one vector for ease of transformation in
    % set_new_costs (the cost-changing models)
    
    ratings = [1 init].*(ones(1,max(toanalyze.display))); %init for each subject
    trialScalar = 1;
    onesubj = repmat(toanalyze(toanalyze.subj==subjnum,:),trialScalar,1);
    
    
    % z-score all components to aid in them being normally distributed
    % much better for comparison across cost magnitudes, model recovery,
    % etc.
    nupdates = zeros(length(onesubj.nupdates),1); nupdates(onesubj.nupdates>0,:) = zscore(onesubj.nupdates(onesubj.nupdates>0,:)); % need to edit nupdates because it has so many zeros from irrelevant task 1
    nmisses = zscore(onesubj.nmisses); nmaintained = zscore(onesubj.maintained); nmatches = zscore(onesubj.nmatches);
    noisiness = zscore(onesubj.noisiness); responses = zscore(onesubj.nresponses); nlures = zscore(onesubj.nlures);
    nerrors = zscore(onesubj.nerrors); nFAs = onesubj.nFAs;
    
    ntrials = sum(~isnan(onesubj.BDM)&~isnan(onesubj.display)); %height(onesubj);
    for trial = 1:ntrials
        
        % cycle over each round of the experiment (round of each task)
        
        torate = onesubj.display(trial);
        % grab which task is currently being rated
        
        % some subject data was saved improperly (primarily towards the end
        % of the experiment, for some reason), resulting in NaNs where
        % there should have been task labels
        % here I exclude these trials from modeling analyses because they
        % throw errors
        if ~isnan(torate) 
            
            rating = ratings(torate) + normrnd(0,epsilon); %mean 0, std epsilon
            % generate a simulated rating:
            % the mean is the true learned rating of that task, using
            % whatever update mechanism has been specified (learning versus
            % cost-chaning).
            % we then add some noise using a gaussian noise process
            % centered on 0, with a std of epsilon (a free parameter)
            
            rating(rating<0)=0; rating(rating>100) = 100;
            real_epsilon = [real_epsilon; subj rating ratings(torate) rating-ratings(torate)];
        else
            rating = NaN;
        end
        if (modeltosim.delta || modeltosim.deltai) & trial > 1
            costs = set_new_costs(costs,delta,trial);
            %figure(10);scatter(trial*ones(1,sum(costs~=0)),costs(costs~=0));hold on
        end
        task = onesubj.task(trial);
        
        cost = costs*[nupdates(trial);nmisses(trial);nmaintained(trial);nmatches(trial);noisiness(trial);responses(trial);nlures(trial);nerrors(trial);nFAs(trial)]; %add all the costs together
        if modeltosim.alpha
            ratings(task) = ratings(task) + alpha*(cost-ratings(task)); %delta rule
        else %alpha fixed or no alpha
            %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
            %ratings(task) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
            ratings(task) = ratings(task)+ 1*(cost-ratings(task));
        end
        simdata = [simdata; subjnum task rating torate nupdates(trial) nmisses(trial) nmaintained(trial) nmatches(trial) noisiness(trial) responses(trial) nlures(trial) nerrors(trial) nFAs(trial)];
    end %of one subj run-through
    
    % when we were having trouble fitting epsilon, I began calculating the
    % true epsilon for each subject, to rule out the possibiliy that
    % randomness/task entropy wasn't contributing too much to our fitting
    % problems. This is no longer an issue, but in case it's useful I've left
    % that functionality in.
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
    
end % of simulating all subjects

end % of function

