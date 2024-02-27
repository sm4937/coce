function [simdata] = simulate_cost_model(modeltosim,allparams,toanalyze)
%simulate_cost_model Simulate fake data with pre-set parameters, pre-set
%model, pre-set parameter values. Primarily for model testing and for model
%validation post-fitting on real subject data.

%simdata, the output variable, is formatted the same way as the real
%subject data, for ease of comparison across dataframes.

% modeltosim is a structure much like modeltofit, which specifies which
% parameters are in play
% allparams is a list of all parameters for that model, for all subjects
% (so should be N x nparams in size)
% toanalyze is a MATLAB table containing the real subject
% responses/behavior/pseudo-randomized task conditions, etc.

global real_epsilon_opt
% this is updated at the very bottom of the script

%global simdata %initialize empty
simdata = [];

subjnums = unique(toanalyze.subj);
nsubjs = length(subjnums);

for subj = 1:nsubjs  %simulate the model for one subject at a time
    real_epsilon = [];
    params = allparams(subj,:);
    % select parameter values for that subject
    subjnum = subjnums(subj);
    
    [uc,epsilon,init,missc,mainc,matchc,noisec,respc,lurec,errorc,fac,rtc,otherlurec,alpha,delta] = set_param_values(params,modeltosim);

    % you'll notice that the same functions are applied here and in the
    % fitting function, getprobs_costlearning
    % this is an attempt to standardize param scaling, model specification,
    % etc. across different points in the model simulation/fitting pipeline
        
    costs = [uc missc mainc matchc noisec respc lurec errorc fac rtc otherlurec];
    % put all costs in one vector for ease of transformation in
    % set_new_costs (the cost-changing models)
    costs_0 = costs; 
    % grab first set of costs for figuring into the delta models
    
    ratings = [1 init].*(ones(1,max(toanalyze.display)));
    % initialize ratings for this subject based on the init free
    % parameter(s) (either one per rated task, or one for all three rated
    % tasks)
    
    trialScalar = 1;
    onesubj = repmat(toanalyze(toanalyze.subj==subjnum,:),trialScalar,1);
    % using this, you can test whether increasing the number of trials by
    % trialScalar amount increases fitting fidelity (it does, but not by a
    % lot, for most models). just increase trialScalar by the amount you
    % want the trial number to be multipled by.
    
    % z-score all components to get them being normally distributed
    % much better for comparison across cost magnitudes, model recovery,
    % etc.
    nupdates = zscore(onesubj.nupdates);
    nmisses = zscore(onesubj.nmisses); 
    nmaintained = zscore(onesubj.maintained);
    nmatches = zscore(onesubj.nmatches);
    noisiness = zscore(onesubj.noisiness); 
    responses = zscore(onesubj.nresponses);
    nlures = zscore(onesubj.nlures);
    nerrors = zscore(onesubj.nerrors); 
    nFAs = zscore(onesubj.nFAs);
    if sum(contains(onesubj.Properties.VariableNames,'RTs'))>0
        RTs = zscore(onesubj.RTs);
        notherlures = zscore(onesubj.notherlures);
        pct_time_elapsed = onesubj.pct_time_elapsed;
        pct_iters_complete = onesubj.pct_iters_complete;
    else
        RTs = zeros(height(onesubj),1);
        notherlures = zeros(height(onesubj),1);
        pct_time_elapsed = zeros(height(onesubj),1);
        pct_iters_complete = zeros(height(onesubj),1);
    end
    
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
            % cost-changing).
            % we then add some noise using a gaussian noise process
            % centered on 0, with a std of epsilon (a free parameter)
            
            rating(rating<0)=0; rating(rating>100) = 100;
            % bound by true rating values (can't be below 0 or above 100
            % in the experiment that subjects do online)
            
            real_epsilon = [real_epsilon; subj rating ratings(torate) rating-ratings(torate)];
            % this real_epsilon tracker is to calculate whether the value
            % of epsilon, due to stochasticity in data simulation, is close
            % in value to the value of epsilon I intended to simulate with.
            
        else
            
            rating = NaN;
            
        end
        
        % if it's a cost-changing (delta) model, update the costs according
        % to delta below
        time = pct_time_elapsed(trial);
        iter = pct_iters_complete(trial);
        
        % this evolving_costs variable can be used to plot the progression
        % of costs according to delta models
        evolving_costs(1,:) = costs_0;
        if sum(contains(modeltosim.paramnames,'delta'))>0 & trial > 1
            
            % which delta model are you running? the default is the linear
            % delta model from bioRxiv preprint 2023
            two_delta_flag = modeltosim.delta_time & modeltosim.delta_practice;
            quadratic_delta_flag = modeltosim.deltaexp;
            
            costs = set_new_costs(costs_0,delta,trial,quadratic_delta_flag,two_delta_flag,time,iter);
            evolving_costs(trial,:) = costs;
        end
        
        % which task did they just complete?
        % not always the same as they just rated
        task = onesubj.task(trial);
        
        %add all the costs together
        cost = costs*[nupdates(trial);nmisses(trial);nmaintained(trial);nmatches(trial);...
            noisiness(trial);responses(trial);nlures(trial);nerrors(trial);nFAs(trial);...
            RTs(trial);notherlures(trial)]; 
        % how costly was it, according to these parameter values?
        
        % learn from that cost (either with learning rate = alpha, or
        % learning rate = 1)
        if modeltosim.alpha
            
            ratings(task) = ratings(task) + alpha*(cost-ratings(task)); %delta rule
        
        else %alpha fixed or no alpha
            
            %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
            %ratings(task) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
            ratings(task) = ratings(task)+ 1*(cost-ratings(task));
        
        end
        
        simdata = [simdata; subjnum task rating torate nupdates(trial) nmisses(trial) nmaintained(trial) nmatches(trial)...
            noisiness(trial) responses(trial) nlures(trial) nerrors(trial) nFAs(trial) RTs(trial) notherlures(trial)...
            time iter];
        % concatenate the simulated data from this subject to the greater
        % simdata matrix
                
    end %of one subj run-through
    
    % when we were having trouble fitting epsilon, I began calculating the
    % true epsilon for each subject, to rule out the possibiliy that
    % randomness/task entropy wasn't contributing too much to our fitting
    % problems. This is no longer an issue, but in case it's useful I've left
    % that functionality in.
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
        
    % do you want to see how the total costs are evolving in response to
    % delta models? turn the plot flag on!
    plot_flag = false;
    
    if plot_flag & sum(contains(modeltosim.paramnames,'delta'))>0
        figure(11)
        scatter(1:ntrials,sum(evolving_costs,2),'LineWidth',1.5)
        set(gcf,'Color','w')
        set(gca,'FontSize',15)
    end
    
end % of simulating all subjects

end % of function

