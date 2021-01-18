function [negllh] = llh_costlearning(params)
%Get negative llh for specific parameter values of simple RL model
%   Runs model with input value of parameters
%   Returns llh

global onesim modeltofit fit_epsilon_opt realsubjectsflag model
model = modeltofit;

if ~realsubjectsflag
% simulated data, not real
% simdata = [simdata; subj task rating torate updates misses mains];
    subj = unique(onesim(:,1));
    stimuli = onesim(:,2);
    realratings = onesim(:,3);
    display = onesim(:,4);
    nupdates = onesim(:,5);
    nmisses = onesim(:,6);
    mains = onesim(:,7);
    nmatches = onesim(:,8);
% real data format
% toanalyze = 1. subj num, 2. nmatches, 3. maintained, 4. nupdates
%   5. nmisses, 6. task displayed, 7. task executed, 8. BDM, 9. accurate
%   BDM for that task
else %fitting real subject data
    subj = unique(onesim.subj);
    stimuli = onesim.task; %?? task?
    realratings = (onesim.BDM); %scale down? just added this 11/08/2020 - see how this works
    display = onesim.display;
    nupdates = zscore(onesim.nupdates); %./max(onesim.nupdates);
    nmisses = zscore(onesim.nmisses); %./max([max(onesim.nmisses) 1]); %if nmisses = 0, divide by 1 instead 
    mains = zscore(onesim.maintained); %./max(onesim.maintained); 
    nmatches = zscore(onesim.nmatches);
end

[uc,epsilon,init,mc,mainc,matchc,alpha,delta] = setParamValues(params);
if ~modeltofit.epsilon
    epsilon = fit_epsilon_opt(subj);
end

ntrials = sum(~isnan(display)&~isnan(realratings)); %length(stimuli);
llh = 0;
%simulate the model for one subject at a time
ratings = init*(ones(1,max(display))); %init for each subject
for trial = 1:ntrials
    torate = display(trial);
    if ~isnan(torate)
        rating = ratings(torate);
        action = realratings(trial);
        % get probability of that rating given the real one
        if ~isnan(action) %then no rating for that iteration of that task, learn but don't pick action
            prob = normpdf(action,rating,epsilon);
            if prob == 0 %this causes an inference error if you don't account for it
                prob = 1e-30;
            end
            llh = llh + log(prob);
        end
    end
    
    if modeltofit.delta & trial > 1
        costs = [uc mc mainc matchc];
        costs = setNewCosts(costs,delta,trial);
        uc = costs(1); mc = costs(2); mainc = costs(3); matchc = costs(4);
    end
    stim = stimuli(trial);
    updates = nupdates(trial);
    misses = nmisses(trial);
    maintenance = mains(trial);
    matches = nmatches(trial);
    cost = uc*updates + mc*misses + mainc*maintenance + matchc*matches; %add all the costs together
    if modeltofit.alpha
        ratings(stim) = ratings(stim) + alpha.*(cost-ratings(stim)); %delta rule
    else %alpha fixed or no alpha
        %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
        ratings(stim) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
    end
    
end %end of task run-through

negllh = -llh;
params;
subj;

end