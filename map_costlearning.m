function [map] = map_costlearning(params)
%Get negative llh for specific parameter values of simple RL model
%   Runs model with input value of parameters
%   Returns llh

global onesim modeltofit mu sigma fit_epsilon_opt noiselessri realsubjectsflag bounds

if ~realsubjectsflag
% simulated data, not real
% simdata = [simdata; subj stim r action alpha];
    subj = unique(onesim(:,1));
    stimuli = onesim(:,2);
    realratings = onesim(:,3);
    nupdates = onesim(:,4);
    nmisses = onesim(:,6);
    mains = onesim(:,8);
% real data format
% toanalyze = 1. subj num, 2. nmatches, 3. maintained, 4. nupdates
%   5. nmisses, 6. task displayed, 7. task executed, 8. BDM, 9. accurate
%   BDM for that task
else %fitting real subject data
    subj = unique(onesim.subj);
    stimuli = onesim.task;
    realratings = (onesim.newBDM)./bounds(2); %scale down? just added this 11/08/2020 - see how this works
    nupdates = onesim.nupdates;
    nmisses = onesim.nmisses; %normalize
    mains = onesim.maintained; 
end

ntrials = length(stimuli);
llh = 0;
init = 0.25.*bounds(2); epsilon = 0.5.*bounds(2); mc = 0; mainc = 0; alpha = 0.35;
nparams = modeltofit.nparams;
if modeltofit.uc
    idx = find(contains(modeltofit.paramnames,'uc'));
    uc = (params(idx));%.*bounds(2);
end
if modeltofit.epsilon
    idx = find(contains(modeltofit.paramnames,'epsilon'));
    epsilon = (params(idx)).*bounds(2);
%     epsilon = 1./(1+exp(-params(idx))); 
    if epsilon == 0; epsilon = 0.00000001; end
end
if modeltofit.init
    idx = find(contains(modeltofit.paramnames,'init'));
    init = (params(idx)).*bounds(2);
end
if modeltofit.mc
    idx = find(contains(modeltofit.paramnames,'mc'));
    mc = (params(idx)).*bounds(2);
end
if modeltofit.mainc 
    idx = find(contains(modeltofit.paramnames,'mainc'));
    mainc = (params(idx)).*bounds(2);
end
if modeltofit.alpha
    idx = find(contains(modeltofit.paramnames,'alpha'));
    alpha = params(idx);
end
%costs = 1+stimuli; %change to learning about just one task

%simulate the model for one subject at a time
ratings = init*(ones(1,max(stimuli))); %init for each subject
for trial = 1:ntrials
    stim = stimuli(trial);
    updates = nupdates(trial);
    misses = nmisses(trial);
    maintenance = mains(trial);
    action = realratings(trial);
    cost = uc*updates;
    if modeltofit.mc
        cost = cost + mc*misses;
    end
    if modeltofit.mainc
        cost = cost + mainc*maintenance;
    end
    if modeltofit.alpha
        ratings(stim) = ratings(stim) + alpha*(cost-ratings(stim)); %delta rule
    else %alpha fixed, compounding cost model
        ratings(stim) = ratings(stim) + alpha*(cost); %not a learning rate of one
    end
    ratings_list(trial) = ratings(stim); 
    epsilon_opt = sqrt((1/trial) * sum(realratings(1:trial) - ratings_list(1:trial)').^2);
    epsilon_evo(trial) = epsilon_opt;
    if ~modeltofit.epsilon
        epsilon = epsilon_opt;
    end
    %dbstop in map_costlearning.m at 74 if isnan(action)
    if ~isnan(action) %then no rating for that iteration of that task, learn but don't pick action
    prob = normpdf(action,ratings(stim),epsilon);
    %prob = normcdf(action-0.01:action:action+0.01,ratings(stim),epsilon);
    if prob == 0 %this causes an inference error if you don't account for it
        prob = 0.00000001;
    end
    llh = llh + log(prob);
    end
end %end of task run-through
%[ratings_list] 
%[noiselessri]
fit_epsilon_opt(subj) = epsilon_opt; %track fit epsilons

negllh = -llh;
lprior = 0;
dbstop in map_costlearning.m at 105 if isinf(lprior)
for p = 1:nparams
    prob = normpdf(params(p),mu(p),sigma(p));
    lprior = lprior + log(prob);
end
map = (negllh)-lprior;
%dbstop in llh_costlearning.m at 57 if (negllh < 0)
plot_flag = false;
if plot_flag
    figure
    plot(1:ntrials,epsilon_evo)
    ylabel('epsilon')
    xlabel('trial')
end

end


