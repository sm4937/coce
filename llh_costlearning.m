function [negllh] = llh_costlearning(params)
%Get negative llh for specific parameter values of simple RL model
%   Runs model with input value of parameters
%   Returns llh

global onesim modeltofit realsubjectsflag bounds
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
    realratings = (onesim.newBDM)./bounds(2); %scale down? added 11/08
    nupdates = onesim.nupdates;%./(max(onesim.nupdates));
    nmisses = onesim.nmisses;%./(max(onesim.nmisses)); %normalize
    mains = onesim.maintained;%./(max(onesim.maintained)); 
end

ntrials = length(stimuli);
llh = 0;
init = 0.25.*bounds(2); epsilon = 0.5.*bounds(2); mc = 0; mainc = 0; alpha = 0.35;
nparams = length(params);
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
    else %alpha is fixed
        ratings(stim) = ratings(stim) + alpha*(cost); %not a learning rate of one
    end
    if ~isnan(action) %then no rating for that iteration of that task, learn but don't pick action
    prob = normpdf(action,ratings(stim),epsilon);
    %prob = normcdf(action-0.01:action:action+0.01,ratings(stim),epsilon);
    if prob == 0 %this causes an inference error if you don't account for it
        prob = 0.00000001;
    end
    llh = llh + log(prob);
    end
end

negllh = -llh;

%dbstop in llh_costlearning.m at 57 if (negllh < 0)

end


