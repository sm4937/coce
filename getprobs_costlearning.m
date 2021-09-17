function [map,negllh,llh] = getprobs_costlearning(params)
%Get minimum a priori score (MAP) for specific parameter values of cost learning model
%   Runs model with input value of parameters
%   Returns map and negllh

global onesim modeltofit mu sigma fit_epsilon_opt realsubjectsflag HBI_flag
if ~realsubjectsflag
% simulated data, not real
% simdata = [simdata; subj task rating torate updates misses mains matches];
    subj = unique(onesim(:,1));
    stimuli = onesim(:,2);
    realratings = onesim(:,3);
    display = onesim(:,4);
    nupdates = onesim(:,5);
    nmisses = onesim(:,6);
    mains = onesim(:,7);
    nmatches = onesim(:,8);
    noisiness = onesim(:,9);
    responses = onesim(:,10);
    nlures = onesim(:,11);
else %fitting real subject data
    subj = unique(onesim.subj);
    stimuli = onesim.task;
    realratings = (onesim.BDM); %scale down? just added this 11/08/2020 - see how this works
    display = onesim.display;
    if modeltofit.alpha
        nupdates = zeros(length(onesim.nupdates),1); nupdates(onesim.nupdates>0,:) = zscore(onesim.nupdates(onesim.nupdates>0,:)); % need to edit nupdates because it has so many zeros from irrelevant task 1
        nmisses = zscore(onesim.nmisses);
        mains = zscore(onesim.maintained); 
        nmatches = zscore(onesim.nmatches);
        noisiness = zscore(onesim.noisiness);
        responses = zscore(onesim.nresponses);
        nlures = zscore(onesim.nlures);
    elseif (modeltofit.delta || modeltofit.deltai)
        nupdates = onesim.nupdates; % need to edit nupdates because it has so many zeros from irrelevant task 1
        nmisses = onesim.nmisses;
        mains = onesim.maintained; 
        nmatches = onesim.nmatches;
        noisiness = onesim.noisiness;
        responses = onesim.nresponses;
        nlures = onesim.nlures;
    end
end

if HBI_flag
    params = applyTrans_parameters(modeltofit,params);
end
[uc,epsilon,init,mc,mainc,matchc,noisec,respc,lurec,alpha,delta] = setParamValues(params,modeltofit);
ntrials = sum(~isnan(realratings)); %length(stimuli);

costs = [uc mc mainc matchc noisec respc lurec];

%simulate the model for one subject at a time
ratings = [1 init].*(ones(1,max(display))); ratings_list = NaN(ntrials,1); %init for each subject 
costs = repmat(costs,ntrials,1);
if (modeltofit.delta || modeltofit.deltai)
    for trial = 2:ntrials
        costs(trial,:) = setNewCosts(costs(trial-1,:),delta,trial);
    end
end
components = [nupdates nmisses mains nmatches noisiness responses nlures];
cost = sum(costs.*components(1:ntrials,:),2); %add all the costs together
for trial = 1:ntrials
    torate = display(trial); stim = stimuli(trial);
    if ~isnan(torate) %skip those trials
        %ratings(ratings<0)=0;
        rating = ratings(torate); ratings_list(trial) = rating; %ratings(stim) 
    end %end of disqualifying data type for displayed task
    
    if modeltofit.alpha
        ratings(stim) = ratings(stim) + alpha.*(cost(trial,:)-ratings(stim)); %delta rule
    else %alpha fixed or no alpha
        %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
        %ratings(stim) = cost(trial); %no learning, no compounding, no delta rule. just a basic regression on last round
        ratings(stim) = ratings(stim)+ 1*(cost(trial,:)-ratings(stim));
    end
end
%calculate epsilon optimally
% epsilon_opt = sqrt((1/ntrials) * sum(realratings(~isnan(realratings)) - ratings_list).^2); fit_epsilon_opt(subj) = epsilon_opt; %track fit epsilons
% if ~model.epsilon; epsilon = epsilon_opt; end 
try
probs = normpdf(realratings(~isnan(realratings)),ratings_list,epsilon); probs(probs==0) = 1e-100; %can't be 0 exactly for llh
catch
    blah = true %what wrong?
end
negllh = -nansum(log(probs));
llh = nansum(log(probs)); map = llh;
if ~HBI_flag %run these calculations if running EM_fit
    priors = normpdf(params',mu',sigma');
    lprior = sum(log(priors));
    map = negllh-lprior;
end
end


