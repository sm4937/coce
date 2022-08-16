function [map,negllh] = map_costlearning(params)
%Get minimum a priori score (MAP) for specific parameter values of cost learning model
%   Runs model with input value of parameters
%   Returns map and negllh

global onesim modeltofit mu sigma fit_epsilon_opt noiselessri realsubjectsflag model

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
% real data format
% toanalyze = 1. subj num, 2. nmatches, 3. maintained, 4. nupdates
%   5. nmisses, 6. task displayed, 7. task executed, 8. BDM, 9. accurate
%   BDM for that task
else %fitting real subject data
    subj = unique(onesim.subj);
    stimuli = onesim.task;
    realratings = (onesim.BDM); %scale down? just added this 11/08/2020 - see how this works
    display = onesim.display;
    nupdates = zeros(length(onesim.nupdates),1); nupdates(onesim.nupdates>0,:) = zscore(onesim.nupdates(onesim.nupdates>0,:)); % need to edit nupdates because it has so many zeros from irrelevant task 1
    nmisses = zscore(onesim.nmisses); %./max([max(onesim.nmisses) 1]); %if nmisses = 0, divide by 1 instead 
    mains = zscore(onesim.maintained); %./max(onesim.maintained); 
    nmatches = zscore(onesim.nmatches);
    noisiness = zscore(onesim.noisiness);
end

model = modeltofit; %make agnostic variable model to send inputs in to param values func without overwriting other stuff outside this loop
[uc,epsilon,init,mc,mainc,matchc,noisec,alpha,delta] = setParamValues(params);
costs = [uc mc mainc matchc noisec];
ntrials = sum(~isnan(realratings)); %length(stimuli);

%simulate the model for one subject at a time
ratings = init*(ones(1,max(display))); ratings_list = NaN(ntrials,1); %init for each subject 
for trial = 1:ntrials
    torate = display(trial);
    if ~isnan(torate) %skip those trials
        %ratings(ratings<0)=0;
        rating = ratings(torate); ratings_list(trial) = rating; %ratings(stim) 
    end %end of disqualifying data type for displayed task
    
    if modeltofit.delta & trial > 1
        costs = setNewCosts(costs,delta,trial);
    end
    stim = stimuli(trial);
    cost = costs*[nupdates(trial);nmisses(trial);mains(trial);nmatches(trial);noisiness(trial)]; %add all the costs together
    if modeltofit.alpha
        ratings(stim) = ratings(stim) + alpha.*(cost-ratings(stim)); %delta rule
    else %alpha fixed or no alpha
        %ratings(stim) = ratings(stim) + alpha*(cost); %compounding cost model
        ratings(stim) = cost; %no learning, no compounding, no delta rule. just a basic regression on last round
    end
end %end of task run-through
%calculate epsilon optimally
epsilon_opt = sqrt((1/ntrials) * sum(realratings(~isnan(realratings)) - ratings_list).^2);
if ~modeltofit.epsilon; epsilon = epsilon_opt; end
probs = normpdf(realratings(~isnan(realratings)),ratings_list,epsilon); probs(probs==0) = 1e-100; %can't be 0 exactly for llh
negllh = -nansum(log(probs));
fit_epsilon_opt(subj) = epsilon_opt; %track fit epsilons

priors = normpdf(params',mu',sigma');
lprior = sum(log(priors));

map = negllh-lprior;
end


