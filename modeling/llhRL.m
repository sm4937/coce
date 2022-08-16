function [negllh] = llhRL(params)
%Get negative llh for specific parameter values of simple RL model
%   Runs model with input value of parameters
%   Returns llh

global onesim modeltofit
% simdata = [simdata; subj stim r action alpha];
stimuli = onesim(:,2);
realratings = onesim(:,3);
ntrials = length(stimuli);
llh = 0;
alpha = params(1); epsilon = 0.7; init = 1;

if modeltofit.epsilon
    epsilon = params(2);
end
if modeltofit.init
    init = params(3)*5;
end
costs = 1+stimuli; %change to learning about just one task

%simulate the model for one subject at a time
ratings = init*(ones(1,length(unique(stimuli)))); %init for each subject
for trial = 1:ntrials
    stim = stimuli(trial);
    cost = costs(trial);
    action = realratings(trial);
    ratings(stim) = ratings(stim) + 0.3*(alpha*cost); %not a learning rate of one
    prob = normpdf(action,ratings(stim),epsilon);
    if prob == 0 %this causes an inference error if you don't account for it
        llh = llh+0;
    else
        llh = llh + log(prob);
    end
end

negllh = -llh;

% Q = init*ones(1,length(unique(stimuli))); %init for each subject
% for trial = 1:ntrials
%     stim = stimuli(trial);
%     action = actions(trial);
%     %llh = llh + log(Q(action,stim));
%     %r = action==stim;
%     %Q(action,stim) = Q(action,stim) + alpha*(r-Q(action,stim));
%     stim = stimuli(trial);
%     cost = costs(trial);
%     r = rand()*5;
%     r = r - alpha*cost; %cost
%     Q(stim) = Q(stim) + 0.3*(r-Q(stim));
% end

end


