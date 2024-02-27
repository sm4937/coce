function [costs] = set_new_costs(costs,delta,trial,quadratic_delta_flag,two_delta_flag,time,iter)
%setNewCosts Use delta to set new operation costs
% Takes costs in, returns new costs
% update parameter values by delta, decrement, increment
ntrials = 32;

excdelta = zeros(1,size(costs,2));
excdelta(costs~=0) = 1;
% executable delta vector, which is 1's where the costs are 0 (no change in
% costs away from 0)
% and the true delta value elsewhere
if length(delta)>1
    excdelta(costs~=0) = excdelta(costs~=0).*delta(1:sum(costs~=0));
    % several deltas, one per cost
else
    excdelta(costs~=0) = excdelta(costs~=0).*delta;
    % one delta, many costs
end
% the executable delta vector of the proper size

if quadratic_delta_flag
    % exponential function of fatigue
    costs = costs.*(trial.^delta);
    % will keep costs at 0 still at 0
    
elseif two_delta_flag
    % first, increment costs according to time spent on the experiment at
    % time of rating
    costs = costs.*(1+(delta(1)*time));
    
    %second, decrement costs according to task practice at time of rating
    costs = costs.*(1-(delta(2)*iter));
    
    % both "time" and "iter" are percentages. either of total time spent on
    % the experiment at time of rating, or total iterations of that task
    % completed at time of rating.
    
    % they are also restricted to be percentages between 0 and 1
    
else
    % Linear, additive cost changing scheme
    % Original equation from bioRxiv preprint (july 2023)
    costs = costs.*(1+(excdelta.*((trial-1)/ntrials)));
    % costs on the right side is always costs_0
    % (first costs, since they're modulated by trial, trial number does the
    % adding up)
    
end


end

