function [costs] = set_new_costs(costs,delta,trial)
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

% Linear, additive cost changing scheme
costs = costs.*(1+(excdelta.*((trial-1)/ntrials)));
% costs on the right side is always costs_0
% (first costs, since they're modulated by trial, trial number does the
% adding up)


end

