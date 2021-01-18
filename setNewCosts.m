function [costs] = setNewCosts(costs,delta,trial)
%setNewCosts Use delta to set new operation costs
% Takes costs in, returns new costs
%update parameter values by delta, decrement, increment

count = trial-1;
%costs = costs + (1/count).*delta;
%costs = costs*exp(-delta*count);
costs = costs*(1+(delta*count));

end

