function [parameters] = applyTrans_parameters(model,parameters)
%applyTrans_parameters Apply transformations to normally distributed
%parameter values that CBM package requires
%   Transforms parameters which requires positive only values, or which are
%   bounded in some way (like learning rates)
if model.epsilon %standard deviation of ratings, can't be negative
    idx = find(contains(model.paramnames,'epsilon'));
    x = parameters(:,idx);
    parameters(:,idx) = exp(x); %transform with exp so no neg values, unbounded positive values
end
if model.alpha %learning rate, like from RL model - must be between 0 and 1
    idx = find(contains(model.paramnames,'alpha'));
    x = parameters(:,idx);
   parameters(:,idx) = 1./(1+exp(-x)); %transform with sigmoid
end

% these two delta parameters should be positive - then one of them is
% forced to be negative (delta_practice reduces costs of cognition)
% if model.delta_time
%     idx = find(contains(model.paramnames,'delta_time'));
%     x = parameters(:,idx);
%     parameters(:,idx) = 1./(1+exp(-x)); %transform 
% end
% if model.delta_practice
%     idx = find(contains(model.paramnames,'delta_practice'));
%     x = parameters(:,idx);
%     parameters(:,idx) = 1./(1+exp(-x)); %transform 
% end
% but recovery is awful, so trying a different approach


end

