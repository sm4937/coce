function [uc,epsilon,init,missc,mainc,matchc,noisec,respc,lurec,errorc,fac,rtc,otherlurec,alpha,delta] = set_param_values(params,model)
%   Scale, set parameter values for cost adaptation models
global scalingvector canbeneg

scalingvector = ones(1,length(params)); canbeneg = false(1,length(params));
ncosts = sum(contains(model.paramnames,'c'));
cost_scaling = 50./ncosts; %works well at 20 also
% i chose this cost scaling scalar because it best reproduced the true
% scales of the fair wage ratings & how they evolved over time (determined
% through simulation, which happens in simulate_cost_model)

epsilon = 10; init = 0; alpha = 0; delta = 0;
uc = 0; missc = 0; mainc = 0; matchc = 0; noisec = 0; 
respc = 0; lurec = 0; errorc = 0; fac = 0;

% new cost components suggested by reviewers for PLoS CB
rtc = 0; otherlurec = 0;

% initialize everything at 0, epsilon at 10 (because epsilon = 0 doesn't
% work because it's a standard deviation parameter), and below, according
% to model specification (contained in the "model" variable) set parameters
% to the values we'll use for simulation/fitting

if model.alpha
    idx = find(contains(model.paramnames,'alpha'));
    canbeneg(idx) = false;
    alpha = params(idx)*scalingvector(idx);
end
% the exponent delta needs to be little in order to get fit with reasomable
% values
if sum(contains(model.paramnames,'deltaexp')) > 0
    idx = find(contains(model.paramnames,'deltaexp'));
    scalingvector(idx) = 0.5; %0.5 in bioRxiv models
    canbeneg(idx) = true;
    delta = params(idx).*scalingvector(idx);
elseif sum(contains(model.paramnames,'delta')) > 0
    idx = find(contains(model.paramnames,'delta'));
    scalingvector(idx) = 0.5; %0.5 in bioRxiv models
    % 1 in supplementary alpha-delta models
    canbeneg(idx) = true;
    delta = params(idx).*scalingvector(idx);
    % delta should now be a vector, not a single number
end

if model.uc
    idx = find(contains(model.paramnames,'uc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    uc = params(idx)*scalingvector(idx);
end
if model.epsilon
    idx = find(contains(model.paramnames,'epsilon'));
    scalingvector(idx) = cost_scaling;
    epsilon = params(idx)*scalingvector(idx);
    if epsilon == 0; epsilon = 0.00000001; end
end
if model.init || model.initi
    idx = find(contains(model.paramnames,'init'));
    scalingvector(idx) = 100;
    init = params(idx).*scalingvector(idx);
    if length(idx) == 1
        init = repmat(init,1,3);
    end
end
if model.missc
    idx = find(contains(model.paramnames,'missc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    missc = params(idx)*scalingvector(idx);
end
if model.mainc
    idx = find(contains(model.paramnames,'mainc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    mainc = params(idx)*scalingvector(idx);
end
if model.matchc
    idx = find(contains(model.paramnames,'matchc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    matchc = params(idx)*scalingvector(idx);
end
if model.noisec
    idx = find(contains(model.paramnames,'noisec'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    noisec = params(idx)*scalingvector(idx);
end
if model.respc
    idx = find(contains(model.paramnames,'respc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    respc = params(idx)*scalingvector(idx);
end
if model.lurec
    idx = find(contains(model.paramnames,'lurec'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    lurec = params(idx)*scalingvector(idx);
end
if model.fac
    idx = find(contains(model.paramnames,'fac'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    fac = params(idx)*scalingvector(idx);
end
if model.errorc
    idx = find(contains(model.paramnames,'errorc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    errorc = params(idx)*scalingvector(idx);
end
if model.rtc
    idx = find(contains(model.paramnames,'rtc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    rtc = params(idx)*scalingvector(idx);
end
if model.otherlurec
    idx = find(contains(model.paramnames,'otherlurec'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    otherlurec = params(idx)*scalingvector(idx);
end

end

