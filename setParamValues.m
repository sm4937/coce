function [uc,epsilon,init,mc,mainc,matchc,noisec,respc,lurec,alpha,delta] = setParamValues(params,model)
%   Scale, set parameter values for cost learning models
global scalingvector canbeneg

labels = fieldnames(model); not_params = sum(contains(labels,'paramnames')|contains(labels,'nparams'));
cost_scaling = 20; %put at 1 for testing HBI_coc code, broke
scalingvector = ones(1,length(params)); canbeneg = false(1,length(params));

epsilon = 10; init = 0; alpha = 0; delta = 0;
uc = 0; mc = 0; mainc = 0; matchc = 0; noisec = 0; respc = 0; lurec = 0;
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
if model.init
    idx = find(contains(model.paramnames,'init'));
    scalingvector(idx) = 100;
    init = params(idx)*scalingvector(idx);
end
if model.mc
    idx = find(contains(model.paramnames,'mc'));
    scalingvector(idx) = cost_scaling; canbeneg(idx) = true;
    mc = params(idx)*scalingvector(idx);
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
if model.alpha
    idx = find(contains(model.paramnames,'alpha'));
    alpha = params(idx)*scalingvector(idx);
end
if model.delta
    idx = find(contains(model.paramnames,'delta'));
    idx = idx(1);
    scalingvector(idx) = 0.20;
    canbeneg(idx) = true;
    delta = params(idx)*scalingvector(idx);
end
if model.deltai
    idx = find(contains(model.paramnames,'delta'));
    %because of the way code upstream of this has been written
    % (model definition code) deltai should always be at 
    % the end of the list of params
    % such that any additional space at the end of a vector of 
    % numbers belongs to the extra unnamed delta parameters
    scalingvector(idx) = 0.20;
    canbeneg(idx) = true;
    delta = params(idx).*scalingvector(idx);
    % delta should now be a vector, not a single number
end

end

