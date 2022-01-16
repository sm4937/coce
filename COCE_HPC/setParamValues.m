function [uc,epsilon,init,missc,mainc,matchc,noisec,respc,lurec,errorc,fac,alpha,delta] = setParamValues(params,model)
%   Scale, set parameter values for cost learning models
global scalingvector canbeneg

scalingvector = ones(1,length(params)); canbeneg = false(1,length(params));
ncosts = sum(contains(model.paramnames,'c'));
cost_scaling = 50./ncosts; %works well at 20

epsilon = 10; init = 0; alpha = 0; delta = 0;
uc = 0; missc = 0; mainc = 0; matchc = 0; noisec = 0; 
respc = 0; lurec = 0; errorc = 0; fac = 0;
if model.alpha
    idx = find(contains(model.paramnames,'alpha'));
    canbeneg(idx) = false;
    alpha = params(idx)*scalingvector(idx);
end
if model.deltai || model.delta
    idx = find(contains(model.paramnames,'delta'));
    scalingvector(idx) = 1; %1; %0.05; 
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

end

