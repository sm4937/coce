function [modellist] = getAllParamCombos(paramlist)
%getAllParamCombos Get all possible model combinations and their
%recoverability
%   paramlist is a cellstr list of possble param names, like
%   {'mc','mainc','lurec','uc'...}
%   modellist is all possible models combining those params, some
%   standalone like "mainc_epsilon_init_alpha", others like
%   "mainc_lurec_epsilon_init_alpha"
%   for COCE, there are three mandatory parameters: epsilon, init, and
%   alpha

nparams = length(paramlist); model_specs = [];
for dim = 1:nparams
    model_specs = [model_specs; nchoosek([1:nparams],dim) NaN(length(nchoosek([1:nparams],dim)),nparams-dim)];
end

for m = 1:size(model_specs,1)
    spec = model_specs(m,:);
    name = 'epsilon_init_alpha';
    for p = 1:sum(~isnan(spec))
        name = [name '_' paramlist{spec(p)}];
    end
    modellist{m} = name;
end

end

