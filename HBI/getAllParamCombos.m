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

if sum(contains(paramlist,{'delta','deltai'}))>0 %choose one update rule for now
    if sum(contains(paramlist,'deltai'))>0 %is it delta or deltai?
        % one delta for all, for one per cost parameter?
        name = 'epsilon_init_deltai';
    else
        name = 'epsilon_init_delta';
    end
    paramlist(contains(paramlist,{'delta','deltai'})) = []; %trim it here
else
    name = 'epsilon_init_alpha';
end

nparams = length(paramlist); model_specs = [];
for dim = 1:nparams
    model_specs = [model_specs; nchoosek([1:nparams],dim) NaN(length(nchoosek([1:nparams],dim)),nparams-dim)];
end

for m = 1:size(model_specs,1)
    spec = model_specs(m,:);
    modelname = name;
    for p = 1:sum(~isnan(spec))
        modelname = [modelname '_' paramlist{spec(p)}];
    end
    modellist{m} = modelname;
end

end

