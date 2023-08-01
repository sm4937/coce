function [modellist] = get_all_param_combos(paramlist)
%getAllParamCombos Get all possible model combinations and their
%recoverability
%   paramlist is a cellstr list of possble param names, like
%   {'mc','mainc','lurec','uc'...}
%   modellist is all possible models combining those params, some
%   standalone like "mainc_epsilon_init_alpha", others like
%   "mainc_lurec_epsilon_init_alpha"
%   for COCE, there is only one mandatory parameter: epsilon
%   all models have one update parameter, either delta, deltai, or alpha
%   all models also have one initialization parameter, init or initi

if sum(contains(paramlist,{'delta','deltai'}))>0 %choose one update rule for now
    if sum(contains(paramlist,'deltai'))>0 %is it delta or deltai?
        % one delta for all, or one per cost parameter?
        name = 'epsilon_init_deltai';
    else
        name = 'epsilon_init_delta';
    end
    paramlist(contains(paramlist,{'delta','deltai'})) = []; %trim it here
else
    name = 'epsilon_init_alpha';
end

if sum(contains(paramlist,'initi'))>0
    name = strrep(name,'init','initi');
    paramlist(contains(paramlist,'initi')) = []; %trim it here
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

