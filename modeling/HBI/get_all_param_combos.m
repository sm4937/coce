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

name = 'epsilon_init';


% these used to be mutually exclusive (alpha vs delta models)
% but here I've removed that constraint. Now I just want to make sure they
% are listed first in the model name/model spec.
if sum(contains(paramlist,'alpha'))>0
    name = [name '_alpha'];
    paramlist(contains(paramlist,'alpha')) = [];
end

% delta models have one shared "fatigue" or "practice" parameter under the
% assumption these things act generally instead of 1 cost at a time. deltai
% models allow these costs to change separately, supposing maybe subjects
% tire of one cost more than the others, or something along those lines.
if sum(contains(paramlist,{'delta','deltai'}))>0 %choose one update rule for now
    if sum(contains(paramlist,'deltai'))>0 %is it delta or deltai?
        % one delta for all, or one per cost parameter?
        name = [name '_deltai'];
    elseif sum(contains(paramlist,{'delta_time','delta_practice'}))==2
        name = [name '_twodeltas'];
    elseif sum(contains(paramlist,'deltaexp'))>0
        name = [name '_deltaexp'];
    else
        name = [name '_delta'];
    end
    paramlist(contains(paramlist,{'delta'})) = []; %trim it here
end

% there's some form of init in all models at this point, so here I'm just
% making sure it's correctly designating init (1 init param) versus initi
% (1 init param per task (so 3)) models.
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

