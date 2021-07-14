function [model] = coc_createModels(name)
%coc_createModels Create model specifications for cost of control
%   Outlines parameters, special functionality, other things that may
%   distinguish models from one another
%   For example, takes in 'alpha_epsilon' and outputs model structure with
%   number of parameters, and which ones they are.

modelStruct = struct;
possibleparams = {'uc','epsilon','init','mc','mainc','alpha','delta','matchc','noisec','respc','lurec'};

params = strsplit(name,'_'); params = strsplit(name,'-'); %split by both in case you've mixed them up somewhere
model.nparams = length(params);
for p = 1:length(possibleparams)
    if sum(contains(params,possibleparams{p}))>0
        eval(['model.' possibleparams{p} ' = true;'])
    else
        eval(['model.' possibleparams{p} ' = false;'])
    end
end
model.paramnames = params;

end

