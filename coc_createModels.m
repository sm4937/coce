function [model] = coc_createModels(name)
%coc_createModels Create model specifications for cost of control
%   Outlines parameters, special functionality, other things that may
%   distinguish models from one another
%   For example, takes in 'alpha_epsilon' and outputs model structure with
%   number of parameters, and which ones they are.

modelStruct = struct;
possibleparams = {'uc','epsilon','init','initi','mc','mainc','alpha','delta','deltai','matchc','noisec','respc','lurec'};

params = strsplit(name,'_'); 
if length(params)<2
    params = strsplit(name,'-');
    %split by both in case you've mixed them up somewhere
end

for p = 1:length(possibleparams)
    if sum(contains(params,possibleparams{p}))>0
        eval(['model.' possibleparams{p} ' = true;'])
    else % if this parameter isn't in the model list, set it to false
        eval(['model.' possibleparams{p} ' = false;'])
    end
end 

if model.initi
    % If initi is in play, create one init parameter for each task
    model.init_1 = true;
    model.init_2 = true;
    model.init_3 = true;
    params = [params 'init_1' 'init_2' 'init_3'];
else
    model.init_1 = false;
    model.init_2 = false;
    model.init_3 = false;
end

% I'm implementing it this way so deltai is always last on the list of
% parameters
% This way, the extra delta parameters which are implied by the 'deltai'
% parameter are able to expand into the end of the total nparams and list
% of parameters when model fitting & simulating is trying to assign values
% to delta parameters which are not explicitly defined here

if sum(contains(params,'deltai'))>0
    model.deltai = true;
    model.delta = false; %this is what you'll use if all parameters 
    % have the same delta, which is false when deltai is in play
    params(contains(params,'deltai')) = [];
    % delete deltai entry
    %model is specified to have diff deltas for diff costs
    if model.uc
        model.deltaupdate = true;
        params{end+1} = 'deltaupdate';
    else
        model.deltaupdate = false;
    end
    if model.mc
        model.deltamiss = true;
        params{end+1} = 'deltamiss';
    else
        model.deltamiss = false;
    end
    if model.mainc
        model.deltamain = true;
        params{end+1} = 'deltamain';
    else
        model.deltamain = false;
    end
    if model.respc
        model.deltaresp = true;
        params{end+1} = 'deltaresp';
    else
        model.deltaresp = false;
    end
    if model.lurec
        model.deltalure = true;
        params{end+1} = 'deltalure';
    else
        model.deltalure = false;
    end
else
    model.deltaupdate = false;
    model.deltamiss = false;
    model.deltamain = false;
    model.deltaresp = false;
    model.deltalure = false;
end

model.nparams = length(params);
model.paramnames = params;

end

