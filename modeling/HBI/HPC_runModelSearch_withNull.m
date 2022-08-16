%% Script to run HBI model search on NYU HPC
% Specifically, want to run with null model as comparison on HPC

% Clean slate
realsubjectsflag = true; HBI_flag = true;

%add relevant paths
addpath('cbm-master/codes/');
addpath('./..'); %other COC code, like model loading
addpath('model-functions/')
load('./../simdata/toanalyze.mat');
% grab task characteristics from real subjects
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];
% Set up some plotting variables

% Dictate models to fit here by specifying which parameters are of interest
paramsofinterest = {'mainc','lurec','missc','fac','respc','initi'};
modelstofit = getAllParamCombos(paramsofinterest);
modelstofit = [modelstofit getAllParamCombos({'mainc','lurec','missc','fac','respc','deltai','initi'})];

% Run model fitting over reliable models only (i.e. over models where the
% fidelity of recovered parameters is significant)
for m = 1:length(modelstofit)
    file = load(['model-details/' modelstofit{m}]);
    recoverability(m) = file.reliable=='y';
end
modelstofit = modelstofit(recoverability);

% Go over existing models to ensure they have corresponding functions for
% their execution. In HBI framework, every model needs its own function.
% In my model fitting code, there is one fitting function which flexibly
% fits all models, so each model function is just named for the model, but
% contains the same call to that general model fitting function
% (getprobscostlearning.m).

function_folder = dir('model-functions'); list_existing = cellstr(string(char(function_folder.name)));
for m = 1:length(modelstofit)
    % now, cycling through models which are recoverable, create function
    % calls for running those models if they don't already exist
    diff = length(list_existing{end})-length(['fit_' modelstofit{m} '.m']);
    if sum(contains(list_existing,['fit_' modelstofit{m} '.m' repmat(' ',1,diff)]))==0 %no function for running this model, yet
        copyfile('dictate_model.m',['model-functions/fit_' modelstofit{m} '.m'])
    end
    %otherwise, do nothing
end

% Format real subject data for feeding in to the HBI framework. HBI is
% looking for 1 struct entry per subject.
nsubjs = length(unique(toanalyze.subj));
for subj = 1:length(unique(toanalyze.subj))
    data{subj} = toanalyze(toanalyze.subj==subj,:);
end

% Define the parameters/inputs to HBI, including MLE fits of each model
v = 6.25;
for m = 1:length(modelstofit)
    model_name = modelstofit{m}; modeltofit = coc_createModels(model_name);
    fnames{m} = [modelstofit{m} '.mat'];
    eval(['func = @fit_' modelstofit{m} ';']); funcs{m} = func;
    priors{m} = struct('mean',zeros(modeltofit.nparams,1),'variance',v); 
    model_labels{m} = strrep(modelstofit{m},'_','-');
    %if fitflag; cbm_lap(data, func, priors{m}, fnames{m}); end
end
fname_hbi = 'HBI_coc_48models_2022.mat'; 
% 63 models includes all alpha/delta combos
% 37 models includes all alpha/deltai combos (reduced cost space)
% 44 models includes all alpha/deltai combos (adding miss costs back in)
% 26 models includes new initi paramater, and excludes update costs entirely.
% 26models_compareMaintenance compares the old way of computing
% maintenance, with the new way (which is just main = n)
% 48models_2022 incorporates new cost schemes w new maintenance costs, multiple deltas,
% and multiple initial fair wages


% % RUN HBI %%

% Inputs are:
%data {cell per subj}, model-specific fitting functions, filenames from
%cbm_lap, %filename for saving full running to
%cbm_hbi(data,funcs,fnames,fname_hbi);
cbm_hbi_null(data,fname_hbi);
save('HPC_workspace.mat')
