%% Execute HBI model fitting through NYU HPC

%% % % FIT SUBJECTS! % % %%

clear all
% Clean slate
global realsubjectsflag HBI_flag
realsubjectsflag = true; HBI_flag = true;

%add relevant paths
addpath('cbm-master/codes/');addpath('model-functions/');
addpath('./..'); %other COC code, like model loading
load('./../simdata/toanalyze.mat');
% grab task characteristics from real subjects

% Dictate models to fit here by specifying which parameters are of interest
paramsofinterest = {'mainc','lurec','missc','fac','respc','initi'};
modelstofit = getAllParamCombos(paramsofinterest);
modelstofit = [modelstofit getAllParamCombos({'mainc','lurec','missc','fac','respc','deltai','initi'})];
% GET ALL POSSIBLE PARAM COMBOS
modelstofit(~contains(modelstofit,'c')) = [];
modelstosim = modelstofit;

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
    cbm_lap(data, func, priors{m}, fnames{m});
end
fname_hbi = 'HBI_HPC_errorcostscheme.mat'; 
% 63 models includes all alpha/delta combos
% 37 models includes all alpha/deltai combos (reduced cost space)
% 44 models includes all alpha/deltai combos (adding miss costs back in)
% 26 models includes new initi paramater, and excludes update costs entirely.
% 26models_compareMaintenance compares the old way of computing
% maintenance, with the new way (which is just main = n)
% HBI_HPC_errorcostscheme includes new miss AND false alarm costs

% % RUN HBI %%

% Inputs are:
%data {cell per subj}, model-specific fitting functions, filenames from
%cbm_lap, %filename for saving full running to
cbm_hbi(data,funcs,fnames,fname_hbi);

% Model fits get saved in a structure called cbm, load that up here and
% grab variables of interest from it.
fits = load(fname_hbi);
cbm   = fits.cbm;
freqs = cbm.output.model_frequency;

% Use CBM toolbox to print parameter means etc.
% 1st input is the file-address of the file saved by cbm_hbi
% 2nd input: a cell input containing model names
% 3rd input: another cell input containing parameter names of the winning model
[~,best] = max(cbm.output.exceedance_prob);
best_model = coc_createModels(modelstofit{best});
for p = 1:length(best_model.paramnames)
    original_name = best_model.paramnames{p};
    param_names{p}= original_name; transform{p} = 'none';
    if strmatch(original_name,'alpha')
        transform{p} = 'sigmoid';
        param_names{p} = ['\' original_name];
    end
    if strmatch(original_name,'epsilon')
        transform{p} = 'exp';
        param_names{p} = ['\' original_name];
    end
    if contains(original_name,'c')
        param_names{p} = strrep(original_name,'c',' cost');
        costs(p) = true;
    else
        costs(p) = false;
    end
end

% Here I create my own structure, best_model, for saving other measures of
% interest, like the name of every model that ran in this model fitting
% regime. This gets called in other analysis scripts.
best_model.lowparams = cbm.output.parameters{best};
best_model.highparams = cbm.output.group_mean{best}; best_model.overallfit = cbm.output;
best_model.overallfit.fitmodels = model_labels; best_model.name = model_labels{best};
modelStruct.best_model = best_model; save('HBI_modelStruct.mat','modelStruct');
