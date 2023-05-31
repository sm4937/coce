%% Proof of model recoverability

% Runs a generate and recover procedure such that model data is simulated
% with set parameters, and then those same parameter values are
% recovered. Gives a sense of the fidelity of model recoverability, which
% is a concern in the Wagers for Work paper because of the low number of
% trials we're fitting per subject (32), and the fact that some of our cost
% components are correlated with each other. Want to ensure accurate
% parameter recovery.

% If you'd like, you can run this script such that the only models you
% use/test in run_model_fitting are recoverable (like I did in the manuscript)
% Requires visual/correlation inspection of each model's recoverability

global realsubjectsflag HBI_flag
realsubjectsflag = false; HBI_flag = true;

%add relevant paths
addpath('../cbm-master/codes/');
addpath('../model-functions/');
addpath(genpath('../../')); %other COC code, like model loading

%if isdir('data')
if isdir('../../../data')
    % ALL experimental data, all 100 subjects live in 'data'
    load('../../../data/filenames.mat')
    prefix = '../../../data/';
    % grab all subjects from those files, then
elseif isdir('../../../example_data')
    files{1} = '../../../example_data/example_subjs.mat';
    prefix = '../../../example_data/';
end
% version 4

load([prefix 'toanalyze.mat'])
% Set up for plotting 
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];

% WHICH PARAMETER COMBINATIONS DO YOU WANT?
% Models are specified by which parameters you want to include, like so:
% "epsilon_initi_alpha_mainc"
% This is a model containing a noise term epsilon (they all do), a separate
% initial rating for each task initi, a learning rate alpha which dictates
% the speed of cost learning, and a maintenance cost only (most frequently
% recovered model in the full dataset of 100 subjects)

% Dictate models to fit here by specifying which parameters are of interest
paramsofinterest = {'uc','mainc','lurec','missc','fac','respc','deltai'};

modelstofit = get_all_param_combos(paramsofinterest);
% Combine into one set of models
% Add delta models into the mix by specifying the parameter 'deltai', which
% provides a different delta for each cost, or the parameter 'delta' which
% provides one delta for all costs.

modelstosim = modelstofit;

function_folder = dir('../model-functions'); list_existing = cellstr(string(char(function_folder.name)));
for m = 1:length(modelstofit)
    % now, cycling through models which are recoverable, create function
    % calls for running those models if they don't already exist
    diff = length(list_existing{end})-length(['fit_' modelstofit{m} '.m']);
    if sum(contains(list_existing,['fit_' modelstofit{m} '.m' repmat(' ',1,diff)]))==0 %no function for running this model, yet
        copyfile('../dictate_model.m',['../model-functions/fit_' modelstofit{m} '.m'])
    end
    %otherwise, do nothing
    eval(['funcs{m}{1} = @fit_' modelstofit{m} ';']); 
    eval(['fnames{m}{1} = "' modelstofit{m} '.mat";']); 
    fnames_typeIIML{m}{1} = fnames{m}{1};
end

subjnums = unique(toanalyze.subj);
nsubjs = length(subjnums); 

nsubjstofit = 50;
% can reduce model fitting space this way, make fitting faster, by reducing
% number of subjects to fit

if nsubjstofit<nsubjs
    repeated_subset = randperm(100,nsubjstofit);
else
    repeated_subset = 1:nsubjs;
end


%Many of these generate/recover results already exist. By making forcefit
%true, you can re-recover the fits, if you'd like (like if you changed
%something in that model). Otherwise, leaving it false allows you to
%examine the gen/rec results without re-running fitting from scratch.
forcefit = false;

% modelstosim = {};
% modelstosim{1} = 'epsilon_initi_alpha_uc_mainc_fac';
% model 22 in run_model_fitting

data = cell([length(modelstofit),1]);
modelStruct = cell([length(modelstofit),1]);
priors = cell([length(modelstofit),1]);
fits = cell([length(modelstofit),1]);
fitparams_allmodels = cell([length(modelstofit),1]);
realparamlist_allmodels = cell([length(modelstofit),1]);

fits_folder = dir('../model-details'); list_existing = cellstr(string(char(fits_folder.name)));
just_fit_list = [];
for m = 1:length(modelstosim)
    
    data{m} = cell([nsubjstofit,1]);
    fits{m} = struct;
    modelStruct{m} = struct;
    fitparams_allmodels{m} = [0];
    realparamlist_allmodels{m} = [0];
    
    model_name = modelstosim{m};
    modeltosim = coc_createModels(model_name); 
    modeltofit = modeltosim;
    diff = length(list_existing{end})-length([model_name '.mat']);
    
    subsetdata = cell(1,nsubjstofit);
    
    if sum(contains(list_existing,[model_name '.mat' repmat(' ',1,diff)]))==0 || forcefit %have you run gen-rec on this model already?
        %if no, run it and save results
        for subj = 1:nsubjstofit 
            subjnum = repeated_subset(subj); %grab real subj number
            params = rand(1,modeltosim.nparams); 
            
            inits = contains(modeltosim.paramnames,'init');
            params(inits) = normrnd(0.5,0.2,1,sum(inits)); %real init spreads
            deltas = contains(modeltosim.paramnames,'delta');
            params(deltas) = normrnd(0.03,1,1,sum(deltas));
            
            onesubj = toanalyze(toanalyze.subj==subjnum,:);
            data{m}{subj} = simulate_cost_model(modeltosim,params,onesubj);
            %names = modeltofit.paramnames;
            %params(contains(names,'delta')) = 0.1*params(contains(names,'delta'));
            %params(contains(names,'init_2')) = normrnd(0.7,0.3); params(contains(names,'c')) = normrnd(1,0.5);
            realparamlist_allmodels{m}(subj,1:length(params)) = params;
        end
        
        %plot simulated dataset to see whether it contains sensical values
        %model_validation_HBI()
        
        priors{m}{1} = struct('mean',zeros(modeltofit.nparams,1),'variance',6.25); 
        
        for sii = 1:length(repeated_subset) %%use subset for this part, to check MLE fits
            subsetdata{sii} = data{m}{sii};
        end
        cbm_lap(subsetdata, funcs{m}{1}, priors{m}{1}, fnames{m}{1});
        
        % Run the full hierarchical fitting and test how it does on recovering true simulated models
        fname_hbi = ['genrec_model' num2str(m) '.mat'];        
        cbm_hbi(subsetdata,funcs{m},fnames_typeIIML{m}{1},fname_hbi);
        %inputs: data {cell per subj}, model-specific fitting functions, filenames from
        %cbm_lap, %filename for saving full running to
        
        fits{m} = load(fname_hbi);
        modelStruct{m} = fits{m}.cbm;
        
        fitparams_allmodels{m} = modelStruct{m}.output.parameters{1};
        fitparams_allmodels{m} = applyTrans_parameters(modeltofit,fitparams_allmodels{m});
        
        just_fit_list = [just_fit_list m];
             
    end

end

%% now, do visual and numerical inspection of parameter fits

for mii = 1:length(just_fit_list)
    figure(1)
    
    m = just_fit_list(mii);
    
    %load(['../model-details/' modelstosim{m} '.mat'])
    
    model = coc_createModels(modelstosim{m});
    fitparams = fitparams_allmodels{m};
    realparamlist = realparamlist_allmodels{m};

    nparams = size(fitparams,2);
    for p = 1:nparams
        subplot(5,3,p)
        scatter(realparamlist(:,p),fitparams(:,p),[],rand(nsubjstofit,3),'Filled')
        hold on
        plot([0 1],[0 1],'--')
        xlabel(['Real ' model.paramnames{p}])
        ylabel('Fit values')
        %xlim([0 1]); ylim([0 1]);
    end
    fig = gcf; fig.Color = 'w';
    MSEs = mean((realparamlist(1:size(fitparams,2))-fitparams).^2);
    disp(['MSE = ' num2str(MSEs)])
    
    [rs,ps] = corr(fitparams,realparamlist);
    r = diag(rs)
    p = diag(ps)
    disp(['Param fits for ' modelstosim{m}])
    % how reliable is the recovery? r's closer to 1 are better, closer
    % to 0 are worse
    
    disp(['Model ' num2str(m)])
    reliable = input('does this model fit look reliable? y/n','s'); close 1
    subset = repeated_subset;
    save(['../model-details/' modelstosim{m} '.mat'],'realparamlist','fitparams','subset','reliable','data','subsetdata')
    
end

%% Now, add in some hierarchical model assignments, re-simulate, test that
% portion of the HBI

true_models = [1; 2; 4; 1; 1; 1; 2; 2; 4; 2];

realparamlist = nan(nsubjstofit,8); 
for subj = 1:nsubjstofit
    subjnum = repeated_subset(subj); %grab real subj number
    modeltosim = coc_createModels(modelstosim{true_models(subj)});
    params = rand(1,modeltosim.nparams);
    onesubj = toanalyze(toanalyze.subj==subjnum,:);
    data{subj} = simulate_cost_model(modeltosim,params,onesubj);
    realparamlist(subj,1:length(params)) = params;
end

v = 6.25;
for m = 1:length(unique(true_models))
    model_name = modelstofit{true_models(m)};
    modeltofit = coc_createModels(model_name);
    fnames{1} = [model_name '.mat']; 
    eval(['funcs{1} = @fit_' model_name ';']); %I'm not sure how to do this in a generalizable way
    % Right now, this doesn't work because it's expecting a different
    % function for every model 
    % But with 31 possible models... I can't hand-write a function for each
    priors{1} = struct('mean',zeros(modeltofit.nparams,1),'variance',v);
    realparamlist = nan(nsubjs,8); 
    for subj = 1:nsubjs
        subjnum = subjnums(subj); %grab real subj number
        params = rand(1,modeltofit.nparams);
        onesubj = toanalyze(toanalyze.subj==subjnum,:);
        data{subj} = simulate_cost_model(modeltofit,params,onesubj);
        realparamlist(subj,1:length(params)) = params;
    end
    cbm_lap(data, funcs{1}, priors{1}, fnames{1});

    % Run the full hierarchical fitting and test how it does on recovering true simulated models
    fname_hbi = 'HBI_coc_lessrecoverablemodels.mat';
    %data {cell per subj}, model-specific fitting functions, filenames from
    %cbm_lap, %filename for saving full running to
    cbm_hbi(data,funcs,fnames,fname_hbi);

    % Analyze fit hierarchical generate/recover
    fits = load(fname_hbi);
    cbm   = fits.cbm;
    freqs = cbm.output.model_frequency;

    figure(1)
    fitparams = cbm.output.parameters{1};
    fitparams = applyTrans_parameters(modeltofit,fitparams);
    nparams = size(fitparams,2);
    for p = 1:nparams
        subplot(4,2,p)
        scatter(realparamlist(:,p),fitparams(:,p),'Filled')
        hold on
        plot([0 0],[1 1],'--')
        xlabel(['Real ' modeltofit.paramnames{p}])
        ylabel('Fit values')
        %xlim([0 1])
    end
    fig = gcf; fig.Color = 'w';
    MSEs = mean((realparamlist(1:size(fitparams,2))-fitparams).^2);
    disp(['MSE = ' num2str(MSEs)])
    
    [r,ps] = corr(fitparams,realparamlist);
    disp(['Param fits for ' model_name])  
    
    reliable = input('does this model fit look reliable? y/n','s'); close 1
    save(['model-details/' model_name],'realparamlist','fitparams','reliable')
end

%% 

model_inspection(m,cbm)

