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

if isfolder('../../data')
    
    % ALL experimental data, all 100 subjects live in 'data'
    prefix = '../../data/';
    % grab all subjects from those files, then
    
elseif isfolder('../../public_data')
    
    prefix = '../../public_data/';
    
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
paramsofinterest = {'uc','mainc','lurec','missc','fac','respc','alpha','initi'};

modelstofit = get_all_param_combos(paramsofinterest);
% Combine into one set of models

modelstosim = modelstofit;

function_folder = dir('model-functions'); list_existing = cellstr(string(char(function_folder.name)));
for m = 1:length(modelstofit)
    % now, cycling through models which are recoverable, create function
    % calls for running those models if they don't already exist
    diff = length(list_existing{end})-length(['fit_' modelstofit{m} '.m']);
    if sum(contains(list_existing,['fit_' modelstofit{m} '.m' repmat(' ',1,diff)]))==0 %no function for running this model, yet
        copyfile('dictate_model.m',['model-functions/fit_' modelstofit{m} '.m'])
    end
    %otherwise, do nothing
    eval(['funcs{m}{1} = @fit_' modelstofit{m} ';']); 
    eval(['fnames{m}{1} = "' modelstofit{m} '.mat";']); 
    fnames_typeIIML{m}{1} = fnames{m}{1};
end

subjnums = unique(toanalyze.subj);
nsubjs = length(subjnums); 

nsubjstofit = 30;
% can reduce model fitting space this way, make fitting faster, by reducing
% number of subjects to fit

if nsubjstofit<nsubjs
    repeated_subset = randperm(100,nsubjstofit);
else
    repeated_subset = 1:nsubjs;
end
subset = repeated_subset;


%Many of these generate/recover results already exist. By making forcefit
%true, you can re-recover the fits, if you'd like (like if you changed
%something in that model). Otherwise, leaving it false allows you to
%examine the gen/rec results without re-running fitting from scratch.
forcefit = true;

% modelstosim = {};
% modelstosim{1} = 'epsilon_initi_alpha_uc_mainc_fac';
% model 22 in run_model_fitting

data = cell([length(modelstofit),1]);
modelStruct = cell([length(modelstofit),1]);
priors = cell([length(modelstofit),1]);
fits = cell([length(modelstofit),1]);
fitparams_allmodels = cell([length(modelstofit),1]);
realparamlist_allmodels = cell([length(modelstofit),1]);

details_folder = dir('model-details'); list_existing = cellstr(string(char(details_folder.name)));
justfit_indices = [];

for m = 1:length(modelstofit)
   
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
    
    %have you run gen-rec on this model already? if no, don't run
    % unless forcefit is on, then re-run
    if sum(contains(list_existing,[model_name '.mat' repmat(' ',1,diff)]))==0 || forcefit 
        %if no, run it and save results
        for subj = 1:nsubjstofit 
            subjnum = repeated_subset(subj); %grab real subj number
            params = rand(1,modeltosim.nparams); 
            
            inits = contains(modeltosim.paramnames,'init');
            params(inits) = normrnd(0.5,0.2,1,sum(inits)); %real init spreads
            deltas = contains(modeltosim.paramnames,'delta');
            params(deltas) = normrnd(0.25,0.5,1,sum(deltas)); 
            %alphas = contains(modeltosim.paramnames,'alpha');
            %params(alphas) = normrnd(0.75,0.01,1,sum(alphas));  
            
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
        fname_hbi = ['genrec_1model.mat'];        
        cbm_hbi(subsetdata,funcs{m},fnames_typeIIML{m}{1},fname_hbi);
        %inputs: data {cell per subj}, model-specific fitting functions, filenames from
        %cbm_lap, %filename for saving full running to
        
        fits{m} = load(fname_hbi);
        modelStruct{m} = fits{m}.cbm;
        
        fitparams_allmodels{m} = modelStruct{m}.output.parameters{1};
        fitparams_allmodels{m} = applyTrans_parameters(modeltofit,fitparams_allmodels{m});
        
        justfit_indices = [justfit_indices m];
        
        reliable = 'm'; %m for maybe!
        realparamlist = realparamlist_allmodels{m}; fitparams = fitparams_allmodels{m};
        save(['model-details/' modelstosim{m} '.mat'],'realparamlist','fitparams','subset','reliable','data','subsetdata')
                     
    end

end


%% MAKE PLOTS OF FITS!
% now, do visual and numerical inspection of parameter fits
% set cost of interest to '' if you want to see all models 
param_of_interest = 'initi';
models_containing_param_of_interest = find(contains(modelstosim,param_of_interest));
%models_containing_param_of_interest = justfit_indices;

for m = 1:length(models_containing_param_of_interest)
    
    mii = models_containing_param_of_interest(m); %models_containing_cost_of_interest(m);
    
    figure(1)
    load(['model-details/' modelstosim{mii} '.mat'])
    
    model = coc_createModels(modelstosim{mii});
    %fitparams = fitparams_allmodels{mii};
    %realparamlist = realparamlist_allmodels{mii};
    
    if sum(contains(model.paramnames,'delta_time'))>0
        idx = find(contains(model.paramnames,'delta_time'));
        model.paramnames{idx} = 'delta_{time}';
        idx = find(contains(model.paramnames,'delta_practice'));
        model.paramnames{idx} = 'delta_{practice}';
    end

    nparams = size(fitparams,2);
    for p = 1:nparams
        subplot(5,3,p+sum(contains(model.paramnames{p},'delta')))
        % temporary: push deltas one row below in order to make them easier
        % to move around on the plot
        scatter(realparamlist(1:length(subset),p),fitparams(1:length(subset),p),[],rand(length(subset),3),'Filled')
        hold on
        plot([0 1],[0 1],'--')
        xlabel(['Real ' model.paramnames{p}])
        ylabel('Fit values')
        %xlim([0 1]); ylim([0 1]);
        ax = gca; ax.FontSize = 18;
    end
    fig = gcf; fig.Color = 'w';
    MSEs = mean((realparamlist(1:size(fitparams,2))-fitparams).^2);
    disp(['MSE = ' num2str(MSEs)])
    
    [rs,ps] = corr(fitparams,realparamlist);
    r = diag(rs)
    p = diag(ps)
    disp(['Param fits for ' modelstosim{mii}])
    % how reliable is the recovery? r's closer to 1 are better, closer
    % to 0 are worse
    
    disp(['Model ' num2str(mii)])
    reliable = input('does this model fit look reliable? y/n','s');
    save(['model-details/' modelstosim{mii} '.mat'],'realparamlist','fitparams','subset','reliable','data','subsetdata')
    close 1
    
end
