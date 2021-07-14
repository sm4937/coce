%% An outer level script for running Piray & Daw's HBI package (cbm)
% Running it on subjects from cost of control data set
% After running with type II ML package (EMfit_sm), thought I'd like to
% account for the fact that some subjects are better fit by different
% models.
% Instructions for doing all this here: https://payampiray.github.io/cbm

% % Starting out with some simple generate/recover on my models of interest
clear all

global realsubjectsflag model_name HBI_flag
realsubjectsflag = false;

%add relevant paths
addpath('cbm-master/codes/');addpath('model-functions/');
addpath('./..'); %other COC code, like model loading
load('./../simdata/toanalyze.mat'); %toanalyze already excluding subjects who fit poorly
model_detail_folder = dir('model-details'); list_existing = cellstr(string(char(model_detail_folder.name)));
% grab task characteristics from real subjects

% Run a test to ensure that individual parameters are being fit reasonably
function_folder = dir('model-functions'); list_existing = cellstr(string(char(function_folder.name))); %list of existing model functions
% WHICH PARAMS DO YOU WANT YOUR MODEL TO CONTAIN?
paramsofinterest = {'mc','mainc','lurec','uc','respc'};
% GET ALL POSSIBLE PARAM COMBOS
modelstosim = getAllParamCombos(paramsofinterest); 
modelstosim = {'epsilon_init_alpha_mainc'}; modelstofit = modelstosim;

subjnums = unique(toanalyze.subj);nsubjs = length(subjnums); 

for m = 1:length(modelstosim)
    modeltosim = coc_createModels(modelstosim{m}); modeltofit = modeltosim;
    diff = length(list_existing{end})-length([modelstosim{m} '.mat']);
    if sum(contains(list_existing,[modelstosim{m} '.mat' repmat(' ',1,diff)]))==0 %have you run gen-rec on this model already?
        %if no, run it and save results
        realparamlist = []; 
        for subj = 1:nsubjs
            subjnum = subjnums(subj); %grab real subj number
            params = rand(1,modeltosim.nparams);
            onesubj = toanalyze(toanalyze.subj==subjnum,:);
            data{subj} = simulate_cost_model(modeltosim,params,onesubj);
            realparamlist(subj,1:length(params)) = params;
        end
        fnames{m} = [modelstosim{m} '.mat'];
        diff = length(list_existing{end})-length(['fit_' modelstofit{m} '.m']);
        if sum(contains(list_existing,['fit_' modelstofit{m} '.m' repmat(' ',1,diff)]))==0 %no function for running this model, yet
            eval(['copyfile dictate_model.m  model-functions/fit_' modelstofit{m} '.m'])
        end
        eval(['func = @fit_' modelstofit{m} ';']); funcs{m} = func;
        priors{m} = struct('mean',zeros(modeltofit.nparams,1),'variance',6.25); 
        cbm_lap(data, func, priors{m}, fnames{m});
        fname = load(fnames{m});
        cbm   = fname.cbm;
        % look at fitted parameters
        fitparams = applyTrans_parameters(modeltofit,cbm.output.parameters);
        save(['model-details/' modelstofit{m}],'fitparams','realparamlist')
    else
        load(['model-details/' modelstofit{m}])
    end
    figure(1)
    for p = 1:modeltofit.nparams
        subplot(4,2,p)
        scatter(realparamlist(:,p),fitparams(:,p),[],rand(nsubjs,3),'Filled');
        hold on
        plot([0 0],[1 1],'--')
        xlabel(['Real ' modeltofit.paramnames{p}])
        ylabel('Fit values')
        xlim([0 1]); ylim([0 1]);
    end
    fig = gcf; fig.Color = 'w';
    [r,p] = corr(realparamlist,fitparams); %have to do some selecting since there are some nans in the simulated parameter values
    rs = diag(r)
    ps = diag(p)
    reliable = input('does this model fit look reliable? y/n','s'); close 1
    save(['model-details/' modelstosim{m}],'realparamlist','fitparams','reliable')
end

%% Now, add in some hierarchical model assignments, re-simulate, test that
% portion of the HBI
clear data params fitparams onesubj
for m = 1:length(modelstofit)
    file = load(['model-details/' modelstofit{m}]);
    recoverability(m) = file.reliable=='y';
end
true_models = find(~recoverability); %which models weren't recovered well?
ns = floor(nsubjs/length(true_models)); % is their recovery improved by 
true_models = repmat(true_models,1,ns); true_models = [true_models repmat(true_models(1),1,nsubjs-length(true_models))];

% realparamlist = nan(nsubjs,8); 
% for subj = 1:nsubjs
%     subjnum = subjnums(subj); %grab real subj number
%     modeltosim = coc_createModels(modelstosim{true_models(subj)});
%     params = rand(1,modeltosim.nparams);
%     onesubj = toanalyze(toanalyze.subj==subjnum,:);
%     data{subj} = simulate_cost_model(modeltosim,params,onesubj);
%     realparamlist(subj,1:length(params)) = params;
% end

v = 6.25;
for m = 2:length(unique(true_models))
    model_name = modelstofit{unique(true_models(m))};
    modeltofit = coc_createModels(model_name);
    fnames{1} = [model_name '.mat']; 
    funcs{1} = @dictate_model; %I'm not sure how to do this in a generalizable way
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
        ylim([0 1])
    end
    fig = gcf; fig.Color = 'w';
    MSEs = mean((realparamlist(1:size(fitparams,2))-fitparams).^2);
    disp(['MSE = ' num2str(MSEs)])
    
    subplot(4,3,p+1)
    bar(freqs)
    hold on
    %scatter(1:length(modelstofit),real_proportion,'ko','Filled')
    legend({'Fit freqs (HBI)','Real freqs'})
    xlabel('Model')
    xtickangle(45)
    
    [r,p] = corr(fitparams,realparamlist)
    disp(['Param fits for ' model_name])  
    
    reliable = input('does this model fit look reliable? y/n','s'); close 1
    save(['model-details/' model_name],'realparamlist','fitparams','reliable')
end
%% % % FIT REAL SUBJECTS! % % %%
%All of the above stuff can remain as-is. Now I'm just going to format the
%subjects' data properly and hope for the best.
clear all
realsubjectsflag = true; HBI_flag = true;
fitflag = false;
%add relevant paths
addpath('cbm-master/codes/');
addpath('./..'); %other COC code, like model loading
load('./../simdata/toanalyze.mat');
% grab task characteristics from real subjects

paramsofinterest = {'mc','mainc','lurec','uc','respc'};
modelstofit = getAllParamCombos(paramsofinterest);
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];
for m = 1:length(modelstofit) %grab only recoverable models
    file = load(['model-details/' modelstofit{m}]);
    recoverability(m) = file.reliable=='y';
end
modelstofit = modelstofit(recoverability);
function_folder = dir('model-functions'); list_existing = cellstr(string(char(function_folder.name)));
for m = 1:length(modelstofit) %now, cycling through models which are recoverable, create function
    % calls for running those models if they don't already exist
    diff = length(list_existing{end})-length(['fit_' modelstofit{m} '.m']);
    if sum(contains(list_existing,['fit_' modelstofit{m} '.m' repmat(' ',1,diff)]))==0 %no function for running this model, yet
        copyfile 'dictate_model.m' ['model-functions/fit_' modelstofit{m} '.m']
    end
    %otherwise, do nothing
end

nsubjs = length(unique(toanalyze.subj));
for subj = 1:length(unique(toanalyze.subj))
    data{subj} = toanalyze(toanalyze.subj==subj,:);
end

v = 6.25;
for m = 1:length(modelstofit)
    model_name = modelstofit{m}; modeltofit = coc_createModels(model_name);
    fnames{m} = [modelstofit{m} '.mat'];
    eval(['func = @fit_' modelstofit{m} ';']); funcs{m} = func;
    priors{m} = struct('mean',zeros(modeltofit.nparams,1),'variance',v); 
    model_labels{m} = strrep(modelstofit{m},'_','-');
    if fitflag; cbm_lap(data, func, priors{m}, fnames{m}); end
end
% _feb23 and _feb23_b have lurec, mc, and respc fit (consistently)
%fname_hbi = 'HBI_smaller_section_of_search_b.mat'; % run 3 model test to look in more detail
fname_hbi = 'HBI_coc_full_model_search_allsubj.mat'; %big model search over tenable
%models, like the 6 param lurec_mc_respc

%data {cell per subj}, model-specific fitting functions, filenames from
%cbm_lap, %filename for saving full running to
if fitflag
    cbm_hbi(data,funcs,fnames,fname_hbi);
end
fits = load(fname_hbi);
cbm   = fits.cbm;
freqs = cbm.output.model_frequency;

% Use CBM toolbox to print parameter means etc.
% 1st input is the file-address of the file saved by cbm_hbi
% 2nd input: a cell input containing model names
% 3rd input: another cell input containing parameter names of the winning model
[~,best] = max(cbm.output.exceedance_prob);
best_model = coc_createModels(modelstofit{best}); morecomplexmodel = coc_createModels(modelstofit{1});
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
    end
end
%best_model.paramnames = param_names; 
best_model.lowparams = cbm.output.parameters{best};
best_model.highparams = cbm.output.group_mean{best}; best_model.overallfit = cbm.output;
best_model.overallfit.fitmodels = model_labels; best_model.name = model_labels{best};
modelStruct.best_model = best_model; save('HBI_modelStruct.mat','modelStruct');
cbm_hbi_plot(fname_hbi, model_labels, param_names(1:best_model.nparams), transform(1:best_model.nparams))
% this function creates a model comparison plot (exceednace probability and model frequency) as well as 
% a plot of transformed parameters of the most frequent model.

% % run a t-test on miss costs versus 0
% % 2nd input: the index of the model of interest in the cbm file
% k = 1; % model including miss costs
% % 3rd input: the test will be done compared with this value (i.e. this value indicates the null hypothesis)
% m = 0; % here the miss costs parameter should be tested against m=0
% % 4th input: the index of the parameter of interest 
% i = 6; % here the weight parameter is the 5th parameter of the hybrid model
% [p,stats] = cbm_hbi_ttest(cbm,k,m,i);

figure
subplot(2,2,1)
bar(freqs)
title('Fit model freq (HBI)')
xticklabels(model_labels)
xtickangle(45)
fig = gcf; fig.Color = 'w';
best_models = find(cbm.output.model_frequency>0.001);
for m = 1:length(best_models)
    modeltofit = coc_createModels(modelstofit{best_models(m)});
    means = applyTrans_parameters(modeltofit,cbm.output.group_mean{best_models(m)});
    subplot(3,2,1+m)
    for p = 1:modeltofit.nparams
        bar(p,means(p),'FaceColor',paramcolors(p,:))
        hold on
    end
    ylim([-2 2])
    xticks(1:length(means))
    errorbar(means,cbm.output.group_hierarchical_errorbar{best_models(m)},'*k')
    xticklabels(modeltofit.paramnames)
    title('Parameter means from subjects best fit by model')
end

subplot(2,2,4); 
costs = find(costs(1:best_model.nparams)); 
means = applyTrans_parameters(modeltofit,cbm.output.group_mean{2});
bar(means(:,costs),'FaceColor',[0 0.7 0])
hold on
errorbar(means(:,costs),cbm.output.group_hierarchical_errorbar{m}(costs),'*k')
%plot(1:length(costs),best_model.lowparams(:,costs),'--k')
fig = gcf; fig.Color = 'w';
xticklabels(param_names(costs))
xtickangle(30)
ylabel('Mean parameter value')
xlabel('Cost parameter')

%% model simulations and validation figures

model_validation_HBI()

