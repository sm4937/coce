%% Calculate AIC and BIC scores for each model of interest
% In response to PLoS CB reviewer 1, including now a supplementary analysis
% of regular old AIC and BIC scores relative to model frequency scores. Do
% they produce the same results?

clear all

addpath('../../../modeling/')

% LOAD MODELING RESULTS FROM HBI FOLDER
load('../HBI_modelStruct_2023.mat');
cbm = best_model.cbm;

% The results from the full hierarchical procedure
freqs = cbm.output.model_frequency;
[~,assignments] = max(cbm.output.responsibility,[],2);
models_at_play = unique(assignments);

n_models = length(best_model.overallfit.fitmodels);

AIC = NaN(length(assignments),length(freqs));
BIC = AIC;

% Now, load individual model files up here
for mii = 1:n_models
    temp = best_model.overallfit.fitmodels{mii};
    model_name = strrep(temp,'-','_');
    
    % load low-level fits which I specifically saved in this directory
    load([model_name '.mat']);
    
    % get llh for each subject
    llh(:,mii) = cbm.math.loglik;
    
    %log evidence also
    log_evidence(:,mii) = cbm.output.log_evidence;
    
    model_structure = coc_createModels(model_name);
    nparams = model_structure.nparams;
    AIC(:,mii) = -2.*llh(:,mii) + 2.*nparams;
    
    % n trials (ratings per subject) = 32
    BIC(:,mii) = -2*llh(:,mii) + nparams.*log(32);
    
    % what about integrated BIC score?
    % Peter: this is reported, for instance, in Quentin Huys' PLoSCB paper
    % on the go/nogo task. Briefly, if you have a population model, like
    % HBI, then the likelihood of the data integrates out
    % (i.e., marginalizes over) the parameters for any particular subject 
    % according to the population prior. Thus, when we're computing
    % the likelihood, we usually sample from the population prior, 
    % and calculate the average likelihood (not the average LOG likelihood). 
    % The iBIC score uses this integrated likelihood, coupled with the 
    % BIC cost for each population parameter. This way, we don't 
    % BIC for the individual subject-parameters, 
    % since we don't optimize for those at the last stage. 
    
    

    
end

[~,best_AIC] = min(AIC,[],2);
[~,best_BIC] = min(BIC,[],2);
[~,best_LE] = max(log_evidence,[],2);

delta_models = find(contains(best_model.overallfit.fitmodels,'delta'));

% plot results to compare
figure()
subplot(4,1,1)
histogram(assignments,'BinWidth',1)
xlim([1 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
ylabel('# subjects')
title('Best model according to model responsibility')

subplot(4,1,2)
histogram(best_AIC,'BinWidth',1)
xlim([1 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
ylabel('# subjects')
title('According to AIC score')

subplot(4,1,3)
histogram(best_BIC,'BinWidth',1)
xlim([1 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
title('According to BIC score')
ylabel('# subjects')
set(gcf,'color','w')

subplot(4,1,4)
histogram(best_LE,'BinWidth',1)
xlim([1 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
title('According to log evidence')
ylabel('# subjects')
set(gcf,'color','w')
xlabel('Model')

%% Okay, what about model frequency, summed log evidences, and mean AIC/BIC score?

[~,best_by_freq] = max(freqs);
[min_AIC,best_by_AIC] = min(nanmean(AIC));
[min_BIC,best_by_BIC] = min(nanmean(BIC));
[max_LE,best_by_LE] = max(sum(log_evidence));

% plot results to compare
figure()
subplot(4,1,1)
bar(freqs)
hold on
scatter(best_by_freq,0.9,'*k','LineWidth',1.25)
xlim([0 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
ylabel('Model frequency')
title('Best model according to model frequencies')

subplot(4,1,2)
bar(nanmean(AIC)-min_AIC)
hold on
scatter(best_by_AIC,50,'*k','LineWidth',1.25)
xlim([0 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
ylabel('mean AIC - minimum AIC')
title('According to mean AIC score')

subplot(4,1,3)
bar(nanmean(BIC)-min_BIC)
hold on
scatter(best_by_BIC,50,'*k','LineWidth',1.25)
xlim([0 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
ylabel('mean BIC - minimum BIC')
title('According to mean BIC score')

subplot(4,1,4)
bar(sum(log_evidence))
hold on
scatter(best_by_LE,-15000,'*k','LineWidth',1.25)
xlim([0 n_models])
set(gca,'xtick',[],'FontSize',17)
xticks([1 delta_models(1)])
xticklabels({'\alpha','\delta'})
title('According to summed log evidences')
ylabel('sum(log evidence)')
set(gcf,'color','w')
xlabel('Model #')


