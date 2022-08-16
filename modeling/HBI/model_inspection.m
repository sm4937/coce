function [outputArg1,outputArg2] = model_inspection(m,cbm)
%model_inspection Quickly draw up model information for specific model, m
% First input, m, is the number of the model in the model fitting regime
% Second input, cbm, is the output of the model fitting procedure run in
% HBI_coc.m, and relying on HBI model fitting package written by Piray &
% Daw.
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];
% Plotting variables

models = cbm.input.model_names;
params = applyTrans_parameters(coc_createModels(models{m}),cbm.output.parameters{m});
nparams = size(params,2);
model_name = models{m};
modeltofit = coc_createModels(models{m});
n = size(params,1);
subset = 1:n;
% if no subset specified in 

load(['model-details/' model_name '.mat']);
load('../../simdata/toanalyze.mat')
% Load generate/recover for this model, too, to see goodness-of-fit
% Plot real versus fit parameters

figure
rows = ceil(nparams/3); 
columns = ceil(nparams/2);
for p = 1:nparams
    subplot(rows,columns,p)
    scatter(realparamlist(subset,p),fitparams(:,p),[],rand(size(fitparams,1),3),'Filled');
    hold on
    plot([0 0],[1 1],'--')
    xlabel(['Real ' modeltofit.paramnames{p}])
    ylabel('Fit values')
    xlim([0 1]); ylim([0 1]);
    ax = gca; ax.FontSize = 14;
end
fig = gcf; fig.Color = 'w';
[r,p] = corr(realparamlist(subset,:),fitparams); %have to do some selecting since there are some nans in the simulated parameter values
rs = diag(r)
ps = diag(p)
disp(['Model ' num2str(m)])
disp([num2str(length(subset)) ' simulated subjects plotted.'])

% Simulate behavior on this model
simdata = [];
for subj = 1:n
    onesubj = toanalyze(toanalyze.subj == subj,:);
    %then simulate their data with that model
    % input to simulate_cost_model the model specs, then the transformed
    % parameters for that subject, then that subject's data
    simdata = [simdata; simulate_cost_model(modeltofit,params(subj,:),onesubj)];
end

%rank subjects by how they're fit by the model
for subj = 1:n
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; %real = (real./25)+1;
    idx = true(length(simulated),1); idx(isnan(real(idx))) = false; 
    MSE(subj) = sqrt(nanmean((real(idx)-simulated(idx)).^2));
end
subjrankings = sortrows([MSE' (1:n)'],1);

figure
for p = 1:nparams 
    subplot(4,2,p)
    for i = 1:length(subjrankings)
        subj = subjrankings(i,2);
        scatter(i,params(subj,p),[],paramcolors(p,:),'Filled')
        hold on
        name = modeltofit.paramnames{p};
        title(['Fit ' name])
        xlabel('Fit rank')
        ylabel('Param value')
        ax = gca; ax.FontSize = 14;
    end
end
fig = gcf; fig.Color = 'w';


end

