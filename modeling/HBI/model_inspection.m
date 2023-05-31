function [outputArg1,outputArg2] = model_inspection(m,cbm)
%model_inspection Quickly draw up model information for specific model, m
% First input, m, is the number of the model in the model fitting regime
% Second input, cbm, is the output of the model fitting procedure run in
% HBI_coc.m, and relying on HBI model fitting package written by Piray &
% Daw.
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0; 0 0.5 0.5; 0 1 1];
% Plotting variables

models = cbm.input.model_names;
params = applyTrans_parameters(coc_createModels(models{m}),cbm.output.parameters{m});
nparams = size(params,2);
model_name = strrep(models{m},'-','_');
modeltofit = coc_createModels(models{m});
n = size(params,1);
subset = 1:n;
% default subject subset, can be overwritten by loading files below

load(['model-details/' model_name '.mat']);
load('../../simdata/toanalyze.mat')
% Load generate/recover for this model, too, to see goodness-of-fit


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
    subplot(5,2,p)
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


% Clean up parameter names for plotting purposes
alpha_indx = contains(modeltofit.paramnames,'alpha');
modeltofit.paramnames{alpha_indx} = '\alpha';

epsilon_indx = contains(modeltofit.paramnames,'epsilon');
modeltofit.paramnames{epsilon_indx} = '\epsilon';

costs_indx = contains(modeltofit.paramnames,'c');
names = modeltofit.paramnames(costs_indx);
for i = 1:length(names)
    % format cost parameter names to be more intelligible
    name = names{i};
    temp = strsplit(name,'c');
    temp = temp{1};
    temp = ['c_{' temp '}'];
    names{i} = temp;
end
modeltofit.paramnames(costs_indx) = names;

init_indx = contains(modeltofit.paramnames,'init');
modeltofit.paramnames(init_indx) = {'init_{1-back}','init_{2-back}','init_{3-detect}'};

% Plot real versus fit parameters
figure
rows = ceil(nparams/3); 
columns = ceil(nparams/2);
for p = 1:nparams
    subplot(rows,columns,p)
    scatter(realparamlist(subset,p),fitparams(:,p),[],rand(size(fitparams,1),3),'Filled');
    hold on
    line_coords = min([max(realparamlist(subset,p)) max(fitparams(:,p))]);
    plot([0 line_coords],[0 line_coords],'k--','LineWidth',2)
    xlabel(['Real ' modeltofit.paramnames{p}])
    ylabel(['Fit ' modeltofit.paramnames{p}])
    title(['Model #' num2str(m)])
    %xlim([0 1]); ylim([0 1]);
    ax = gca; ax.FontSize = 14;
end
fig = gcf; fig.Color = 'w';
[r,p] = corr(realparamlist(subset,:),fitparams); %have to do some selecting since there are some nans in the simulated parameter values
rs = diag(r)
ps = diag(p)
disp(['Model ' num2str(m)])
disp([num2str(length(subset)) ' simulated subjects plotted.'])


end

