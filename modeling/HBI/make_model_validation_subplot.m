%% Plot simulated data from HBI model fits
% Using all models, not just the "winning" mainc model (model freq 70%)
% Takes a subplot label as the only argument

% % Load up real subject's data
load('simdata/toanalyze.mat')
%scale things appropriately
toanalyze.BDM = (toanalyze.BDM./25) + 1;

taskcolors = [0.75 0.75 0.75;0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75];
tasklabels = {'1-detect','1-back','3-detect','2-back'};

responsibility = cbm.output.responsibility; %modelstofit = best_model.overallfit.fitmodels;
lowparams = cbm.output.parameters; %for accessibility, grab important info from cbi structures
NFCsubjs = [21 36 22 74 79 73];
% I grabbed two random subjects from each NFC group
% Plot their individual fits by the model

[~,best] = max(cbm.output.exceedance_prob);

% %simulate data from fit parameters, best fitting model for each subject
simdata = []; simdata_onemodel = simdata;
for subj = 1:n
    onesubj = toanalyze(toanalyze.subj == subj,:);
    [score,num] = max(responsibility(subj,:)); %identify model which best fit this subject, in particular
    %then simulate their data with that model
    subj_model = coc_createModels(modelstofit{num});
    % input to simulate_cost_model the model specs, then the transformed
    % parameters for that subject, then that subject's data
    simdata = [simdata; simulate_cost_model(subj_model,applyTrans_parameters(subj_model,lowparams{num}(subj,:)),onesubj)];
end

simdata(:,3) = (simdata(:,3)./25)+1;
ntasks = length(unique(simdata(:,2)));

sim_task1 = []; sim_task2 = []; sim_task3 = [];
onemodel_sim_task1 = []; onemodel_sim_task2 = []; onemodel_sim_task3 = [];
rtask1 = []; rtask2 = []; rtask3 = [];
for task = 1:(ntasks-1)
    for subj = 1:n %Cycle over subjects
        onesubj = simdata(simdata(:,1)==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onesubj(:,4)==task+1)) = onesubj(onesubj(:,4)==task+1,3)';
        % Pick out task iterations and their corresponding fair wage
        % ratings
        eval(['sim_task' num2str(task) ' = [sim_task' num2str(task) '; curve];'])
        
        onereal = toanalyze(toanalyze.subj==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onereal.display==task+1)) = onereal.BDM(onereal.display==task+1,:)';
        eval(['rtask' num2str(task) ' = [rtask' num2str(task) '; curve];'])
    end
end

errorbar(nanmean(rtask1),nanstd(rtask1)/sqrt(n),'Color',taskcolors(2,:),'LineWidth',2)
hold on
errorbar(nanmean(rtask3),nanstd(rtask3)/sqrt(n),'Color',taskcolors(3,:),'LineWidth',2)
errorbar(nanmean(rtask2),nanstd(rtask2)/sqrt(n),'Color',taskcolors(4,:),'LineWidth',2)
errorbar(nanmean(sim_task1),nanstd(sim_task1)/sqrt(n),'--','Color',taskcolors(2,:)+0.1,'LineWidth',2)
errorbar(nanmean(sim_task3),nanstd(sim_task3)/sqrt(n),'--','Color',taskcolors(3,:)+0.1,'LineWidth',2)
errorbar(nanmean(sim_task2),nanstd(sim_task2)/sqrt(n),'--','Color',taskcolors(4,:)+0.1,'LineWidth',2)
legend({'1-back','3-detect','2-back'})
ax = gca; ax.FontSize = 14;
ylim([1.9 3.3])
xlabel('Iteration')
ylabel('Fair wage')
