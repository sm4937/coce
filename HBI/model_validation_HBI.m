%% Plot simulated data from HBI model fits
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
for subj = 1:nsubjs
    onesubj = toanalyze(toanalyze.subj == subj,:);
    [score,num] = max(responsibility(subj,:)); %identify model which best fit this subject, in particular
    %then simulate their data with that model
    subj_model = coc_createModels(modelstofit{num});
    % input to simulate_cost_model the model specs, then the transformed
    % parameters for that subject, then that subject's data
    simdata = [simdata; simulate_cost_model(subj_model,applyTrans_parameters(subj_model,lowparams{num}(subj,:)),onesubj)];
    
    one_best_model = coc_createModels(modelstofit{best});
    simdata_onemodel = [simdata_onemodel; simulate_cost_model(one_best_model,applyTrans_parameters(one_best_model,lowparams{best}(subj,:)),onesubj)];
    % simulate data from the model w highest exceedance probability, alone
    % does accounting for other possible models include model validation?
end

% for subj = 1:nsubjs
%     simdata = [simdata; data{subj}];
% end

% % Load up real subject's data
load('../simdata/toanalyze.mat')
%scale things appropriately
simdata(:,3) = (simdata(:,3)./25)+1;
simdata_onemodel(:,3) = (simdata_onemodel(:,3)./25)+1;
toanalyze.BDM = (toanalyze.BDM./25) + 1;
ntasks = length(unique(simdata(:,2)));

sim_task1 = []; sim_task2 = []; sim_task3 = [];
onemodel_sim_task1 = []; onemodel_sim_task2 = []; onemodel_sim_task3 = [];
rtask1 = []; rtask2 = []; rtask3 = [];
for task = 1:(ntasks-1)
    for subj = 1:nsubjs %Cycle over subjects
        onesubj = simdata(simdata(:,1)==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onesubj(:,4)==task+1)) = onesubj(onesubj(:,4)==task+1,3)';
        % Pick out task iterations and their corresponding fair wage
        % ratings
        eval(['sim_task' num2str(task) ' = [sim_task' num2str(task) '; curve];'])
        
        onesubj = simdata_onemodel(simdata_onemodel(:,1)==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onesubj(:,4)==task+1)) = onesubj(onesubj(:,4)==task+1,3)';
        % Pick out task iterations and their corresponding fair wage
        % ratings
        eval(['onemodel_sim_task' num2str(task) ' = [onemodel_sim_task' num2str(task) '; curve];'])
        
        onereal = toanalyze(toanalyze.subj==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onereal.display==task+1)) = onereal.BDM(onereal.display==task+1,:)';
        eval(['rtask' num2str(task) ' = [rtask' num2str(task) '; curve];'])
    end
end

figure
subplot(2,2,1)
% Plot real subject data for comparison
errorbar([nanmean(toanalyze.BDM(toanalyze.display==2,:)) nanmean(toanalyze.BDM(toanalyze.display==4,:)) nanmean(toanalyze.BDM(toanalyze.display==3,:))],[nanstd(toanalyze.BDM(toanalyze.display==2,:)) nanstd(toanalyze.BDM(toanalyze.display==4,:)) nanstd(toanalyze.BDM(toanalyze.display==3,:))]/sqrt(nsubjs), ...
    'k','LineWidth',1.5,'DisplayName','Real');
hold on
errorbar([nanmean(simdata_onemodel(simdata_onemodel(:,4)==2,3)) nanmean(simdata_onemodel(simdata_onemodel(:,4)==4,3)) nanmean(simdata_onemodel(simdata_onemodel(:,4)==3,3))],[nanstd(simdata_onemodel(simdata_onemodel(:,4)==2,3)) nanstd(simdata_onemodel(simdata_onemodel(:,4)==4,3)) nanstd(simdata_onemodel(simdata_onemodel(:,4)==3,3))]./sqrt(nsubjs), ...
    'k--','LineWidth',1.5,'DisplayName','Simulated')
xticks([1:3])
xlim([0.5 3.5])
ylim([1.9 3.3])
xticklabels({'1-back','3-detect','2-back'})
xlabel('Task')
ylabel('Mean fair wage')
legend('Location','Best')
title('Real data vs. data sim with only best model')
fig = gcf; fig.Color = 'w';
ax = gca; ax.FontSize = 14;

subplot(2,2,3)
errorbar(nanmean(rtask1),nanstd(rtask1)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
hold on
errorbar(nanmean(rtask3),nanstd(rtask3)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
errorbar(nanmean(rtask2),nanstd(rtask2)/sqrt(nsubjs),'Color',taskcolors(4,:),'LineWidth',2)
errorbar(nanmean(onemodel_sim_task1),nanstd(onemodel_sim_task1)/sqrt(nsubjs),'--','Color',taskcolors(2,:)+0.1,'LineWidth',2)
errorbar(nanmean(onemodel_sim_task3),nanstd(onemodel_sim_task3)/sqrt(nsubjs),'--','Color',taskcolors(3,:)+0.1,'LineWidth',2)
errorbar(nanmean(onemodel_sim_task2),nanstd(onemodel_sim_task2)/sqrt(nsubjs),'--','Color',taskcolors(4,:)+0.1,'LineWidth',2)
legend({'1-back','3-detect','2-back'})
ax = gca; ax.FontSize = 14;
ylim([1.9 3.3])
xlabel('Iteration')
ylabel('Fair wage')

subplot(2,2,2)
% Plot real subject data for comparison
errorbar([nanmean(toanalyze.BDM(toanalyze.display==2,:)) nanmean(toanalyze.BDM(toanalyze.display==4,:)) nanmean(toanalyze.BDM(toanalyze.display==3,:))],[nanstd(toanalyze.BDM(toanalyze.display==2,:)) nanstd(toanalyze.BDM(toanalyze.display==4,:)) nanstd(toanalyze.BDM(toanalyze.display==3,:))]/sqrt(nsubjs), ...
    'k','LineWidth',1.5,'DisplayName','Real');
hold on
errorbar([nanmean(simdata(simdata(:,4)==2,3)) nanmean(simdata(simdata(:,4)==4,3)) nanmean(simdata(simdata(:,4)==3,3))],[nanstd(simdata(simdata(:,4)==2,3)) nanstd(simdata(simdata(:,4)==4,3)) nanstd(simdata(simdata(:,4)==3,3))]./sqrt(nsubjs), ...
    'k--','LineWidth',1.5,'DisplayName','Simulated')
xticks([1:3])
xlim([0.5 3.5])
ylim([1.9 3.3])
xticklabels({'1-back','3-detect','2-back'})
xlabel('Task')
ylabel('Mean fair wage')
legend('Location','Best')
title('Real data vs. data sim with all models')
fig = gcf; fig.Color = 'w';
ax = gca; ax.FontSize = 14;


subplot(2,2,4)
errorbar(nanmean(rtask1),nanstd(rtask1)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
hold on
errorbar(nanmean(rtask3),nanstd(rtask3)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
errorbar(nanmean(rtask2),nanstd(rtask2)/sqrt(nsubjs),'Color',taskcolors(4,:),'LineWidth',2)
errorbar(nanmean(sim_task1),nanstd(sim_task1)/sqrt(nsubjs),'--','Color',taskcolors(2,:)+0.1,'LineWidth',2)
errorbar(nanmean(sim_task3),nanstd(sim_task3)/sqrt(nsubjs),'--','Color',taskcolors(3,:)+0.1,'LineWidth',2)
errorbar(nanmean(sim_task2),nanstd(sim_task2)/sqrt(nsubjs),'--','Color',taskcolors(4,:)+0.1,'LineWidth',2)
legend({'1-back','3-detect','2-back'})
ax = gca; ax.FontSize = 14;
ylim([1.9 3.3])
xlabel('Iteration')
ylabel('Fair wage')


subjs = NFCsubjs;
tasksymbols = {'s','o','d'}; %differentiate task with scatter symbol
% and trial number with how dark the scatter symbol is
ntrials = sum(simdata(:,1)==1); trialcolors = linspace(0.1,1,ntrials); trialcolors = repmat(trialcolors',1,3); 

figure
taskcount = false(1,3);
corr_mat = [];
for i = 1:nsubjs
    subj = i;
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3); 
    onesim_onemodel = simdata_onemodel(simdata_onemodel(:,1)==subj,:);
    simulated_onemodel = onesim_onemodel(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; 
    for trial = 1:sum(~isnan(real))
        task = onereal.display(trial);
        if ~isnan(task)
            if sum(NFCsubjs==subj)>0
                subplot(3,2,find(NFCsubjs==subj))
                tasksymbol = tasksymbols{task-1};
                if ~taskcount(task-1)
                    scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled','DisplayName',tasklabels{task-1})
                    taskcount(task-1) = true; %have you already labeled this task?
                else
                    scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled')
                end
                hold on
                line
            end
            % now that plot made
            try
                corr_mat = [corr_mat; real(trial) simulated(trial) simulated_onemodel(trial)];
            catch
                % do nothing
                % if there are some simulated trials without ratings bc the
                % task index wasn't saved, then move on
            end
        end
    end
    %legend('location','best')
    %legend({'Real','Simulated'})
    xlabel('real values'); ylabel('simulated');
    ax = gca; ax.FontSize = 14;
end
fig = gcf; fig.Color = 'w';

% Run r-squared analysis to compare best model to best-model-per-subject
% validation

[r_onemodel,p_onemodel] = corr(corr_mat(:,1),corr_mat(:,3));
[r_allmodels,p_allmodels] = corr(corr_mat(:,1),corr_mat(:,2));

disp(['R-squared for one model (mainc): ' num2str(r_onemodel^2)])
disp(['R-squared including all models: ' num2str(r_allmodels^2)])
