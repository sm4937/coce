%% Plot simulated data from HBI model fits
taskcolors = [0.75 0.75 0.75;0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75];
paramcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];
tasklabels = {'1-detect','1-back','3-detect','2-back'};

responsibility = cbm.output.responsibility; modelstofit = best_model.overallfit.fitmodels;
lowparams = cbm.output.parameters; %for accessibility, grab important info from cbi structures

% %simulate data from fit parameters, best fitting model for each subject
simdata = [];
for subj = 1:nsubjs
    onesubj = toanalyze(toanalyze.subj == subj,:);
    [score,num] = max(responsibility(subj,:)); %identify model which best fit this subject, in particular
    %then simulate their data with that model
    subj_model = coc_createModels(modelstofit{num});
    % input to simulate_cost_model the model specs, then the transformed
    % parameters for that subject, then that subject's data
    simdata = [simdata; simulate_cost_model(subj_model,applyTrans_parameters(subj_model,lowparams{num}(subj,:)),onesubj)];
end
% for subj = 1:nsubjs
%     simdata = [simdata; data{subj}];
% end

% % Load up real subject's data
load('simdata/toanalyze.mat')



figure
subplot(2,2,1)
% Plot real subject data for comparison



%scale things appropriately
simdata(:,3) = (simdata(:,3)./25)+1;
ntasks = length(unique(simdata(:,2)));

subplot(2,2,2)
errorbar([nanmean(simdata(simdata(:,4)==2,3)) nanmean(simdata(simdata(:,4)==4,3)) nanmean(simdata(simdata(:,4)==3,3))],[nanstd(simdata(simdata(:,4)==2,3)) nanstd(simdata(simdata(:,4)==4,3)) nanstd(simdata(simdata(:,4)==3,3))]./sqrt(nsubjs),'k','LineWidth',1.5)
hold on 
% for subj = 1:nsubjs
%     onesubj = simdata(simdata(:,1)==subj,:);
%     errorbar([nanmean(onesubj(onesubj(:,2)==2,3)) nanmean(onesubj(onesubj(:,2)==4,3)) nanmean(onesubj(onesubj(:,2)==3,3))],[nanstd(onesubj(onesubj(:,2)==2,3)) nanstd(onesubj(onesubj(:,2)==4,3)) nanstd(onesubj(onesubj(:,2)==3,3))]/sqrt(length(onesubj)),'LineWidth',1);
% end
title('Mean rating of tasks (sim)')
xticks([1:3])
xlim([0.5 3.5])
ylim([1 5])
xticklabels({'1-back','3-detect','2-back'})
xlabel('Task')
ylabel('Mean rating (sim)')

subplot(2,2,4)
task1 = []; task2 = []; task3 = [];
for task = 1:(ntasks-1)
    for subj = 1:nsubjs
        onesubj = simdata(simdata(:,1)==subj,:);
        curve = NaN(1,11);
        curve(1:sum(onesubj(:,4)==task+1)) = onesubj(onesubj(:,4)==task+1,3)';
        eval(['task' num2str(task) ' = [task' num2str(task) '; curve];'])
    end
end
errorbar(nanmean(task1),nanstd(task1)/sqrt(nsubjs),'Color',taskcolors(1,:),'LineWidth',2)
hold on
errorbar(nanmean(task2),nanstd(task2)/sqrt(nsubjs),'Color',taskcolors(2,:),'LineWidth',2)
errorbar(nanmean(task3),nanstd(task3)/sqrt(nsubjs),'Color',taskcolors(3,:),'LineWidth',2)
legend({'1-back','2-back','3-detect'})
title('Group learning curves for each task (sim)')
fig = gcf; fig.Color = 'w';

subjs = randperm(nsubjs,6);
tasksymbols = {'s','o','d'}; %differentiate task with scatter symbol
% and trial number with how dark the scatter symbol is
ntrials = sum(simdata(:,1)==1); trialcolors = linspace(0.1,1,ntrials); trialcolors = repmat(trialcolors',1,3); 

figure
taskcount = false(1,3);
for i = 1:length(subjs)
    subj = subjs(i);
    subplot(3,2,i)
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3); 
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; real = (real./25)+1;
    for trial = 1:sum(~isnan(real))
        task = onereal.display(trial);
        if ~isnan(task)
            tasksymbol = tasksymbols{task-1};
            if ~taskcount(task-1)
                scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled','DisplayName',tasklabels{task-1})
                taskcount(task-1) = true; %have you already labeled this task?
            else
                scatter(real(trial),simulated(trial),[],trialcolors(trial,:),tasksymbol,'Filled')
            end
            hold on
            line
            %re-formatting for Peter
%             scatter(trial,real(trial),[],taskcolors(task,:),tasksymbol,'Filled')
%             hold on
%             scatter(trial,simulated(trial),[],taskcolors(task,:),'o','Filled')
        end
    end
    title(['Subj ' num2str(subj)])
    %legend('location','best')
    %legend({'Real','Simulated'})
    xlabel('real values'); ylabel('simulated');
end
fig = gcf; fig.Color = 'w';

%rank subjects by how they're fit by the model
for subj = 1:nsubjs
    onesim = simdata(simdata(:,1)==subj,:);
    simulated = onesim(:,3);
    onereal = toanalyze(toanalyze.subj==subj,:);
    real = onereal.BDM; real = (real./25)+1;
    idx = true(length(simulated),1); idx(isnan(real(idx))) = false; 
    MSE(subj) = sqrt(nanmean((real(idx)-simulated(idx)).^2));
end
subjrankings = sortrows([MSE' (1:nsubjs)'],1);

figure
nparams = length(best_model.paramnames);
for p = 1:nparams 
    subplot(4,2,p)
    for i = 1:length(subjrankings)
        subj = subjrankings(i,2);
        scatter(i,best_model.lowparams(subj,p),[],paramcolors(p,:),'Filled')
        hold on
        name = best_model.paramnames{p};
        title(['Fit ' name])
        xlabel('Fit rank')
        ylabel('Param value')
    end
end
fig = gcf; fig.Color = 'w';