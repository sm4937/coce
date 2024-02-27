%% Plot simulated data from HBI model fits
taskcolors = [0.75 0.75 0.75;0 0.75 0.75; 0.75 0.75 0; 0.75 0 0.75];
tasklabels = {'1-detect','1-back','3-detect','2-back'};

responsibility = cbm.output.responsibility; %modelstofit = best_model.overallfit.fitmodels;
lowparams = cbm.output.parameters; %for accessibility, grab important info from cbi structures
% I grabbed two random subjects from each NFC group
% Plot their individual fits by the model

[~,best] = max(cbm.output.protected_exceedance_prob);

% First, determine relevant models
[score,assignments] = max(responsibility,[],2); %identify model which best fit this subject, in particular
costs_of_interest = {'uc','mainc','lurec','fac'};

% % compile list of subjects' best models, sorted by which costs they contain
list_by_cost_parameter = NaN(nsubjs,length(costs_of_interest));
for cii = 1:length(costs_of_interest)
    
    temp = find(contains(cbm.input.model_names,costs_of_interest{cii}))';
    % this will pull subjects out who have that cost (cii) in their
    % best-fitting model
    adherents = find(sum(assignments==temp',2));
    list_by_cost_parameter(1:length(adherents),cii) = adherents;
    
end

% % run r-squared analyses using these models, to compare them
for cii = 1:length(costs_of_interest)
    
    subjects_with_cost = list_by_cost_parameter(:,cii);
    % trim extra entries (stored as NaN);
    subjects_with_cost = subjects_with_cost(~isnan(subjects_with_cost));
    cost_name = costs_of_interest{cii};
    
    for sii = 1:length(subjects_with_cost)
        subj = subjects_with_cost(sii);
        % get data out
        onesubj = toanalyze(toanalyze.subj == subj,:);
        
        % assign subject to the correct model, containing this cost
        num = assignments(subj);
        subj_model = coc_createModels(modelstofit{num});
        eval([cost_name '_models{subj} = subj_model;'])
        eval([cost_name '_models_params{subj} = applyTrans_parameters(subj_model,lowparams{num}(subj,:));'])
    end
    
    % Run r-squared analysis to compare r-squared provided by models
    % containing each cost
    
    nboot = 10;
    eval(['mean_r_squared_' cost_name ' = get_mean_r_squared(nboot, subjects_with_cost, ' cost_name '_models, ' cost_name '_models_params, toanalyze)'])
    
end

%% The previous analysis is strange, because there are just MORE parameters
% in the best models containing the NON-update costs (most of them also include
% update costs!)

% Let's do something a little less intuitive, but which perhaps makes more
% sense (in terms of comparing r-squareds across costs).
model = modelstofit{30};
all_costs_model = coc_createModels(model);
all_costs_model_repeated = repmat({all_costs_model},nsubjs,1);

for sii = 1:nsubjs
    % get param values per subject in correct format
    all_costs_model_params{sii} = applyTrans_parameters(all_costs_model,lowparams{30}(sii,:));
end

nboot = 500;
mean_r_squared_all_costs = get_mean_r_squared(nboot, 1:nsubjs, all_costs_model_repeated, all_costs_model_params, toanalyze);

% eliminate each cost in turn
for cii = 1:length(costs_of_interest)
    cost_name = costs_of_interest{cii};
    to_zero_out = find(contains(all_costs_model.paramnames,cost_name));
    for sii = 1:nsubjs
        params = all_costs_model_params{sii};
        params(:,to_zero_out) = 0;
        eval(['no_' cost_name '_params{sii} = params;'])
    end
    
    % Run r-squared analysis to compare r-squared when you remove 1 cost at a
    % time (by zeroing out the associated cost parameter)
    eval(['mean_r_squared_no_' cost_name ' = get_mean_r_squared(nboot, 1:nsubjs, all_costs_model_repeated, no_' cost_name '_params, toanalyze)'])
    
end

