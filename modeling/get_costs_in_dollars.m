function [costs,costs_per_unit] = get_costs_in_dollars(cost_parameters,cost_labels,toanalyze)
%get_costs_in_dollars Calculates cost (in dollars) of each costly task
%component (costs of cognition & of errors)
%   cost_parameters are the magnitudes of the cost parameters, as fit by
%   the models
%   costs are the transformed values

cost_scaling = 50; %value in set_param_values for 1 parameter models
cost_scaling = cost_scaling./25 - 1; % scale to match BDM ratings

costs = cost_parameters.*cost_scaling;

col_names = toanalyze.Properties.VariableNames;
relevant = contains(col_names,cost_labels,'IgnoreCase',true);
relevant(contains(col_names,'otherlures','IgnoreCase',true)) = false;
columns = table2array(toanalyze(:,relevant));

%since z-scoring happened within subject, this transform should, too

nsubjs = length(unique(toanalyze.subj));

for i = 1:nsubjs
    % cycle over subjects, extract distribution of cost components for each
    one_subj_columns = columns(toanalyze.task>1 & toanalyze.subj==i,:);
    % constrain to tasks with ratings
    
    %transform back from z-scores to obtain mean steps from minimum,
    %thereby mean cost of updating, interference, etc., normalized for
    %comparison across them
    floor_over_rounds = min(one_subj_columns);
    std_over_rounds = nanstd(one_subj_columns);
    steps_above_minimum = (one_subj_columns - floor_over_rounds)./std_over_rounds;
    mean_steps(i,:) = nanmean(steps_above_minimum);
    costs_per_unit(i,:) = costs./std_over_rounds;

end

costs = costs.*mean_steps;



end

