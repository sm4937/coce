% TO DO:
% The summing is probably meant to be done to get the marginal
% distributions from the joint distribution over all 3 cost params of
% interest. Not before.


% Make a list of all param names, in the order the modelling procedure sees
% them.
% The indices here should align with the parameters, group_means &
% group_hierarchical_errorbars outputs. Basically all outputs in the cbm
% struct.
all_params_all_models = {'uc','epsilon','initi_1','initi_2','initi_3','missc', ...
    'mainc','matchc','noisec','respc','lurec','errorc','fac','alpha', ...
    'delta_uc','delta_missc','delta_mainc','delta_matchc','delta_noisec', ...
    'delta_respc','delta_lurec','delta_errorc','delta_fac'};
nparams = length(all_params_all_models);

% all possible values for all possible parameters
xs = [-10:0.05:10];

% Write Theta for all the possible parameters across models with empirical
% prior P(Theta)
% Approximate P(Theta)
p_theta = ones(nparams,length(xs));
% Uniform prior for now, bounded at extreme values
    
% take subject s, with posterior model weightings rho^s_j for model j
% [from Payam]

% Get \rho^s_j
rho = cbm.output.responsibility;

% say model j involves parameters Theta_j but not Theta_{j'}

parameters = false(length(all_params_all_models),length(modelstofit));
for j = 1:length(modelstofit)
    model_j = modelstofit{j};
    parameters_j = strsplit(model_j,'_');
    for p = 1:length(parameters_j)
        parameters(contains(all_params_all_models,parameters_j(p)),j) = true;
    end
    if sum(contains(parameters_j,'alpha'))
        parameters(contains(all_params_all_models,'delta'),j) = false;
    end
end

parameters_j_prime = ~parameters;
% and subject s has posterior distribution P(Theta^s_j|D^s) for model j
% where D^s is the data for subject s [from Payam]
% The means are the fit values for each subject
% The stds of these dists are sqrt(cbm.math.qhquad.Ainvdiag{j});
% Saved below in:
% dists_j(:,:,p) = normpdf(xs,fit_values(:,p),subj_stds(p,:)');

% Cycle over each model
marginal = cell(length(all_params_all_models),1);
for j = 1:length(modelstofit)
    fit_values = cbm.output.parameters{j};
    subj_stds = cbm.math.qhquad.Ainvdiag{j};

    % get P(Theta^s_j|D^s)
    dists_j = [];
    for p = 1:sum(parameters(:,j))
        % Individual posterior probability distributions over parameter
        % values
        % As in, each subject has their own posterior and I'm saving it in
        % dists_j
        % subjmeans are in cbm.math.qhquad.theta{j}(p)
        dists_j(:,:,p) = normpdf(xs,fit_values(:,p),subj_stds(p,:)');
        
        p_param = p_theta(parameters(:,j),:);
        p_theta_j = p_param(p,:);
        % This is the "empirical prior", the flat prior
        
        overall_index = find(parameters(:,j));
        overall_index = overall_index(p);
        % In the model space, where does this parameter value live?
        
        pop_prior_j = normpdf(xs,cbm.math.thetabar{j}(p),cbm.math.Sdiag{j}(p));
        % confirm with CBM paper soon that this is really the model's
        % group-level posterior distribution
        
        marginal{overall_index} = [marginal{overall_index}; rho(:,j) .* dists_j(:,:,p) .* p_theta_j];
        % so concatenate, from all models, the distribution of parameter
        % values for each individual parameter
    end
    % approximating the distribution of these parameters this way for now
    % so posterior distribution over all parameter values for all subjects
    
    % then, we might be interested in
    % P(Theta | D^s) =
    %   \sum_j P(Theta,j|D^s)
    % = \sum_j [\rho^s_j P(Theta^s_j|D^s)*P(Theta_{j'})]
        
end

% In this summing step, you can select only subjects in each NFC group
% must make an index which reflects the length of marginal
% So in this case it would be like:
% split = tertileSplit(data.NFC);
split = randi(3,100,1);
for p = 1:length(all_params_all_models)
    highidx = repmat(split==3,size(marginal{p},1)/nsubjs,1);
    temp_high = sum(marginal{p}(highidx,:),1);
    P_theta_given_D_high{p} = temp./sum(temp);
    % Normalize here, which I'm not sure I should be doing
end

% Okay, so of particular interest are mainc, lurec, fac
top_costs = {'mainc','lurec','fac'};
overall_index = contains(all_params_all_models,top_costs);
overall_index(contains(all_params_all_models,'delta')) = false;
top_costs_indx = find(overall_index);

% Now incorporate:
% using the population prior P(Theta_{j'}) for those parameters not in
% model j and the posterior

% high NFC only for now, didn't get to the other groups yet
mainc_dist = sum(P_theta_given_D_high{top_costs_indx(1)},1) .* ones(length(xs),1,1);
lurec_dist = sum(P_theta_given_D_high{top_costs_indx(2)},1) .* ones(1,length(xs),1);
fac_dist = sum(P_theta_given_D_high{top_costs_indx(3)},1) .* ones(1,1,length(xs));

joint_dist = mainc_dist .* lurec_dist .* fac_dist;
joint_dist = joint_dist./(sum(joint_dist(:)));

% Now I have a joint distribution over 

% so then, if you have subjects s_1...s_k in your low tertile group (say),
% you can work out the approximate posterior distribution
% 
% P(Theta ; low_tertile) =
%    \prod_{i=1}^k \sum_j [\rho^{s_i}_j P(Theta^{s_i}_j|D^{s_i})*P(Theta_{j'})]
% 
% in the low_tertile [in practice you should calculate it carefully using
% logs - in the three dimensions of maintenance, lure and false-alarm cast]
% 
% and you could then compare these posteriors for the three groups.
% 
% In practice, if a subject isn't well fit by a model (rho^s_j is low);
% then that subject won't pull parameters that are only in that model away
% from their population prior very much - which is just the property you
% want. 

%% Plot distributions of mainc, lurec, and fac, to compare their distributions to one another.

modelcolors = [1 0 0; 1 0.5 0; 1 0 0.5; 0 0 1; 0 0.5 1; 0 0.7 0; 1 1 0];

params_for_comparison = find(contains(all_params_all_models,{'mainc','lurec','fac'}));
params_for_comparison(end-2:end) = []; % delete "delta" versions of these params

params_from_HBI(1) = find(contains(strsplit(modelstofit{best_models(1)},'_'),'mainc'));
params_from_HBI(2) = find(contains(strsplit(modelstofit{best_models(2)},'_'),'lurec'));
params_from_HBI(3) = find(contains(strsplit(modelstofit{best_models(3)},'_'),'fac'));

figure; hold on;
plot(xs',P_theta_given_D{params_for_comparison(1)},'Color',modelcolors(1,:),'LineWidth',1.5,'DisplayName','Maintenance')
plot(xs',P_theta_given_D{params_for_comparison(2)},'Color',modelcolors(2,:),'LineWidth',1.5,'DisplayName','Interference')
plot(xs',P_theta_given_D{params_for_comparison(3)},'Color',modelcolors(3,:),'LineWidth',1.5,'DisplayName','False Alarm')
legend('Location','Best')
fig = gcf; ax = gca; 
fig.Color = 'w'; ax.FontSize = 14;

for p = 1:length(params_for_comparison)
    group_level_means_sarah(p) = sum(xs.*P_theta_given_D{params_for_comparison(p)});
    HBI_output = cbm.output.group_mean{best_models(p)};
    group_level_means_payam(p) = HBI_output(params_from_HBI(p));
end

disp('Your means are: ')
disp(group_level_means_sarah)
disp('The means from the HBI package are: ')
disp(group_level_means_payam)

% These are different but they don't account for all models containing
% those parameter values. So that needs to be amended, I think, in the
% future. 