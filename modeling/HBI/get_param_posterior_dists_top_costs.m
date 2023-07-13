function [big_posterior,joint,xs] =  get_param_posterior_dists_top_costs(best_model)

%% Obtain posterior distributions over parameter values, taking into account
% different solutions for different models, and different fits for
% different subjects

modelstofit = best_model.overallfit.fitmodels;
cbm = best_model.cbm;
nsubjs = size(best_model.lowparams,1);

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
xs = [-1.5:0.01:1.5]; %previously 0.01 was step size
% bins too big?

% Write Theta for all the possible parameters across models with empirical
% prior P(Theta)
% Approximate P(Theta)
p_theta = ones(nparams,length(xs));
% Uniform prior for now, bounded at extreme values

% say model j involves parameters Theta_j but not Theta_{j'}

parameters = false(length(all_params_all_models),length(modelstofit));
for j = 1:length(modelstofit)
    model_j = modelstofit{j};
    parameters_j = strsplit(model_j,'-');
    for p = 1:length(parameters_j)
        parameters(contains(all_params_all_models,parameters_j(p)),j) = true;
    end
    if sum(contains(parameters_j,'alpha'))
        parameters(contains(all_params_all_models,'delta'),j) = false;
    end
end


% and subject s has posterior distribution P(Theta^s_j|D^s) for model j
% where D^s is the data for subject s [from Payam]
% The means are the fit values for each subject
% The stds of these dists are sqrt(cbm.math.qhquad.Ainvdiag{j});
% Saved below in:
% dists_j(:,:,p) = normpdf(xs,fit_values(:,p),subj_stds(p,:)');

% Okay, so of particular interest are uc, mainc, lurec, fac?
top_costs = {'uc','mainc','lurec'};
overall_index = contains(all_params_all_models,top_costs);
overall_index(contains(all_params_all_models,'delta')) = false;
top_costs_indx = find(overall_index);


p_theta_j_prime = zeros(length(top_costs),length(xs));

big_posterior = zeros(length(xs),length(xs),length(xs));

% Get \rho^s_j
rho = cbm.output.responsibility;

% cycle over each model of interest (fit models)
for j = 1:length(modelstofit)

    rho_j = sum(rho(:,j))/100;

    % does this model in particular contain each of the costs we care
    % about?
    % get model group-level prior
    model = coc_createModels(modelstofit{j});
    idx_1 = find(contains(model.paramnames,top_costs{1})); 
    contains_1 = ~isempty(idx_1);
    if contains_1
        oneD = normpdf(xs,cbm.output.group_mean{j}(idx_1),cbm.math.Sdiag{j}(idx_1));
        %oneD = normpdf(xs,cbm.output.group_mean{j}(idx_1),cbm.output.group_hierarchical_errorbar{j}(idx_1)); 
        p_theta_j_prime(1,:) = p_theta_j_prime(1,:) + rho_j*oneD;
    end

    % 2nd dimension 
    idx_2 = find(contains(model.paramnames,top_costs{2})); 
    contains_2 = ~isempty(idx_2);
    if contains_2
        twoD = normpdf(xs,cbm.output.group_mean{j}(idx_2),cbm.math.Sdiag{j}(idx_2));
        %twoD = normpdf(xs,cbm.output.group_mean{j}(idx_2),cbm.output.group_hierarchical_errorbar{j}(idx_2));        
        p_theta_j_prime(2,:) = p_theta_j_prime(2,:) + rho_j*twoD;
    end

    % 3rd dimension 
    idx_3 = find(contains(model.paramnames,top_costs{3})); 
    contains_3 = ~isempty(idx_3);
    if contains_3
        threeD = normpdf(xs,cbm.output.group_mean{j}(idx_3),cbm.math.Sdiag{j}(idx_3));
        %threeD = normpdf(xs,cbm.output.group_mean{j}(idx_3),cbm.output.group_hierarchical_errorbar{j}(idx_3));        
        p_theta_j_prime(3,:) = p_theta_j_prime(3,:) + rho_j*threeD;
    end

end

p_theta_j_prime = p_theta_j_prime./sum(p_theta_j_prime,2);
% normalize so these priors sum to 1 
% this provides an aggregrate prior from all models
plotflag = false;
if plotflag
    figure
    plot(p_theta_j_prime(1,:),'DisplayName',top_costs{1}); hold on; plot(p_theta_j_prime(2,:),'DisplayName',top_costs{2}); plot(p_theta_j_prime(3,:),'DisplayName',top_costs{3}); legend('Location','Best')
    title('Group-level priors')
end

%Cycle over each model
joint = cell(nsubjs,1);

% take subject s, with posterior model weightings rho^s_j for model j
% [from Payam]
for s = 1:nsubjs

    joint{s} = zeros(length(xs),length(xs),length(xs));
    % initialize to blank

    for j = 1:length(modelstofit)
        if rho(s,j) ~= 0 % if this model didn't load on this subject at all, skip all
            % this computation
            % not sure if rho is ever truly 0 but we'll see
            fit_values = cbm.output.parameters{j};

            % fit_values = log(fit_values);

            % do I need to log transform these fit_values? are
            % they exp transforms of the numbers used in
            % fitting?
            subj_stds = sqrt(cbm.math.qhquad.Ainvdiag{j});

            contained_costs = false(1,length(top_costs)); 
            contained_costs(1) = parameters(top_costs_indx(1),j);
            contained_costs(2) = parameters(top_costs_indx(2),j);
            contained_costs(3) = parameters(top_costs_indx(3),j);

            % grab population-level priors for parameters not in model j 
            pop_prior_j_prime = ones(length(top_costs),length(xs));
            % make this 1? make this nan?

            dists_j = pop_prior_j_prime;

            for j_prime = 1:length(top_costs)
                if ~contained_costs(j_prime)
                    % select only costs NOT in this model
                    pop_prior_j_prime(j_prime,:) =  pop_prior_j_prime(j_prime,:) .* p_theta_j_prime(j_prime,:);

                else % get P(Theta^s_j|D^s)
                     % Individual posterior probability distributions over parameter
                    % values
                    % As in, each subject has their own posterior and I'm saving it in
                    % dists_j
                    % subjmeans are in cbm.math.qhquad.theta{j}(p)
                    cost_name = top_costs{j_prime};
                    model = coc_createModels(modelstofit{j});
                    idx_cost = find(contains(model.paramnames,cost_name)); 
                    dists_j(j_prime,:) = normpdf(xs,fit_values(s,idx_cost),subj_stds(idx_cost,s)');
                    dists_j(j_prime,:) = dists_j(j_prime,:)./sum(dists_j(j_prime,:));
                end
            end
            % confirm with CBM paper soon that this is really the model's
            % group-level posterior distribution

            % the underlying priors are not the issue BUT
            % the conversion to the 3D arrays is the issue
            plotflag = false;
            if plotflag & rand()<0.5
                figure(1)
                plot(pop_prior_j_prime(1,:)+rand(),'DisplayName','MainC','LineWidth',2); 
                hold on; 
                % jitter to see all 6 lines just to confirm
                plot(pop_prior_j_prime(2,:)+rand(),'DisplayName','LureC','LineWidth',2); 
                plot(pop_prior_j_prime(3,:)+rand(),'DisplayName','FA C','LineWidth',2);  legend('Location','Best')
                plot(dists_j(1,:)+rand(),'DisplayName','Main C','LineWidth',2)
                plot(dists_j(2,:)+rand(),'DisplayName','Lure C','LineWidth',2)
                plot(dists_j(3,:)+rand(),'DisplayName','FA C','LineWidth',2)
                hold off;
            end

            pop_prior_j_prime_3D = pop_prior_j_prime(1,:) .* pop_prior_j_prime(2,:)' .* reshape(pop_prior_j_prime(3,:),[1 1 length(xs)]);
            dists_j_3D = dists_j(1,:) .* dists_j(2,:)' .* reshape(dists_j(3,:),[1 1 length(xs)]);

            plotflag = false;
            if plotflag
                figure(5)
                plot(squeeze(sum(sum(pop_prior_j_prime_3D,1),3)),'DisplayName','Lure C Prior','LineWidth',1.5)
                hold on; plot(squeeze(sum(sum(pop_prior_j_prime_3D,2),3)),'DisplayName','Main C Prior','LineWidth',1.5)
                plot(squeeze(sum(sum(pop_prior_j_prime_3D,1),2)),'DisplayName','False alarm C Prior','LineWidth',1.5)
                plot(squeeze(sum(sum(dists_j_3D,1),3)),'DisplayName','Lure C','LineWidth',1.5)
                plot(squeeze(sum(sum(dists_j_3D,2),3)),'DisplayName','Main C','LineWidth',1.5)
                plot(squeeze(sum(sum(dists_j_3D,1),2)),'DisplayName','False alarm C','LineWidth',1.5)
                legend('Location','Best')
                hold off;
            end


            % We are interested in
            % P(Theta | D^s) =
            %   \sum_j P(Theta,j|D^s)
            % = \sum_j [\rho^s_j P(Theta^s_j|D^s)*P(Theta_{j'})]

            joint{s} = joint{s} + (rho(s,j) .* dists_j_3D .* pop_prior_j_prime_3D);
            % Sum over all models j
            % THIS MIGHT BE WRONG - MOVE RHO TO EXTERNAL PART OF CODE
            % WHERE SUBJECTS GET PUT INTO NFC GROUPS?

        end % of excluding s, j combinations with 0 as rho value

    end % end of model loop

    plotflag = false;
    if plotflag
        figure(5); 
        plot(xs,squeeze(sum(sum(joint{s},1),3)),'DisplayName','Lure C','LineWidth',1.5)
        hold on;
        plot(xs,squeeze(sum(sum(joint{s},2),3)),'DisplayName','Main C','LineWidth',1.5)
        plot(xs,squeeze(sum(sum(joint{s},1),2)),'DisplayName','False alarm C','LineWidth',1.5)
        hold off
        legend('location','best')
        pause(0.1)
        title(['Posteriors at Subject #' num2str(s)])
    end

    non_zero_default = 0.000000000000000000001;
    joint{s}(joint{s}==0) = non_zero_default;
    big_posterior = big_posterior + log(joint{s});

end % end of subject loop

save('joint_cost_distributions_2023.mat','joint','big_posterior','xs','-v7.3')






% xs = exp(xs);
% group_level_means_trans(1) = sum(xs'.*mainc_dist);
% group_level_means_trans(2) = sum(xs.*lurec_dist);
% group_level_means_trans(3) = sum(xs'.*fac_dist);
% disp('The transformed means are: ')
% disp(group_level_means_sarah)

% These are different but they don't account for all models containing
% those parameter values. So that needs to be amended, I think, in the
% future. 