function [big_posterior,joint,xs] =  get_param_posterior_dists_all_costs(best_model)

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
xs = [-1.5:0.05:1.5]; %previously 0.01 was step size
% now we're in 6D...
% bins too big?
% doesn't work otherwise though

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

%top_costs = {'uc','missc','mainc','respc','lurec','fac'};
top_costs = {'uc','mainc','lurec','fac'};
overall_index = contains(all_params_all_models,top_costs);
overall_index(contains(all_params_all_models,'delta')) = false;
top_costs_indx = find(overall_index);


p_theta_j_prime = zeros(length(top_costs),length(xs));

size_posterior = repmat(length(xs),1,length(top_costs));
big_posterior = zeros(size_posterior);

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
    
    % 4th  
    idx_4 = find(contains(model.paramnames,top_costs{4})); 
    contains_4 = ~isempty(idx_4);
    if contains_4
        fourD = normpdf(xs,cbm.output.group_mean{j}(idx_4),cbm.math.Sdiag{j}(idx_4));
        %threeD = normpdf(xs,cbm.output.group_mean{j}(idx_3),cbm.output.group_hierarchical_errorbar{j}(idx_3));        
        p_theta_j_prime(4,:) = p_theta_j_prime(4,:) + rho_j*fourD;
    end
    
end

p_theta_j_prime = p_theta_j_prime./sum(p_theta_j_prime,2);
% normalize so these priors sum to 1 
% this provides an aggregrate prior from all models

plotflag = false;
if plotflag
    figure
    for cii = 1:length(top_costs)
        plot(p_theta_j_prime(cii,:),'DisplayName',top_costs{cii}); hold on;
    end
    legend('Location','Best')
    title('Group-level priors')
end

%Cycle over each model
joint = cell(nsubjs,1);

% take subject s, with posterior model weightings rho^s_j for model j
% [from Payam]
for s = 1:nsubjs

    joint{s} = zeros(size_posterior);
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
            for cii = 1:length(top_costs)
                contained_costs(cii) = parameters(top_costs_indx(cii),j);
            end

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
            
%             pop_prior_j_prime_6D = pop_prior_j_prime(1,:) .* pop_prior_j_prime(2,:)' .* reshape(pop_prior_j_prime(3,:),[1 1 length(xs)]) .* reshape(pop_prior_j_prime(4,:),[1 1 1 length(xs)]);
%             dists_j_6D = dists_j(1,:) .* dists_j(2,:)' .* reshape(dists_j(3,:),[1 1 length(xs)]) .* reshape(dists_j(4,:),[1 1 1 length(xs)]);
            
            pop_prior_j_prime_6D = pop_prior_j_prime(1,:) .* pop_prior_j_prime(2,:)';
            dists_j_6D = dists_j(1,:) .* dists_j(2,:)';

            for cii = 3:length(top_costs)
                shape_vec = ones(1,length(top_costs));
                shape_vec(cii) = length(xs);
                pop_prior_j_prime_6D = pop_prior_j_prime_6D .* reshape(pop_prior_j_prime(cii,:),shape_vec);
                dists_j_6D = dists_j_6D .* reshape(dists_j(cii,:),shape_vec);
            end

            plotflag = false;
            if plotflag
                figure(7)
                plot(squeeze(sum(sum(sum(pop_prior_j_prime_6D,1),3),4)),'DisplayName','Main C Prior','LineWidth',1.5)
                hold on;
                plot(squeeze(sum(sum(sum(pop_prior_j_prime_6D,2),3),4)),'DisplayName','Update C Prior','LineWidth',1.5)
                plot(squeeze(sum(sum(sum(pop_prior_j_prime_6D,1),2),4)),'DisplayName','Lure C Prior','LineWidth',1.5)
                plot(squeeze(sum(sum(sum(pop_prior_j_prime_6D,1),2),3)),'DisplayName','FA C Prior','LineWidth',1.5)

                plot(squeeze(sum(sum(sum(dists_j_6D,1),3),4)),'DisplayName','Main C','LineWidth',1.5)
                plot(squeeze(sum(sum(sum(dists_j_6D,2),3),4)),'DisplayName','Update C','LineWidth',1.5)
                plot(squeeze(sum(sum(sum(dists_j_6D,1),2),4)),'DisplayName','Lure C','LineWidth',1.5)
                plot(squeeze(sum(sum(sum(dists_j_6D,1),2),3)),'DisplayName','FA C','LineWidth',1.5)
                legend('Location','Best')
                hold off;
            end


            % We are interested in
            % P(Theta | D^s) =
            %   \sum_j P(Theta,j|D^s)
            % = \sum_j [\rho^s_j P(Theta^s_j|D^s)*P(Theta_{j'})]

            joint{s} = joint{s} + (rho(s,j) .* dists_j_6D .* pop_prior_j_prime_6D);
            % Sum over all models j
            % THIS MIGHT BE WRONG - MOVE RHO TO EXTERNAL PART OF CODE
            % WHERE SUBJECTS GET PUT INTO NFC GROUPS?

        end % of excluding s, j combinations with 0 as rho value

    end % end of model loop

    plotflag = false;
    if plotflag
        
        figure(6); 
        for cii = 1:length(top_costs)
            
            margin_list = [1:length(top_costs)];
            margin_list(cii) = [];
            one_marginal = joint{s};
            %one_marginal = nansum(nansum(nansum(nansum(nansum(joint{s},margin_list(5)), ...
                %margin_list(4)),margin_list(3)),margin_list(2)),margin_list(1));
            one_marginal = nansum(nansum(nansum(joint{s},margin_list(3)), ...
                margin_list(2)),margin_list(1));
            plot(xs,squeeze(one_marginal),'DisplayName',top_costs{cii},'LineWidth',1.5)
            hold on;
            
        end
        legend('location','best')
        pause(0.1)
        title(['Posteriors at Subject #' num2str(s)])
        hold off
    end

    non_zero_default = 0.000000000000000000001;
    joint{s}(joint{s}==0) = non_zero_default;
    big_posterior = big_posterior + log(joint{s});

end % end of subject loop

save('joint_cost_distributions_2023.mat','joint','big_posterior','xs','-v7.3')

