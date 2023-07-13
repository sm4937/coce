%% Sandbox for trying to get covariances of parameters across WFW models
% Sarah Master & Pat Little 06/05/2023

% load relevant low-level fits, containing A info
fname_hbi = 'different_deltas_versus_winning_models_strict.mat';
N = 100;

%Model fits get saved in a structure called cbm, load that up here and
% grab variables of interest from it.
fits = load(fname_hbi);
high_level_fits = fits.cbm;

% formula for a gaussian
gaussian = @(x,mu,sigma) (1./sqrt(2*pi*sigma))*exp((-(x-mu).^2)/(2*sigma));

% must be fed in IN ORDER
get_new_mean = @(mus,sigmas) (mus.*(1/sigmas))./(sum(1./sigmas));
get_new_variance = @(sigmas) (1./(sum(1./sigmas)));

% get group means for each model, from joint distributions
mus(1) = sum(xs.*mainc_dist);
mus(2) = sum(xs'.*uc_dist);
mus(3) = sum(xs'.*lurec_dist);
mus(4) = sum(xs'.*fac_dist);

low_level_fit_dir = '/Users/sarah/Desktop/HPC_copies/quick_delta_tests_low_level';
list = dir(low_level_fit_dir);
names = char(list.name);

cbm_modellist = high_level_fits.input.fcbm_maps;
costs_of_interest = {'uc','mainc','fac','lurec','missc','respc'};

% initialize variables
cov_mat_weighted = [];
cov_mat_not_weighted = [];

main_models = ['epsilon_initi_alpha_uc.mat              ';...
    'epsilon_initi_alpha_lurec.mat           ';...
    'epsilon_initi_alpha_uc_mainc_lurec.mat  '];
% edit down just to 3 main models in play
names = main_models;

% pull empiricial covariances out from each individual model
for mii = 1:length(names)
    name = names(mii,:);
    rho = 0; r = 0;
    costs_in_this_model = false(1,6);

    modelfile = strsplit(names(mii,:),' ');
    modelidx = find(contains(cbm_modellist,modelfile{1}));
    modelname = strrep(modelfile{1},'.mat','');
    model = coc_createModels(modelname);
    
    load([low_level_fit_dir '/' modelfile{1}])
    low_level_fits = cbm;
    
    model_rho = high_level_fits.output.responsibility(:,modelidx);
    sum_rho = sum(model_rho);
    
    theta_m = high_level_fits.output.group_mean{modelidx};
    
    for cii = 1:6
        costs_in_this_model(cii) = contains(modelfile{1},costs_of_interest{cii});
    end
    
    cost_indices = find(contains(model.paramnames,costs_of_interest(costs_in_this_model)));
    variance_m_c = zeros(length(cost_indices),1);
    variance_m = zeros(1,length(model.paramnames));
    sigma_m_s_ab = zeros(5,5);

    
    for s = 1:N %cycle over subjects to produce variance and co-variance within-model
        
        rho_m_s = model_rho(s); %individual responsibility for this model

        %model.paramnames %params in model, in order
        %theta_m % population parameter means, all parameters
        
        theta_m_s = high_level_fits.output.parameters{modelidx}(s,:); %individual parameter value
        sigma_m_s = diag(inv(low_level_fits.math.A{s}))'; % individual variance over parameter c
        variance_m = variance_m + rho_m_s*((theta_m_s-theta_m).^2 + sigma_m_s);
        
        for cii = 1:length(cost_indices) %cycle over costs specifically
            
            c = cost_indices(cii);
            theta_m_c = theta_m(c); % population parameter mean
            
            theta_m_s_c = high_level_fits.output.parameters{modelidx}(s,c); %individual parameter value
            sigma_m_s_c = inv(low_level_fits.math.A{s}(c,c)); % individual variance over parameter c
            variance_m_c(cii,:) = variance_m_c(cii,:) + rho_m_s*((theta_m_s_c-theta_m_c).^2 + sigma_m_s_c);
            
        end % of cycling over costs
        
        if length(cost_indices) > 1
            % cycle over individual costs to find covariances between them
            
            all_combos = nchoosek(cost_indices,2);
            cov_mat_onemodel = inv(low_level_fits.math.A{s});
            
            for perm = 1:size(all_combos,1)
                % use i,j for cov_mat_onemodel
                i = all_combos(perm,1); j = all_combos(perm,2);
                % use different indices (x,y )for cov_mat
                % (pulled from cost list order up top)
                x = find(contains(costs_of_interest,model.paramnames(i)));
                y = find(contains(costs_of_interest,model.paramnames(j)));
                theta_m_s_a = high_level_fits.output.parameters{modelidx}(s,i);
                theta_m_s_b = high_level_fits.output.parameters{modelidx}(s,j);
                sigma_m_s_ab(x,y) = sigma_m_s_ab(x,y) + rho_m_s * ((theta_m_s_a-theta_m(i))*(theta_m_s_b-theta_m(j)) + cov_mat_onemodel(i,j));
            end % of looping over cost combinations
           
        end % of identifying models with multiple cost parameters
              
    end % of cycling over subjects
    
    variance_m_c = variance_m_c./sum_rho;
    covariance_c = sigma_m_s_ab;
    covariance_c = covariance_c./sum_rho;
    
    cii_index = find(costs_in_this_model);
    % make a print-out describing the variance of each cost in each model
    for cii = 1:length(variance_m_c)
        disp([costs_of_interest{cii_index(cii)} ' model ' num2str(mii) ' variance = ' num2str(variance_m_c(cii,:))])
    end
    
end % cycling over models



%% old obsolete code


% % pull empiricial covariances out from each individual model
% for mii = 1:length(names)
%     name = names(mii,:);
%     rho = 0; r = 0;
%     costs_in_this_model = false(1,6);
%     if contains(name,'.mat') % is this a model file? if yes, continue
%         modelfile = strsplit(names(mii,:),' ');
%         modelidx = find(contains(cbm_modellist,modelfile{1}));
%         modelname = strrep(modelfile{1},'.mat','');
%         model = coc_createModels(modelname); 
%         
%         load([low_level_fit_dir '/' modelfile{1}])
%         low_level_fits = cbm;
%         
%         r = high_level_fits.output.model_frequency(modelidx);
%         rho = high_level_fits.output.responsibility(:,modelidx);
%         
%         for cii = 1:6
%             costs_in_this_model(cii) = contains(modelfile{1},costs_of_interest{cii});
%         end
%         
%         cost_indices = find(contains(model.paramnames,costs_of_interest(costs_in_this_model)));
%         % ennummerate all possible cost combinations
%         if length(cost_indices) > 1 & ~isempty(modelidx)
%             all_combos = nchoosek(cost_indices,2);
%             all_combos = [all_combos; fliplr(all_combos)];
%             
%             for sii = 1:N
%                 cov_mat = NaN(6,6);
%                 rho_i = rho(sii);
%                 cov_mat_onemodel = inv(low_level_fits.math.A{sii});
%                 
%                 for perm = 1:size(all_combos,1)
%                     % use i,j for cov_mat_onemodel
%                     i = all_combos(perm,1); j = all_combos(perm,2);
%                     % use different indices (x,y )for cov_mat 
%                     % (pulled from cost list order up top)
%                     x = find(contains(costs_of_interest,model.paramnames(i)));
%                     y = find(contains(costs_of_interest,model.paramnames(j)));
%                     cov_mat(x,y) = cov_mat_onemodel(i,j);
%                 end
%                 
%                 cov_mat_weighted = cat(3,cov_mat_weighted,rho_i.*cov_mat);
%                 cov_mat_not_weighted = cat(3,cov_mat_not_weighted,cov_mat);
%             end % of looping over subjects
%             
%             %cov_mat_weighted = cov_mat_weighted./sqrt(r);
%             %cov_mat_not_weighted = cov_mat_not_weighted./sqrt(r);
%             % by dividing by sqrt(r) we get S (?) or the group hierarchical
%             % errorbar
%             
%         end % of making sure there is covariance info for 2 costs
%         
%     end % of making sure it's a model file, instead of another file
%     
% end
% 
% 
% costs_of_interest = {'uc','mainc','fac','lurec','missc','respc'};
% cov_mat = nansum(cov_mat_weighted,3)./N;
% cov_mat = round(cov_mat,4);
% weighted = table;
% for cii = 1:length(costs_of_interest)
%     eval(['weighted.' costs_of_interest{cii} ' = cov_mat(:,cii);'])
%     weighted.Properties.RowNames{cii,:} = costs_of_interest{cii};
% end
% weighted;
% 
% cov_mat = nanmean(cov_mat_not_weighted,3);
% cov_mat = round(cov_mat,4);
% not_weighted = table;
% for cii = 1:length(costs_of_interest)
%     eval(['not_weighted.' costs_of_interest{cii} ' = cov_mat(:,cii);'])
%     not_weighted.Properties.RowNames{cii,:} = costs_of_interest{cii};
% end
% not_weighted;
