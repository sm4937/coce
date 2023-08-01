%% Sandbox for trying to get covariances of parameters across WFW models
% Sarah Master & Pat Little 06/05/2023

% THIS SCRIPT GETS RUN WITHIN PAPER_GRAPHS_AND_STATS.M, can't be run on its
% own because it DEPENDS ON OUTPUTS FROM PAPER_GRAPHS_AND_STATS.M


% load relevant low-level fits, containing A info
fname_hbi = 'final_fits_strict.mat';
N = 100;

%Model fits get saved in a structure called cbm, load that up here and
% grab variables of interest from it.
fits = load(fname_hbi);
high_level_fits = fits.cbm;

% some variables, like within-model parameter variances & covariances are
% not available in the high-level model fit files,
% so I'm saving an old version of the fits (single-subject level) and
% pulling from that directory to get these numbers
low_level_fit_dir = 'low-level-fits/';
list = dir(low_level_fit_dir);
names = char(list.name);

cbm_modellist = high_level_fits.input.fcbm_maps;
costs_of_interest = {'uc','mainc','fac','lurec','missc','respc'};

main_models = ['epsilon_initi_alpha_uc.mat              ';...
    'epsilon_initi_alpha_lurec.mat           ';...
    'epsilon_initi_alpha_uc_mainc_lurec.mat  '];
% edit down just to 3 main models in play
names = main_models;

% pull empiricial covariances out from each individual model
% this code does NOT combine these inferences across models, just across
% subjects
for mii = 1:size(names,1)
    name = names(mii,:);
    rho = 0; r = 0;
    costs_in_this_model = false(1,6);

    % in loop, grab model info from model number (m)
    modelfile = strsplit(names(mii,:),' ');
    modelidx = find(contains(cbm_modellist,modelfile{1}));
    modelname = strrep(modelfile{1},'.mat','');
    model = coc_createModels(modelname);
    
    load([low_level_fit_dir '/' modelfile{1}])
    low_level_fits = cbm;
    
    % pull responsibility (important for calculating variance & such)
    model_rho = high_level_fits.output.responsibility(:,modelidx);
    sum_rho = sum(model_rho);
    
    %pull mean of the population prior (group-level mean)
    theta_m = high_level_fits.output.group_mean{modelidx};
    
    for cii = 1:6
        costs_in_this_model(cii) = contains(modelfile{1},costs_of_interest{cii});
    end
    
    %which costs are in this model?
    cost_indices = find(contains(model.paramnames,costs_of_interest(costs_in_this_model)));
    
    % initialize math variables
    variance_m_c = zeros(length(cost_indices),1);
    variance_m = zeros(1,length(model.paramnames));
    sigma_m_s_ab = zeros(5,5);

    %cycle over subjects to produce variance and co-variance within-model
    for s = 1:N 
        
        %individual subject model responsibility
        rho_m_s = model_rho(s); 

        %model.paramnames %params in model, in order
        %theta_m % population parameter means, all parameters
        
        % from high-level model fits, pull individual subject (s) values
        theta_m_s = high_level_fits.output.parameters{modelidx}(s,:); %individual parameter value
        sigma_m_s = diag(inv(low_level_fits.math.A{s}))'; % individual variance over parameter c
        variance_m = variance_m + rho_m_s*((theta_m_s-theta_m).^2 + sigma_m_s);
        
        % % CALCULATE VARIANCE WITHIN-MODEL for EACH COST SEPARATELY % % 
        % each cost has different index and different outputs within the
        % model and subject level
        for cii = 1:length(cost_indices) %cycle over costs specifically
            
            c = cost_indices(cii);
            theta_m_c = theta_m(c); % population parameter mean
            
            % here is where I'm pulling cost-specific, subject-specific,
            % and model-specific values, and doing the math on them to
            % combine them across subjects 
            theta_m_s_c = high_level_fits.output.parameters{modelidx}(s,c); %individual parameter value
            sigma_m_s_c = inv(low_level_fits.math.A{s}(c,c)); % individual variance over parameter c
            variance_m_c(cii,:) = variance_m_c(cii,:) + rho_m_s*((theta_m_s_c-theta_m_c).^2 + sigma_m_s_c);
            
        end % of cycling over costs
        
        
        % % CALCULATE COVARIANCE WITHIN-MODEL and WITHIN-SUBJECT ACROSS COSTS! % % 
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
    
    % normalize by total model responsibility
    variance_m_c = variance_m_c./sum_rho;
    covariance_c = sigma_m_s_ab;
    covariance_c = covariance_c./sum_rho;
    
    % make a print-out describing the variance of each cost in each model
    cii_index = find(costs_in_this_model);
    for cii = 1:length(variance_m_c)
        disp([costs_of_interest{cii_index(cii)} ' model ' num2str(mii) ' variance = ' num2str(variance_m_c(cii,:))])
    end
    
end % of cycling over models


%% old obsolete code

% % initialize variables
% cov_mat_weighted = [];
% cov_mat_not_weighted = [];

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
