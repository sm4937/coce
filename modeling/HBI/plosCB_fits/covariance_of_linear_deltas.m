%% Sandbox for trying to get covariances of parameters across WFW models
% Sarah Master & Pat Little 06/05/2023

% THIS SCRIPT GETS RUN WITHIN PAPER_GRAPHS_AND_STATS.M, can't be run on its
% own because it DEPENDS ON OUTPUTS FROM PAPER_GRAPHS_AND_STATS.M


% load relevant low-level fits, containing A info
fname_hbi = 'twodelta_analyses.mat';
N = 100;

%Model fits get saved in a structure called cbm, load that up here and
% grab variables of interest from it.
fits = load(fname_hbi);
high_level_fits = fits.cbm;

% some variables, like within-model parameter variances & covariances are
% not available in the high-level model fit files,
% so I'm saving an old version of the fits (single-subject level) and
% pulling from that directory to get these numbers
low_level_fit_dir = '../low-level-fits/';

cbm_modellist = high_level_fits.input.fcbm_maps;
costs_of_interest = {'delta_time','delta_practice'};

main_models =  ['epsilon_initi_twodeltas_lurec.mat'; ...
'epsilon_initi_twodeltas_mainc.mat'; ...
'epsilon_initi_twodeltas_uc.mat   '];
% edit down just to 3 main models in play
names = main_models;

% pull empiricial covariances out from each individual model
% this code does NOT combine these inferences across models, just across
% subjects
for mii = 1:size(names,1)
    name = names(mii,:);
    rho = 0; r = 0;

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
    
    %which costs are in this model?
    cost_indices = find(contains(model.paramnames,costs_of_interest));
    
    % initialize math variables
    variance_m_c = zeros(length(cost_indices),1);
    variance_m = zeros(1,length(model.paramnames));
    sigma_m_s_ab = 0;

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
            
            cov_mat_onemodel = inv(low_level_fits.math.A{s});
            
            x = cost_indices(1);
            y = cost_indices(2);
            theta_m_s_a = high_level_fits.output.parameters{modelidx}(s,x);
            theta_m_s_b = high_level_fits.output.parameters{modelidx}(s,y);
            sigma_m_s_ab = sigma_m_s_ab + rho_m_s * ((theta_m_s_a-theta_m(x))*(theta_m_s_b-theta_m(y)) + cov_mat_onemodel(x,y));
            
        end % of identifying models with multiple cost parameters
        
    end % of cycling over subjects
    
    % normalize by total model responsibility
    variance_m_c = variance_m_c./sum_rho;
    covariance_c = sigma_m_s_ab;
    covariance_c = covariance_c./sum_rho;
    
    % make a print-out describing the variance of each cost in each model
    for cii = 1:2
        disp([costs_of_interest{cii} ' model ' num2str(mii) ' variance = ' num2str(variance_m_c(cii,:))])
    end
    
end % of cycling over models

