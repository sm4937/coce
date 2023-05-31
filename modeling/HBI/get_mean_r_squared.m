function [mean_r_squared] = get_mean_r_squared(nboot,subj_list,model_list,params,toanalyze)
% takes in model (or models) to simulate with, one per subject
% simulates dataset nboot number of times
% correlates all subject data each time
% then means that, at the end of the bootstrapping bit

disp('Calculating mean r-squared iteratively')

nsubjs = length(subj_list);

for boot = 1:nboot

    if mod(boot,10)==0
        disp(['Iter # ' num2str(boot)])
    end
    
    % %simulate data from fit parameters, & specified model
    % (model_list{subj})
    % go subject by subject because each has own model, own fit parameters
    simwages = []; truewages = [];
    for s = 1:nsubjs
        subj = subj_list(s);
        onesubj = toanalyze(toanalyze.subj == subj,:);
        %then simulate their data with that model
        subj_model = model_list{subj};
        % input to simulate_cost_model the model specs, then the transformed
        % parameters for that subject, then that subject's data
        
        simdata = simulate_cost_model(subj_model,params{subj},onesubj);
        % NOTE! PARAMETERS ALREADY TRANSFORMED OUTSIDE THIS FUNCTION!
        
        %save in one long vector, for all subjects' sim data
        simwages = [simwages; (simdata(:,3)./25)+1;  NaN(32-length(simdata(:,3)),1)]; %cushion for those trials at the end of
       
        % some subjects' data where the task identity wasn't saved for
        % whatever reason so ratings have no accompanying task
        % love running online experiments w/ weird internet connection bugs
        % :)
        
        truewages = [truewages; toanalyze.BDM(toanalyze.subj==subj,:)];
    end
    
    % clean up NaNs
    todelete = sum(isnan(simwages),2)>0 | sum(isnan(truewages),2)>0;
    simwages(todelete,:) = [];
    truewages(todelete,:) = [];
    
    [r,p] = corr(simwages,truewages);
    r_squared(boot,:) = r^2;
end

mean_r_squared = nanmean(r_squared,1);
% mean down the vector, produce mean r^2 over all boots of algorithm

end
