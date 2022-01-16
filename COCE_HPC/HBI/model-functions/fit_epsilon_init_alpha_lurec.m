function [llh] = fit_epsilon_init_alpha_lurec(params,data)
%DICTATE_MODEL Dictate which inputs to put into get probs cost learning
%   Then output llh to be used by CBM fitting procedure
global HBI_flag onesim modeltofit 
HBI_flag = true;
onesim = data;
model_name = strrep(mfilename,'fit_','');
modeltofit = coc_createModels(model_name);
[~,~,llh] = getprobs_costlearning(params);
end

