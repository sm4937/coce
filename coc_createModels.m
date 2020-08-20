function [model] = coc_createModels(name)
%coc_createModels Create model specifications for cost of control
%   Outlines parameters, special functionality, other things that may
%   distinguish models from one another
%   For example, takes in 'alpha_epsilon' and outputs model structure with
%   number of parameters, and which ones they are.

modelStruct = struct;

modelStruct.uc.uc = true;
modelStruct.uc.epsilon = false;
modelStruct.uc.init = false;
modelStruct.uc.mc = false;
modelStruct.uc.mainc = false;
modelStruct.uc.alpha = false;
modelStruct.uc.nparams = 1;
modelStruct.uc.paramnames = {'uc'};

modelStruct.uc_epsilon.uc = true;
modelStruct.uc_epsilon.epsilon = true;
modelStruct.uc_epsilon.init = false;
modelStruct.uc_epsilon.mc = false;
modelStruct.uc_epsilon.mainc = false;
modelStruct.uc_epsilon.alpha = false;
modelStruct.uc_epsilon.nparams = 2;
modelStruct.uc_epsilon.paramnames = {'uc','epsilon'};

modelStruct.uc_epsilon_init.uc = true;
modelStruct.uc_epsilon_init.epsilon = true;
modelStruct.uc_epsilon_init.init = true;
modelStruct.uc_epsilon_init.mc = false;
modelStruct.uc_epsilon_init.mainc = false;
modelStruct.uc_epsilon_init.alpha = false;
modelStruct.uc_epsilon_init.paramnames = {'uc','epsilon','init'};
modelStruct.uc_epsilon_init.nparams = length(modelStruct.uc_epsilon_init.paramnames);

modelStruct.uc_init.uc = true;
modelStruct.uc_init.epsilon = false;
modelStruct.uc_init.init = true;
modelStruct.uc_init.mc = false;
modelStruct.uc_init.mainc = false;
modelStruct.uc_init.alpha = false;
modelStruct.uc_init.paramnames = {'uc','init'};
modelStruct.uc_init.nparams = length(modelStruct.uc_init.paramnames);

%UC = update cost
%epsilon = noise in BDM reporting
%init = initial cost of each task

modelStruct.uc_epsilon_init_mc.uc = true;
modelStruct.uc_epsilon_init_mc.epsilon = true;
modelStruct.uc_epsilon_init_mc.init = true;
modelStruct.uc_epsilon_init_mc.mc = true;
modelStruct.uc_epsilon_init_mc.mainc = false;
modelStruct.uc_epsilon_init_mc.alpha = false;
modelStruct.uc_epsilon_init_mc.nparams = 4;
modelStruct.uc_epsilon_init_mc.paramnames = {'uc','epsilon','init','mc'};

modelStruct.uc_init_mc.uc = true;
modelStruct.uc_init_mc.epsilon = false;
modelStruct.uc_init_mc.init = true;
modelStruct.uc_init_mc.mc = true;
modelStruct.uc_init_mc.mainc = false;
modelStruct.uc_init_mc.alpha = false;
modelStruct.uc_init_mc.nparams = 3;
modelStruct.uc_init_mc.paramnames = {'uc','init','mc'};

modelStruct.uc_epsilon_init_mainc.uc = true;
modelStruct.uc_epsilon_init_mainc.epsilon = true;
modelStruct.uc_epsilon_init_mainc.init = true;
modelStruct.uc_epsilon_init_mainc.mc = false;
modelStruct.uc_epsilon_init_mainc.mainc = true;
modelStruct.uc_epsilon_init_mainc.alpha = false;
modelStruct.uc_epsilon_init_mainc.nparams = 4;
modelStruct.uc_epsilon_init_mainc.paramnames = {'uc','epsilon','init','mainc'};

modelStruct.uc_init_mainc.uc = true;
modelStruct.uc_init_mainc.epsilon = false;
modelStruct.uc_init_mainc.init = true;
modelStruct.uc_init_mainc.mc = false;
modelStruct.uc_init_mainc.mainc = true;
modelStruct.uc_init_mainc.alpha = false;
modelStruct.uc_init_mainc.nparams = 3;
modelStruct.uc_init_mainc.paramnames = {'uc','init','mainc'};

modelStruct.uc_epsilon_init_mc_mainc.uc = true;
modelStruct.uc_epsilon_init_mc_mainc.epsilon = true;
modelStruct.uc_epsilon_init_mc_mainc.init = true;
modelStruct.uc_epsilon_init_mc_mainc.mc = true;
modelStruct.uc_epsilon_init_mc_mainc.mainc = true;
modelStruct.uc_epsilon_init_mc_mainc.alpha = false;
modelStruct.uc_epsilon_init_mc_mainc.nparams = 5;
modelStruct.uc_epsilon_init_mc_mainc.paramnames = {'uc','epsilon','init','mc','mainc'};

modelStruct.uc_epsilon_init_mc_mainc_alpha.uc = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.epsilon = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.init = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.mc = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.mainc = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.alpha = true;
modelStruct.uc_epsilon_init_mc_mainc_alpha.nparams = 6;
modelStruct.uc_epsilon_init_mc_mainc_alpha.paramnames = {'uc','epsilon','init','mc','mainc','alpha'};

eval(['model = modelStruct.' name])

end

