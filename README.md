# the costs of cognitive effort (COCE)
the majority of this is Matlab analysis and modeling code, written to analyze data collected on Amazon Mechanical Turk using the task contained within the /cost folder.
The task code is javascript via JsPsych and custom JsPsych plugins.

in /cost, the experimental code (including images, custom plugins, and pre-built jspsych functions)
cost.html is the main experimental page (initializes all, loads up entire experimental timeline)
Index.html is an obsolete example for a main experiment page (the old version of cost.html)
There is also code for collecting Need for Cognition and Perfectionism scores with the Need for Cognition Scale (NFC) and Short Almost Perfect Scale (SAPS).

in /data_loading_and_scoring, I have a bunch of functions I wrote for loading in my mturk data from .json format, or otherwise formatting data. It has code for scoring the NFC and SAPS questionnaires. It also has code for formatting every subject's data into the same format, from the completely unprocessed MTurk files, for ease of comparison on the group level, stats running, etc.

in /example_data, there are processed data files, including the file toAnalyze,
which is a table full of values important for running model fitting within /modeling.
the version provided in the public repository contains example data from one batch
worth of subjects. toAnalyze is created within modeling/run_simple_sims.m

the outer level of /modeling has modeling functions that are generally necessary for running the cost learning/changing models described in the paper. they are used both in the HBI & type II ML model fitting directories, as they are the core of my models. you'll see from looking inside simulate_cost_model.m or getprobs_costlearning.m, that the scripts are written to be super flexible to the model specified such that they can handle many situations, including different parameter values, different included parameters (alpha vs. delta, for example), different fitting algorithms, different subjects, etc. simulate_cost_model.m simulates data given specific parameter values/model specifications, while getprobs_costlearning takes in the model & subject data and returns the probability of that subject's choices given those parameter values & that model. getprobs_costlearning.m can provide log likelihoods, negative log likelihoods, and maximum a priori values for each parameter set it's given.

in /modeling/HBI, code for running the Hierarchical Bayesian Inference model fitting package by Piray & Daw. Included in this folder is multiple files per model which describe how well the model fits real subject data, as well as simulated subject data. These files are referenced by the outer shell HBI_coc.m, which can load any model you specify to display how well it fits with the current model settings.

in modeling/typeII_ML_fitting is code I wrote for running type II maximum likelihood (expectation maximization (EM) algorithm), incorporating group-level priors into individual level parameter fits. this is not the method we ended up using, so it's not entirely up to date with the final results, but it may be helpful to someone... so I'm including it. can't promise it's perfect or all that readable.

in /plotting I have some functions for making nice figures. there are two plugins in there, one called superbar (from https://www.mathworks.com/matlabcentral/fileexchange/57499-superbar), and the other called violin (https://www.mathworks.com/matlabcentral/fileexchange/45134-violin-plot).

there is code inside /obsolete which could be used to re-process data from old versions of this task. the published (fingers crossed) version is the final one, v4:
v1: cost of control with the BDM, a detection task, plain n-back training, and a combined n-back and detection task (based on Kool & Botvinick 2014)
v2: cost of control with the BDM, and two versions of an n-switch task
v3: cost of control with the BDM, and 1- and 2-back tasks
v4: cost of control with the BDM, a 1- and 2-back task, and a 3-detect task (final version with N = 100 subjects collected)
