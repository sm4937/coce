# the costs of cognitive effort (COCE)
the majority of this is Matlab analysis and modeling code, written to analyze data collected on Amazon Mechanical Turk using the task contained within the /cost folder.
The task code is javascript via JsPsych and custom JsPsych plugins.

in /cost, the experimental code (including images, custom plugins, and pre-built jspsych functions)
cost.html is the main experimental page (initializes all, loads up entire experimental timeline)
There is also code for collecting Need for Cognition and Perfectionism scores with the Need for Cognition Scale (NFC) and Short Almost Perfect Scale (SAPS).

in /simdata, there are processed data files, including the file toAnalyze,
which is a table full of values important for running model fitting within /modelling.
the version provided in the public repository contains example data from one batch
worth of subjects. toAnalyze is created within modeling/run_simple_sims.m

in /HBI, code for running the Hierarchical Bayesian Inference model fitting package by Piray & Daw. Included in this folder is multiple files per model which describe how well the model fits real subject data, as well as simulated subject data. These files are referenced by the outer shell HBI_coc.m, which can load any model you specify to display how well it fits with the current model settings.



there is code inside /obsolete which could be used to re-process data from old versions of this task. the published (fingers crossed) version is the final one, v4:
v1: cost of control with the BDM, a detection task, plain n-back training, and a combined n-back and detection task (based on Kool & Botvinick 2014)
v2: cost of control with the BDM, and two versions of an n-switch task
v3: cost of control with the BDM, and 1- and 2-back tasks
v4: cost of control with the BDM, a 1- and 2-back task, and a 3-detect task (final version with N = 100 subjects collected)
