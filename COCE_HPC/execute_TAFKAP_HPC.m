%% Run TAFKAP on data, get estimation accuracy out
function execute_TAFKAP_HPC(subjnum)

% Example subject
%subjlist = [4 5 6];

%for sii = 1:length(subjlist)
    %subjnum = subjlist(sii); 
    % % Grab task data
    load(['task_subj' num2str(subjnum) '.mat'])
    eval(['task = task_subj' num2str(subjnum) ';'])

    %Specify which file you're interested in
    spec = '_megaROIs';

    if contains(ls,['TAFKAP_subj' num2str(subjnum) spec '.mat'])
        % if TAFKAP has been run on this subject already, load that output
        load(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')

    else %TAFKAP has not been run on this subject yet
        % load in BOLD data (already meaned)
        load(['late_delay_mean_PSC_subj' num2str(subjnum) '_VEselected' spec '.mat'])

        BOLD = delay_mean_PSC;
        TAFKAP_output = run_TAFKAP_eowm_HPC(BOLD,task);
        save(['TAFKAP_subj' num2str(subjnum) spec '.mat'], 'TAFKAP_output')

    end % of deciding to run TAFKAP or not

%end
