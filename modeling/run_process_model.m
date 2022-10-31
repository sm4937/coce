%% Load in subject data for use in simulations
% this code stands alone
% it loads subject data all on its own, formats it, runs it through our
% simple WM/attention process model, and outputs the operations that the
% model uses to complete each task.
% these cost components are then fed into the modelling code in /modelling,
% through loading toanalyze.mat, a file that gets saved at the end of the
% second section of this script.

% in addition, there's some analysis in here of the overall distributions
% of the cost components being extracted in this script (seeing whether
% they are correlated, normally distributed, etc., etc.)
% this level of analysis is more important on the computational modeling
% side of things, and doesn't aim to understand subject behavior, hence it
% stands alone here & is not tied to behavioral analysis in
% paper_graphs_and_stats.m or run_supplementary_analyses.m.

clear all

if isdir('../data')
    % ALL experimental data, all 100 subjects live in 'data'
    load('../data/filenames.mat')
    load('../data/fullsubjnumbers.mat')
    data_directory = '../data/';
    % grab all subjects from those files, then
elseif isdir('../example_data')
    files{1} = 'example_subjs.mat';
    load('../example_data/fullsubjnumbers.mat')
    data_directory = '../example_data/';
end
% specify behavioral files collected on version 4 of WFW behavioral
% task

%grab usable subjects, pre-saved in paper_graphs_and_stats.m
% (subject numbers are in "list", which are randomly generated at the beginning of
% the MTurk task & about 7 digits long)

preloadflag = false;
% do you want to run all these analyses? they take some time

if preloadflag
    
    tosim = [];
    % initialize task data saving variable, tosim
    % to sim is an input to the process model which is run in the cell
    % below
    nsubjs_perfile = zeros(length(files),1);
    % identify which file is good to present as example data
    
    for file = 1:length(files) %loop over behavioral files
        
        load([data_directory files{file}])
        clear tasklist oneperson
        % load stuff up, clean workspace for each subject
        
        for row = 1:length(list) % loop over subjects
            subj = list(row);
            % grab subject number
            data = Untitled(Untitled.subjnum==subj,:);
            % grab that subject's data, if it's in the current file
            
            if height(data)>0 %that subj might not be in that file
                % if it is, process it into usable form for process model
                % below
                nsubjs_perfile(file) = nsubjs_perfile(file)+1;
                data = sortrows(data,find(ismember(data.Properties.VariableNames,'trial_index')));
                % subject data are sometimes saved out of order, but
                % sorting by "trial_index" variable fixes this issue
                % must be a PHP or JsPsych bug in the original task code
                
                oneperson = table;
                % table for saving data from this one subj
                
                display = double(string(data.task_displayed)); % I have to put this early
                % because subj 98 has duplicate values here, an extra trial
                % which should be trimmed at the start
                display(end-15:end) = []; %trim post-task difficulty rating trials
                if sum(~isnan(display))>32 %one subj has 33 display trials due to rare data saving problem
                    duplicate = find(diff(data.trial_index(~isnan(data.trial_index)))==0);
                    % there's duplicate trial in the mix
                    data(duplicate+1,:) = [];
                    % remove it
                end
                
                tasknum = double(string(data.tasknum));
                idx = ~isnan(tasknum); %clean this up a little before sim stuff
                % this is a KEY VARIABLE ^^
                % this indexes all task trials, not BDM trials, not
                % questionnaires, etc.
                % this gets used to trim the FULL data file for each
                % subject into just the data relevant to
                % modeling/simulation work (i.e. task key presses, which
                % stim is on screen when...)
                
                oneperson.block = tasknum(idx);
                oneperson.subj = repmat(row,height(oneperson),1);
                oneperson.n = data.n(idx,:);
                oneperson.n = double(string(oneperson.n)); %account for "diverse" data types
                oneperson.n(1:10) = 0; %the first 10 trials are 0-back practice
                % one subj has a missing n-back separation trial for some
                % reason
                oneperson.detect = data.detect(idx)==categorical(cellstr('1'));
                oneperson.nback = data.nback(idx)==categorical(cellstr('1'));
                values = double(string(data.value));
                display = double(string(data.task_displayed));
                display(end-15:end) = []; %trim post-task difficulty rating trials
                % trim practice BDMs from here, too
                oneperson.BDM = NaN(height(oneperson),1); oneperson.display = NaN(height(oneperson),1);
                % just lots of cleanup
                
                % find where one block ends and the next begins
                d = sum([false;diff(oneperson.block)>0])-sum(~isnan(values));
                valuestoinsert = [values(~isnan(values)); NaN(d,1)];
                d = sum([false;diff(oneperson.block)>0])-sum(~isnan(display));
                displaytoinsert = [display(~isnan(display)); NaN(d,1)];
                
                if sum(abs(diff(oneperson.n))>0)>3 %find subjects with missing key values from BDM trials (1 each)
                    missing = find(abs(diff(oneperson.n))>0); missing = missing(end); missing = oneperson.block(missing+1); %ignore practice block
                    old = valuestoinsert; valuestoinsert = old(1:missing-1,:);
                    valuestoinsert(missing) = NaN; valuestoinsert = [valuestoinsert; old(missing+1:32)];
                    old = displaytoinsert; displaytoinsert = old(1:missing-1,:);
                    displaytoinsert(missing) = NaN; displaytoinsert = [displaytoinsert; old(missing+1:32)];
                end
                
                oneperson.BDM([false;diff(oneperson.block)>0]) = valuestoinsert;
                oneperson.display([false;diff(oneperson.block)>0]) = displaytoinsert;
                oneperson.display(oneperson.display==7) = 3;
                oneperson.BDM = (oneperson.BDM-1)*25; %from 0 to 100, to clarify math stuff
                oneperson.task = NaN(height(oneperson),1);
                tasklist = data.task(idx);
                oneperson.task(tasklist==categorical(cellstr('detection'))) = 0;
                oneperson.task(tasklist==categorical(cellstr('n-back'))&oneperson.n==1) = 1;
                oneperson.task(tasklist==categorical(cellstr('n-back'))&oneperson.n==2) = 2;
                oneperson.task(tasklist==categorical(cellstr('ndetection'))) = 3;
                oneperson.correct = false(height(oneperson),1);
                oneperson.correct = data.correct(idx)==categorical(cellstr('1'));
                oneperson.stimnum = data.stimnum(idx);
                oneperson.nmatches = double(string(data.nmatches(idx)));
                oneperson.nmatches(isnan(oneperson.nmatches)) = 0; %replace NaN's with 0
                oneperson.keypress = data.key_press(idx);
                tosim = [tosim; oneperson];
                % save all of the data you've pulled from
            end
        end
    end %end of loading/formatting
    
    save([data_directory 'simdata.mat'],'tosim')
    
else % if NOT running all of this stuff, and just preloading, then:
    
    load([data_directory 'simdata.mat'])
    
end

%% Pull possible costs from subjects' behavior and task presented to them
% Goes through the tasks round by round to tally the costs incurred by each
% subject on each round. Does it according to task rules, but assuming
% an internal process which is common across tasks.
% This is what we call our SIMPLE PROCESS MODEL!

practiceblocks = tosim(tosim.block==-1,:);
% Sadly this is not relevant yet, as many of my subjects' practice data
% were not saved trial-by-trial as is necessary for running the process
% model. One day it will come in handy for further iterations of this task.
% But for now, it's irrelevant.

components = []; trim = []; %initialize empty variables

noise = 0; %turn this to a value that is NOT 0 to see how cost components change with
% theoretical WM storage noise.
% for now, we're just registering whether people forget in their behavioral
% responses, not fitting any free parameters from within their process
% models (like WM storage noise)
% this could be an expansion of our framework

for subj = 1:length(unique(tosim.subj)) % loop over subjects captured in tosim
    oneperson = tosim(tosim.subj==subj,:);
    completed(subj,:) = ceil([sum(oneperson.task==0) sum(oneperson.task==1) sum(oneperson.task==3) sum(oneperson.task==2)]./15);
    for block = 0:max(oneperson.block)
        blockdata = oneperson(oneperson.block==block,:);
        nupdates = 0; nlures = 0; %initialize cost variables
        maintained = NaN(1,2); %keep track of storage over time, to get a sense of maintained info over time
        task = blockdata.task(2)+1;
        n = blockdata.n(2); matches = blockdata.nmatches(2);
        if isnan(n) %again, with this duplicate trial for subj 98, causing problems
            n = blockdata.n(3); matches = blockdata.nmatches(3);
        end
        storage = NaN(1,n); %empty storage for n-back
        for trial = 2:height(blockdata)
            stim = blockdata.stimnum(trial);
            if n == 0 %do detection, pretty simple
                maintained(trial,:) = 0;
            else %n-backs and n-detects
                if task==3 %only in 2-back can there be WM lures inside buffer
                    if storage(2)==stim
                        nlures = nlures + 1; %not a real 2-back, but a lure trial
                    end
                end
                if sum(stim == storage)==n %full storage remaining the same
                    % do nothing, no cost, no update
                else %update storage with some noise
                    storage(1) = [];
                    die = rand(); %get a random number
                    if die < noise %noise param wins
                        storage(end+1) = NaN;
                    else %noise param doesn't win
                        storage(end+1) = stim; %clear out old info, add in new info, shift stuff over
                    end
                    nupdates = nupdates + 1; %cost of update/gating in/out
                end
                maintained(trial,:) = storage;
            end
        end
        % grab noisiness for whole block instead of trial
        noisiness = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==false)/max([sum(blockdata.detect|blockdata.nback) 1]);
        nmatches = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==true);
        nmisses = sum((blockdata.detect|blockdata.nback)==true&blockdata.correct==false);
        % count times they were supposed to respond & didn't
        nFAs = sum((blockdata.detect|blockdata.nback)==false & ~isnan(blockdata.keypress));
        nerrors = nmisses + nFAs;
        % count times they responded & weren't supposed to (false alarm)
        nresponses = sum(~isnan(blockdata.keypress(2:end)));
        maintained(1,:) = []; maintained = nanmean(sum(maintained>0,2));
        % this is a weird measure, it's maintainence adjusted for errors,
        % this should be confounded with errors
        % let's just set it to N and see if that's different in any way
        maintained = n;

        components = [components; subj nmatches maintained nupdates nmisses blockdata.display(1)+1 task blockdata.BDM(1) noisiness nresponses nlures nerrors nFAs];
    end %end of block by block loop

    %     if sum(completed(subj,2:end)>2)<3 %3 iterations or more for each rated task?
    %         disp(['excluding subj ' num2str(subj)])
    %         components(components(:,1)==subj,:) = [];
    %         trim = [trim;subj];
    %     end
    % before our model fitting was fully hierarchical, I included code here to
    % exclude subjects with low amounts of task performance on one or more
    % tasks ( this is the if sum(completed...) block above ). now that the
    % model fitting is fully hierarchical, I have removed this bit of code,
    % because the nature of the model fitting allows for subjects with low
    % amounts of task data to still be accurately fit by the models.

end %of loop over subjects
%completed

% save it all in toanalyze.mat
toanalyze = table;
toanalyze.subj = components(:,1); toanalyze.nmatches = components(:,2);
toanalyze.maintained = components(:,3); toanalyze.nupdates = components(:,4);
toanalyze.nmisses = components(:,5); toanalyze.display = components(:,6);
toanalyze.task = components(:,7); toanalyze.BDM = components(:,8);
toanalyze.noisiness = components(:,9); toanalyze.nresponses = components(:,10);
toanalyze.nlures = components(:,11); toanalyze.nerrors = components(:,12);
toanalyze.nFAs = components(:,13);
save([data_directory 'toanalyze.mat'],'toanalyze','trim')

% how many times did the 2-back actually get completed?
% this ALL subject measure is relevant because model-fitting is
% hierarchical
% if all people completely avoided the 2-back, then we'd have a model
% recovery problem. but it looks like there's a good spread in the numbers
% here.
figure
histogram(completed(:,end))
fig = gcf; fig.Color = 'w';
title('Completed 2-backs across subjects')

