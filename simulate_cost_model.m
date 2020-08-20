function [simdata] = simulate_cost_model(modeltosim,allparams,toanalyze)
%simulate_cost_model Simulate specified cost model with specified
%parameters
%Get simdata out, which is plottable and analyzable
    %lower level parameters
    global bounds
    init = 0.25.*bounds(2); epsilon = 0.7.*bounds(2); mc = 0; mainc = 0; alpha = 0.35;
    nparams = modeltosim.nparams; nsubjs = size(allparams,1);
    
    %global simdata %initialize empty
    simdata = [];
    tasks = unique(toanalyze.task);
    for subj = 1:nsubjs  %simulate the model for one subject at a time
    real_epsilon = []; params = allparams(subj,:);
    if modeltosim.uc
        idx = find(contains(modeltosim.paramnames,'uc'));
        uc = params(idx);
        %epsilon = 0.2; %note to self: make sure hist of params is normally distrib
    end
    if modeltosim.epsilon
        idx = find(contains(modeltosim.paramnames,'epsilon'));
        epsilon = params(idx).*bounds(2);
        %epsilon = 0.2; %note to self: make sure hist of params is normally distrib
    end
    if modeltosim.init
        idx = find(contains(modeltosim.paramnames,'init'));
        init = params(idx).*bounds(2);
    end
    if modeltosim.mc
        idx = find(contains(modeltosim.paramnames,'mc'));
        mc = params(idx).*bounds(2);
    end
    if modeltosim.mainc
        idx = find(contains(modeltosim.paramnames,'mainc'));
        mainc = params(idx).*bounds(2);
    end
    if modeltosim.alpha
        idx = find(contains(modeltosim.paramnames,'alpha'));
        alpha = params(idx);
    end
    ratings = init*(ones(1,length(unique(tasks)))); %init for each subject
    onesubj = toanalyze(toanalyze.subj==subj,:);
    ntrials = height(onesubj);
    for trial = 1:ntrials
        task = onesubj.task(trial);
        updates = onesubj.nupdates(trial);%./(max(onesubj.nupdates));
        misses = onesubj.nmisses(trial);%./(max(onesubj.nmisses)); %normalize between 0 and 1
        %misses = ceil(rand().*4);
        mains = onesubj.maintained(trial);%./(max(onesubj.maintained));
        %mains = ceil(rand().*4);
        cost = uc*updates;
        if modeltosim.mc
            cost = cost + mc*misses;
        end
        if modeltosim.mainc
            cost = cost+ mainc*mains;
        end
        if modeltosim.alpha
            ratings(task) = ratings(task) + alpha*(cost-ratings(task)); %delta rule
        else %alpha is fixed
            ratings(task) = ratings(task) + alpha*(cost); %not a learning rate of one
        end
        noiselessri(trial) = ratings(task);
        rating = ratings(task) + normrnd(0,epsilon); %mean 0, std epsilon
        rating(rating<0)=0; rating(rating>bounds(2)) = bounds(2);
        real_epsilon = [real_epsilon; subj rating ratings(task) rating-ratings(task)];
        simdata = [simdata; subj task rating updates uc misses mc mains mainc];
    end %of one subj run-through
    
    real_epsilon_opt(subj) = sqrt((1/ntrials) * sum(real_epsilon(:,4)).^2);
    end %of simulating all subjects
end

