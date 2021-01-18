function [modeltofit] = EMfit_sm(data,modeltofit,starting)
%EMfit_sm Run parameter and model score fitting procedure
%   data is organized like this:
%   1. subj num, 2. nmatches, 3. maintained, 4. nupdates
%   5. nmisses, 6. task displayed, 7. task executed, 8. BDM, 9. accurate
%   BDM for that task

%create global variables
global onesim fit_epsilon_opt mu sigma realsubjectsflag scalingvector canbeneg realparamlist
scalingvector = ones(1,modeltofit.nparams); canbeneg = false(1,modeltofit.nparams);

% set fmincon options, settings for model fitting
options = optimoptions('fmincon','Display','off');
nparams = modeltofit.nparams;
nsubjs = size(unique(data(:,1)),1);

%recover!
pmin=-5*ones(1,nparams);pmax=5*ones(1,nparams); %keep within reasonable ranges
if modeltofit.epsilon; idx = find(contains(modeltofit.paramnames,'epsilon')); pmin(idx) = 0; pmax(idx) = Inf; end %epsilon cannot be negative, should be unbounded at top
if modeltofit.alpha; idx = find(contains(modeltofit.paramnames,'alpha')); pmin(idx) = 0; pmax(idx) = 1; end %set alpha specially
%appliedscaling = scalingvector(1:nparams); 
appliedneg = canbeneg(1:nparams); %scale to current model
%starting with just alpha
bestmu = rand(1,nparams); bestsigma = rand(1,nparams); %start from random low numbers
highlconvergence = []; modelscores = [];

maxiters = 40; iters = 0; convergence = false; converge_threshold = 0.005; stability = 4;
while iters <= maxiters & ~convergence %allow the algorithm to converge early
    iters = iters + 1
    mus = rand(10,nparams); sigmas = rand(10,nparams); mus(1,:) = bestmu; sigmas(1,:) = bestsigma;
    %mus = mus.*appliedscaling; sigmas = sigmas.*appliedscaling; 
    chancevector = rand(10,nparams)<0.5;
    mus(chancevector&repmat(appliedneg,10,1)) = -mus(chancevector&repmat(appliedneg,10,1)); %throw some negative starting points in there too
    sigmas(sigmas<0.05)=0.05; %don't want to run into infinitely tall/0 width param distributions
    bestdist = []; paramsbytry = []; smbytry = []; 
    for tries = 1:10 %of different mus/sigmas
        bestforsubj = [];
        mu = mus(tries,:); sigma = sigmas(tries,:);
        dbstop in EMfit_sm.m at 40 if sum(isinf(sigma))>0
        for subj = 1:nsubjs %go over each subject
            if istable(data)
                idx = data.subj == subj;
            else
                idx = data(:,1) == subj;
            end
            onesim = data(idx,:);
            fitparams = [];
            for z = 1:5 %z random starts for each subject
                params = normrnd(mu,sigma); %initialization points, random
                if z == 1 & (size(bestforsubj,1)==nsubjs)
                    params = bestforsubj(subj,2:2+(nparams-1));
                    % first guess: last best fitting values, if they exist
                end
                if modeltofit.epsilon; idx = find(contains(modeltofit.paramnames,'epsilon'));
                   if params(idx)<0; params(idx) = 0.05; end
                end %bound epsilon
                % need to bound param guesses to what they're permitted to
                % be
                if modeltofit.alpha; idx = find(contains(modeltofit.paramnames,'alpha'));
                    if params(idx) < 0; params(idx) = 0;
                    elseif params(idx) > 1; params(idx) = 1;end
                end %alpha cannot be negative or above 1 either
                % grab parameters from less restrictive model if already
                % fit, start from there
                if ~isempty(starting)&iters==1
                    params = starting(subj,:);
                    missing = abs(nparams-size(starting,2));
                    params = [params zeros(1,missing)];
                end
                try
                    [bestparam,map,~,~,~,~,hess]=fmincon(@getprobs_costlearning,params,[],[],[],[],pmin,pmax,[],options);
                catch
                    bug_in_inside_flag = true
                end
                secondmoment = diag(inv(hess))';
                fitparams = [fitparams; subj z bestparam secondmoment map]; %for the iters over subject, which iters produced what
            end
            [bestmap,which] = min(fitparams(:,end));
            bestforsubj = [bestforsubj; subj fitparams(which,3:3+(nparams-1)) fitparams(which,3+nparams:end-1) bestmap];
            %track subj num, params, sigmas, map score
        end
        %re-infer mu and sigma?
        currentparams = bestforsubj(:,2:2+(nparams-1));
        currentsecondmoments = bestforsubj(:,(2+nparams):end-1);
        paramsbytry{tries} = currentparams;
        smbytry{tries} = currentsecondmoments;
        mapsubjs = mean(bestforsubj(:,end)); %multiply all subj maps or add all llh's? Or mean, to avoid outlier problems?
        bestdist = [bestdist; mus(tries,:) sigmas(tries,:) mapsubjs];
    end %of cycling through random/best mu and sigma starting points
    [MAP,which] = min(bestdist(:,end));
    bestmu = mean(paramsbytry{which}); 
    bestsigma = sqrt(mean(paramsbytry{which}.^2 + smbytry{which})-(bestmu.^2)); %variance
    modelscores = [modelscores; MAP];
    highlconvergence(iters,:) = bestmu; movement = abs(diff(highlconvergence));
    if iters>stability
        convergence = sum(sum(movement((end-stability+1):end,:)<converge_threshold))==nparams*stability; %threshold of a tiny change in mu allows algorithm to stop
    end
    figure(1)
    plot(1:size(highlconvergence,1),highlconvergence,'DisplayName','Convergence of mu')
    xlabel('Iteration')
    fig = gcf; fig.Color = 'w';
    legend(modeltofit.paramnames)
    title('Convergence of high-level parameters')
end %of trying to find the best mu and sigma
%score model fit for all subjects individually
for subj = 1:nsubjs
    mu = bestmu; sigma = bestsigma;
    if istable(data)
        idx = data.subj == subj;
    else
        idx = data(:,1) == subj;
    end
    onesim = data(idx,:);
    params = bestforsubj(:,2:2+(nparams-1));
    [bestfit(subj),group_llhs(subj)] = getprobs_costlearning(params(subj,:));
end

[modeltofit.map] = min(modelscores);
modeltofit.lowparams = paramsbytry{which}; %bestforsubj(:,2:2+(nparams-1));
modeltofit.highparams = [bestmu bestsigma];
modeltofit.fitbysubj = bestfit;
modeltofit.secondmoments = smbytry{which};
modeltofit.convergence = convergence;
% If still developing model, makes sure your model fits from your true
% param values are the minimum in terms of llh.
figure
plot(1:size(highlconvergence,1),highlconvergence,'DisplayName','Convergence of mu')
xlabel('Iteration')
fig = gcf; fig.Color = 'w';
legend(modeltofit.paramnames)
title('Convergence of high-level parameters')
end

%% Tests to put back in if models fit poorly (i.e. you get weird hessian values and therefore weird sigmas)
%dbstop in MLEM_sm.m at 161 if sum(~isreal(bestsigma))>0
% while sum(~isreal(bestsigma))>0 %ill-fitting?
%     thrownout = [thrownout; MAP];
%     bestdist(which,:) = []; %erase last minimum, poorly fitting 
%     idxes = 1:tries; idxes(which) = [];
%     [MAP,which] = min(bestdist(:,end)); %replace with second-best
%     bestmu = mean(paramsbytry{idxes(which)});
%     bestsigma = sqrt(mean(paramsbytry{idxes(which)}.^2 + smbytry{idxes(which)})-(bestmu.^2)); %variance
% end
