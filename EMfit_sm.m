function [modeltofit] = EMfit_sm(data,modeltofit)
%EMfit_sm Run parameter and model score fitting procedure
%   data is organized like this:
%   1. subj num, 2. nmatches, 3. maintained, 4. nupdates
%   5. nmisses, 6. task displayed, 7. task executed, 8. BDM, 9. accurate
%   BDM for that task

%create global variables
global modeltofit onesim fit_epsilon_opt noiselessri mu sigma realsubjectsflag bounds

% set fmincon options, settings for model fitting
options = optimoptions('fmincon','Display','off');
toexecute = 'normrnd(mu,sigma) ';
nparams = modeltofit.nparams;
nsubjs = size(unique(data(:,1)),1);

%recover!
pmin=zeros(1,nparams);pmax=ones(1,nparams);
niters = 10; 
%starting with just alpha
bestmu = rand(1,nparams); bestsigma = rand(1,nparams);  
highlconvergence = []; modelscores = [];
count = 0; thrownout = 500; %count number of outliers being fit 
for iters = 1:niters %whole EM algorithm
    mus = rand(10,nparams); sigmas = rand(10,nparams); 
    mus(1,:) = bestmu; sigmas(1,:) = bestsigma; 
    sigmas(sigmas<0.05)=0.05;
    bestdist = []; paramsbytry = []; smbytry = []; lowlconvergence = []; bestmus = []; bestsigmas = [];
    for tries = 1:10 %of different mus/sigmas
        bestforsubj = [];
        mu = mus(tries,:); sigma = sigmas(tries,:);
        dbstop in EMfit_sm.m at 33 if sum(isinf(sigma))>0
        for subj = 1:nsubjs
            if istable(data)
                idx = data.subj == subj;
            else
                idx = data(:,1) == subj;
            end
            onesim = data(idx,:);
            fitparams = [];
            for z = 1 %z random starts for each subject
                params = [];
                for i = 1:nparams
                    params(i) = normrnd(mu(i),sigma(i));
                end
                params(params<0) = 0; params(params>1) = 1;
                %params = realparamlist(subj,:);

                [bestparam,map,~,~,~,~,hess]=fmincon(@map_costlearning,params,[],[],[],[],pmin,pmax,[],options);
                secondmoment = diag(inv(hess))';
                %dbstop in MLEM_sm.m at 210 if sum(secondmoment < 0)>0
                fitparams = [fitparams; subj z bestparam secondmoment map]; %for the iters over subject, which iters produced what
            end
            [bestmap,which] = min(fitparams(:,end));
            bestforsubj = [bestforsubj; subj fitparams(which,3:3+(nparams-1)) fitparams(which,3+nparams:end-1) bestmap];
            %track subj num, alpha, llh
        end
        if bestforsubj(:,3) == 0 & modeltofit.epsilon
            count = count + 1; %count, per model, how many times epsilon goes to 0
        end
        %re-infer mu and sigma?
        currentparams = bestforsubj(:,2:2+(nparams-1));
        currentsecondmoments = bestforsubj(:,(2+nparams):end-1);
        paramsbytry{tries} = currentparams;
        smbytry{tries} = currentsecondmoments;
        %newmu = mean(currentparams);
        %newsigma = std(currentparams);
        %meanLLH = mean(bestforsubj(:,end));
        mapsubjs = mean(bestforsubj(:,end)); %multiply all subj maps or add all llh's
        bestdist = [bestdist; mus(tries,:) sigmas(tries,:) mapsubjs];
    end %of cycling through random/best mu and sigma starting points
    [MAP,which] = min(bestdist(:,end));
    bestmu = mean(paramsbytry{which}); 
    bestsigma = sqrt(mean(paramsbytry{which}.^2 + smbytry{which})-(bestmu.^2)); %variance
    %dbstop in MLEM_sm.m at 161 if sum(~isreal(bestsigma))>0
    while sum(~isreal(bestsigma))>0 %ill-fitting?
        thrownout = [thrownout; MAP];
        bestdist(which,:) = []; %erase last minimum, poorly fitting 
        idxes = 1:tries; idxes(which) = [];
        [MAP,which] = min(bestdist(:,end)); %replace with second-best
        bestmu = mean(paramsbytry{idxes(which)});
        bestsigma = sqrt(mean(paramsbytry{idxes(which)}.^2 + smbytry{idxes(which)})-(bestmu.^2)); %variance
    end
    for subj = 1:nsubjs
        mu = bestmu; sigma = bestsigma;
        if istable(data)
            idx = data.subj == subj;
        else
            idx = data(:,1) == subj;
        end
        onesim = data(idx,:);
        params = bestforsubj(:,2:2+(nparams-1));
        bestfit(subj) = map_costlearning(params(subj,:));
    end
    modelscores = [modelscores; MAP];
    highlconvergence(iters,:) = [bestmu bestsigma];
end %of trying to find the best mu and sigma
% % MODEL COMPARISON SCORES % %
%calculate relative model fit by getting mean LLH with parameters sampled
%from best fitting distribution, with fit parameters
k = nsubjs; llhfromsample = []; %holding modeltofit constant
for sample = 1:k
    params = [];
    for i = 1:nparams
        params(i) = normrnd(bestmu(i),bestsigma(i));
    end
    llhfromsample(sample) = llh_costlearning(params);
end

[modeltofit.map,which] = min(modelscores);
modeltofit.lowparams = bestforsubj(:,2:2+(nparams-1));
modeltofit.highparams = [bestmu bestsigma];
modeltofit.fitbysubj = bestfit;
modeltofit.llh = mean(llhfromsample);
%save all of this somewhere useful
if min(modelscores)>min(thrownout)
    oop = true
end

end

