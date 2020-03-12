%% statistics
% build model of BDM value by
% 1. accuracy on last run-through of that task
% 2. n-switches on last run-through of that task (too correlated with 1?)
% 3. NFC
% 4. SAPS

long = [];
for subj = 1:n
    hard = find(data.task_progression(subj,:)==tasks(4));
    next = hard+1; 
    %easy = find(data.task_progression(subj,:)==tasks(3));
    hard(next==default_length+1) = []; next(next==default_length+1) = []; %trim last trials, no final BDM to measure response to perf
    short = [data.perf(subj,hard)' data.nswitches(subj,hard)' repmat(data.NFC(subj),length(hard),1) repmat(data.SAPS(subj),length(hard),1) data.allnswitchcosts{subj,:}(:,hard)' data.values(subj,next)'];
    long = [long; short];
end

inter = [long(:,1:4) long(:,2).*long(:,3) long(:,5)];

fitlm(long(:,1:5),long(:,6)) %linear regression of 4 predictors, predicting BDM value

%fitlm(inter(:,1:5),inter(:,6)) %linear regression of 5 predictors, including an interaction term, predicting BDM value

%% This is fine, but let's run a random effects model of this, too, where subject is the random effect
% model of BDM value by
% 1. accuracy on last run-through of that task
% 2. n-switches on last run-through of that task (too correlated with 1?)
% 3. NFC
% 4. SAPS
% with categorical variables of 5. subject
% and 6. session

long = [];
for subj = 1:n
    hard = data.task_progression(subj,:)==tasks(4); hard(inattentive(subj,:)) = false; %remove trials where subjects were not paying attention
    hard = find(hard); %turn into workable index
    next = hard+1; 
    hard(next==default_length+1) = []; next(next==default_length+1) = []; %trim last trials, no final BDM to measure response to perf
    short = [data.perf(subj,hard)' data.nswitches(subj,hard)' repmat(data.NFC(subj),length(hard),1) repmat(data.SAPS(subj),length(hard),1) repmat(subj,length(hard),1) repmat(data.session(subj),length(hard),1) data.values(subj,next)'];
    long = [long; short];
end

tbl = table;
tbl.performance = long(:,1); 
tbl.nswitches = long(:,2); %-min(long(:,2)); %or 0 baseline them?
%-mean(long(:,2)); %negative correlation of intercept and slope, need to
%center the n-switches column to prevent ill-fitting model?
%tbl.NFC = long(:,3); tbl.SAPS = long(:,4); % to include these or just to relate them to measures from mixed effects model?
tbl.NFC = long(:,3);
tbl.subj = categorical(long(:,5)); %tbl.sess = categorical(long(:,6)); 
% For now, excluding session as a random variable.
% "If you only have two or three levels, the model will struggle to partition the variance - it will give you an output, but not necessarily one you can trust." only 4 sessions for now
tbl.y = long(:,7);
%fitlme uses program A as a reference and creates the necessary dummy variables I[.]. Since the model already has an intercept, fitlme only creates dummy variables for programs B, C, and D. 
%This is also known as the 'reference' method of coding dummy variables. 

no_random = fitlme(tbl,'y ~ 1 + nswitches');

subj_random_intercept = fitlme(tbl,'y ~ 1 + nswitches + ( 1 | subj)'); %random slope of subj, interaction of nswitches x subj
% best fitting model accoridng to AIC and BIC, and less complex
[B,~,rEffects] = randomEffects(subj_random_intercept);
% best fitting model by AIC, BIC (489, 503), and less complex than the next
% model down
[r,p] = corr(B(~isnan(data.NFC)),data.NFC(~isnan(data.NFC)),'Type','Spearman'); %no relationship of random intercept term and NFC score

subj_random_slope = fitlme(tbl,'y ~ nswitches + (-1 + nswitches | subj)');
[B, ~, rEffects] = randomEffects(subj_random_slope);

subj_intercept_slope = fitlme(tbl,'y ~ 1 + nswitches + (nswitches|subj)');
% random slope and random intercept for each subject
[~,~,rEffects] = randomEffects(subj_intercept_slope);
estimates_inter = fixedEffects(subj_intercept_slope);
% not enough data for this ?

tbl.slopes = NaN(height(tbl),1); tbl.intercepts = NaN(height(tbl),1);
for subj = 1:n
    slopes = rEffects.Estimate(2:2:n*2);
    intercepts = rEffects.Estimate(1:2:n*2);
    slope = slopes(subj)+estimates_inter(2);
    inter = intercepts(subj)+estimates_inter(1);
    indxs = tbl.subj == categorical({num2str(subj)});
    tbl.slopes(indxs) = slope;
    tbl.intercepts(indxs) = inter;
end %pick out slopes and intercepts for each subject - maybe plotting the wrong thing on the negative correlation graph?

% figure();
% plotResiduals(subj_random_slope,'fitted')

subj_indie_inter_slope = fitlme(tbl,'y ~ nswitches + (1 | subj) + (-1 + nswitches | subj)');

comparison = compare(subj_indie_inter_slope,subj_intercept_slope); 
p_inter_wins = comparison.pValue < 0.05;
if p_inter_wins
    disp('Interaction model wins')
else
    disp('Interaction model does not win')
end
%the interaction model contains the non-interacting model because the non-interacting model is just a parameter-specific case of the interacting one

% figure,
% scatter(rEffects.Estimate(1:2:end),rEffects.Estimate(2:2:end)) 
% title('Random Effects','FontSize',15)
% xlabel('Intercept','FontSize',15)
% ylabel('Slope','FontSize',15)

% The estimated column in the random-effects table shows the estimated best
% linear unbiased predictor (BLUP) vector of random effects.
% The problem is that the slope and intercepts are negatively correlated
% which is a sign of a bad fit p = 0.008, r = -0.54
% plus AIC and BIC are worse at 493 and 514, even though the model is more
% complex

[~,~,rEffects] = randomEffects(subj_intercept_slope);
estimates_inter = fixedEffects(subj_intercept_slope);
figure
subplot(2,3,1)
scatter(rEffects.Estimate(1:2:n*2)+estimates_inter(1),rEffects.Estimate(2:2:n*2)+estimates_inter(2))
title('Slope vs. Intercept in Interacting Slope/Inter Model')
[r,p] = corr(rEffects.Estimate(1:2:n*2),rEffects.Estimate(2:2:n*2));
if p<0.05
    lsline
end
ylabel('Slope')
xlabel('Intercept')

subplot(2,3,2)
scatter(1:n,rEffects.Estimate(1:2:n*2)+estimates_inter(1))
title('Fit Intercepts in Interacting Slope/Inter Model')
xlabel('Subject')
ylabel('Fit Intercepts')

subplot(2,3,3)
scatter(1:n,rEffects.Estimate(2:2:n*2)+estimates_inter(2))
title('Fit Slopes in Interacting Slope/Inter Model')
xlabel('Subject')
ylabel('Fit Slopes')

[~,~,rEffects] = randomEffects(subj_indie_inter_slope);
estimates_indie = fixedEffects(subj_indie_inter_slope);

subplot(2,3,4)
scatter(rEffects.Estimate(1:n)+estimates_indie(1),rEffects.Estimate(n+1:n*2)+estimates_indie(2))
[r,p] = corr(rEffects.Estimate(1:n),rEffects.Estimate(n+1:n*2));
if p<0.05
    lsline
end
title('Intercept v Slope in Independent Slope/Inter Model')
xlabel('Intercept')
ylabel('Slope')

subplot(2,3,5)
scatter(1:n,rEffects.Estimate(1:n)+estimates_indie(1))
title('Fit Intercepts in Independent Slope/Inter Model')
xlabel('Subject')
ylabel('Intercept')

subplot(2,3,6)
scatter(1:n,rEffects.Estimate(n+1:n*2)+estimates_indie(2))
title('Fit Slopes in Independent Slope/Inter Model')
xlabel('Subject')
ylabel('Slope')

AICs = [-no_random.LogLikelihood.*2+2.*0, ...
    -subj_random_intercept.LogLikelihood.*2+(2.*numel(subj_random_intercept.covarianceParameters{1})), ...
    -subj_random_slope.LogLikelihood.*2+(2.*numel(subj_random_slope.covarianceParameters{1})), ...
    -subj_intercept_slope.LogLikelihood.*2+(2.*numel(subj_intercept_slope.covarianceParameters{1})), ... 
    -subj_indie_inter_slope.LogLikelihood.*2+(2.*numel(subj_indie_inter_slope.covarianceParameters))];
figure
bar(AICs)
hold on
[score,which] = min(AICs);
plot(which,score,'*k','LineWidth',1.5)
title('AIC by model, fit to real data')
ylabel('AIC')
xlabel('Model')
xticklabels({'No Random Effects','Random Intercept','Random Slope','Random Slope and Intercept','Random Slope and Intercept, non-interacting'})
xtickangle(45)

figure
[~,~,rEffects] = randomEffects(subj_random_intercept);
estimates = fixedEffects(subj_random_intercept);
medi = median(data.NFC(~isnan(data.NFC)));
for row = 1:n
    if data.NFC(row)>medi
       color = 'or';
    else
       color = 'ob';
    end
    scatter(row,rEffects.Estimate(row)+estimates(1),'k','Filled')
    hold on
end
xlabel('Subject Number')
ylabel('Fit Intercept')
ax = gca; fig = gcf;
ax.FontSize = 12; fig.Color = 'w';
title('Winning model, subject by intercept value')

%% make fake data and try to fit it
% generate and recover
n_initializations = 5;
info = [];
slope_means = [3 2 7 1 9];
for j = 1:n_initializations
% random intercept model
n = 22; %number of subjects
ntrials = 15;
% 
% M = rand(n,1); %more or less random slopes around 0
% B = rand(n,1)*3; 
% X = ceil(rand(ntrials,1)*15);

X = repmat(randperm(15,15)',ntrials/15,1);
%X = randperm(ntrials,ntrials);
%X = zscore(X); % center for scaling?
means = [10 5];
M = slope_means(j)+randn(n,1);
B = means(2)+randn(n,1); %M-5; %

ys = []; xs = []; subjs = [];
for subj = 1:n
   %y = M(subj).*X(subj,1:ntrials) + B(subj); 
   X = reshape(X(randperm(length(X),length(X))),length(X),1); %shuffle, make column
   %re-initialize X every time
   noise = rand(length(X),1)*2;
   y = M(subj).*X + B(subj) + noise;
   ys = [ys; y]; 
   xs = [xs; X(1:ntrials)];
   subjs = [subjs; repmat(subj,ntrials,1)];
end

tbl = table;
tbl.x = xs;
tbl.y = ys;
tbl.subj = subjs;

no_random = fitlme(tbl,'y ~ 1 + x'); %no random effects

subj_random_intercept = fitlme(tbl,'y ~ 1 + x + ( 1 | subj)'); %random intercept of subj
% best fitting model according to AIC and BIC, and less complex
[B_hat,~,rEffects] = randomEffects(subj_random_intercept);

[r,p] = corr(B_hat,B); %can successfully recover intercepts this way
figure 
subplot(5,1,1)
scatter(B_hat,B)
xlabel('Fit Intercept')
ylabel('Real Intercept')
title('Random Intercept Model')

subj_random_slope = fitlme(tbl,'y ~ x + (-1 + x | subj)');
[M_hat, ~, rEffects] = randomEffects(subj_random_slope);
subplot(5,1,2)
%M_hat = rEffects.Estimate(1:2:end);
scatter(M_hat,M)
xlabel('Fit Slope')
ylabel('Real Slope')
title('Random Slope Model')

% figure();
% plotResiduals(subj_random_slope,'fitted')

n_repetitions = 1;
plot_flag = n_repetitions==1; %n_repetitions==5; %how many times do you want to repeat fitting with same data? should be the same every time
both_models = [];

clearvars rEffects
for loop = 1:n_repetitions

subj_intercept_slope = fitlme(tbl,'y ~ 1 + x + (x|subj)');
% random slope and random intercept for each subject
[~,~,rEffects] = randomEffects(subj_intercept_slope);
estimates_inter = fixedEffects(subj_intercept_slope);
if plot_flag
    subplot(3,3,1)
    scatter(rEffects.Estimate(2:2:end)+estimates_inter(2),M) %fix plots by adding mean estimate !!!!!! %%%%%
    xlabel('Fit Slope')
    ylabel('Real Slope')
    title('Random Intercept and Slope Model')
    
    subplot(3,3,2)
    scatter(rEffects.Estimate(1:2:end)+estimates_inter(1),B)
    xlabel('Fit Intercept')
    ylabel('Real Intercept')
    title('Random Intercept and Slope Model')    

    subplot(3,3,3)
    scatter(rEffects.Estimate(1:2:end)+estimates_inter(1),rEffects.Estimate(2:2:end)+estimates_inter(2))
    xlabel('Fitted Intercept')
    ylabel('Fitted Slope')
    title('Fitted Random Effects of Intercept and Slope Model')
end

subj_indie_inter_slope = fitlme(tbl,'y ~ x + (1 | subj) + (-1 + x | subj)');
[~,~,rEffects2] = randomEffects(subj_indie_inter_slope);
estimates_indie = fixedEffects(subj_indie_inter_slope);
if plot_flag
    subplot(3,3,5)
    scatter(rEffects2.Estimate(1:n)+estimates_indie(1),B)
    xlabel('Fit Intercept')
    ylabel('Real Intercept')
    title('Independent Slope and Intercept Model')

    subplot(3,3,4)
    scatter(rEffects2.Estimate(n+1:n*2)+estimates_indie(2),M) %slope is second value for each subject's rEffects
    xlabel('Fit Slope')
    ylabel('Real Slope')
    title('Independent Slope and Intercept Model')

    subplot(3,3,6)
    scatter(rEffects2.Estimate(n+1:n*2)+estimates_indie(2),rEffects2.Estimate(1:n)+estimates_indie(1)) %slope is second value for each subject's rEffects
    xlabel('Fit Slope')
    ylabel('Fit Intercept')
    title('Independent Slope and Intercept Model')
    
    subplot(3,3,7)
    scatter(rEffects2.Estimate(1:n),rEffects.Estimate(1:2:end))
    title('Fit Intercepts - do they agree?')
    xlabel('Non-interacting fit')
    ylabel('Interacting fit')
    
    subplot(3,3,8)
    scatter(rEffects2.Estimate(n+1:n*2),rEffects.Estimate(2:2:end))
    title('Fit Slopes - do they agree?')
    xlabel('Non-interacting fit')
    ylabel('Interacting fit')
end

AICs = [-no_random.LogLikelihood.*2+2.*0, ...
    -subj_random_intercept.LogLikelihood.*2 + (2.*numel(subj_random_intercept.covarianceParameters{1})), ...
    -subj_random_slope.LogLikelihood.*2 + (2.*numel(subj_random_slope.covarianceParameters{1})), ...
    -subj_intercept_slope.LogLikelihood.*2 + (2.*numel(subj_intercept_slope.covarianceParameters{1})), ...
    -subj_indie_inter_slope.LogLikelihood.*2 + (2.*numel(subj_indie_inter_slope.covarianceParameters{1}))];

both_models = [both_models; AICs(4:5)];

comparison = compare(subj_indie_inter_slope,subj_intercept_slope); 
p_inter_wins = comparison.pValue < 0.05;
if p_inter_wins
    disp('Interaction model wins')
else
    disp('Interaction model does not win')
end

%model comparison
if plot_flag
    figure
    bar(AICs)
    title('AIC by Model')
    ylabel('AIC')
    xlabel('Model')
    xticklabels({'No Random Effects','Random Intercept','Random Slope','Random Slope and Intercept','Random Slope and Intercept, non-interacting'})
    xtickangle(45)
end

end

if plot_flag %stop plotting this
    figure
    %AICs of the two models with slope and intercept terms 
    subplot(1,2,1)
    bar(nanmean(both_models,1))
    hold on
    errorbar(nanmean(both_models,1),nanstd(both_models,[],1)/sqrt(n_repetitions),'*','LineWidth',1)
    title([num2str(n_repetitions) ' runs of model fits for random intercept and slope models'])
    ylabel('Mean AIC')
    xlabel('Model')
    xticklabels({'Interacting','Non-Interacting'})
    
    subplot(1,2,2)
    plot(1:n_repetitions,both_models(:,1),'m','LineWidth',1.5)
    hold on
    plot(1:n_repetitions,both_models(:,2),'c','LineWidth',1.5)
    title('Fit by repetition')
    legend({'Interacting','Non-Interacting'})
end

info = [info; M B repmat(AICs(4:5),length(M),1) repmat(j,length(M),1)];

end

%save(['info_' num2str(j) 'reps_' num2str(n) 'subjects_' num2str(ntrials) 'trials_correlated.mat'],'info')

%% Plot results of this generate and recover
%load('info_51reps_22subjects_15trials_correlated.mat')
%load('info_51reps_22subjects_15trials_uncorrelated.mat')
inits = sum(info(:,5)==1);
n_initializations = length(unique(info(:,5)));
figure
bar([nanmean(info(:,3)) nanmean(info(:,4))])
hold on
errorbar([nanmean(info(:,3)) nanmean(info(:,4))],[nanstd(info(:,3)) nanstd(info(:,4))]/sqrt(n_initializations),'*','LineWidth',1)
diff = nanmean(info(:,3)-info(:,4));
title(['mean difference: ' num2str(diff)])
%ylim([-6000 -5000])
ylabel('Mean AIC')
xlabel('Model')
xticklabels({'Interacting','Non-Interacting'})

figure
subplot(1,2,1)
AICs = info(1:inits:end,3:4);
plot(1:n_initializations,AICs(:,1),'b','LineWidth',1.5)
hold on
plot(1:n_initializations,AICs(:,2),'r','LineWidth',1.5)
[~,bests(1)] = min(abs(AICs(:,1))); [~,bests(2)] = min(abs(AICs(:,2)));
plot(bests(1),AICs(bests(1),1),'*k')
plot(bests(2),AICs(bests(2),2),'*k')
legend({'Interacting','Non-Interacting'})
title('AICs for each model by repetition')

subplot(1,2,2)
values = [];
for init = 1:n_initializations
   [r,p] = corr(info(init:(init+inits-1),1),info(init:(init+inits-1),2));
   values = [values; r p];
end
plot(1:n_initializations,values(:,1),'k','LineWidth',1.5)
idx = values(:,2)<0.05;
hold on
plot(find(idx),values(idx,1))
title('Correlation of random slopes and intercepts')

figure
slopes = info(:,1);
AIC = info(:,3);
AIC2 = info(:,4);
matrix = sortrows([slopes AIC AIC2],1);
scatter(matrix(:,1), matrix(:,2)) 
hold on
scatter(matrix(:,1), matrix(:,3))
%check out correlation of M and AIC score- probably highly related since no noise added in generation procedure
xlabel('Slope')
ylabel('AIC score')

%% Simulate slopes based on NFC

alpha = 0.50; %baseline slope
NFC = 1; %fake, as if NFC were between 0 and 1
B = -10:10;
S = alpha*(1./(1+exp(-B.*NFC)));
figure
scatter(B,S)
ylabel('Slope')
title(['Formula alpha*sigmoid(-B*NFC), NFC: ' num2str(NFC)])
xlabel('Beta Value')
xticks(B)
xlim([B(1) B(end)])
xticklabels(B)
ylim([-1 1])

NFCs = data.NFC;
maxNFC = 5;
NFCs = NFCs./5;

Beta = 3;
alpha = 3;
slopes = alpha*(1-(1./(1+exp(-Beta.*NFCs))));
%slopes = alpha + (Beta.*NFCs);
nswitches = 6:15;
Bs = rand(n,1)*3;
y = slopes.*nswitches + Bs;

figure
subplot(2,1,2)
toplot = 5:10;
for subj = toplot
    scatter(nswitches,y(subj,:),'o','Filled')
    hold on
end
xlabel('N Switches')
ylabel('BDM value')
title('Simulated effect of nswitches according to NFC')
labels = [repmat('NFC: ',length(toplot),1) num2str(NFCs(toplot))];
legend({labels})
ax = gca; 
ax.FontSize = 12;

subplot(2,1,1)
scatter(slopes,NFCs)
xlabel('Slope')
ylabel('NFC score normalized')
title('Slope of subcomponent effect by NFC score')
ax = gca; fig = gcf;
fig.Color = 'w'; ax.FontSize = 12;

