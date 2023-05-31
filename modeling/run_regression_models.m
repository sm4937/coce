function [outputArg1,outputArg2] = run_regression_models(X,Y,labels,wordy_flag)
%run_regression_models: Given an X, and a Y, run the most complex model,
% assess its significance, then eliminate predictors to see what is really
% driving the significance

% X is the matrix of predictors
% Y is a vector of the measure to be predicted
% labels is a list of which predictors are in X, in order
% wordy_flag is a boolean, which you toggle on/off depending how detailed
% you want the model output descriptions to be on the command line

%wordy_flag = true;

% this func does not assume these have been cleaned of NaN's, so just in
% case, cleans them of NaNs here 
% removing all rows containing NaNs across X and Y
invalid = sum(isnan(X),2)>0; Y(invalid,:) = []; X(invalid,:) = [];
invalid = sum(isnan(Y),2)>0; X(invalid,:) = []; Y(invalid,:) = [];

% run regression with MATLAB 'regress' function
[betas,BINT,~,~,stats] = regress(Y,X);
p_fullmodel = stats(3);
p_medmodel = 1;
p_reducedmodel = 1;

% calculates mean squared error between regression predicted outputs and
% true Y values
predicted = X*betas;
distance = predicted-Y;
MSE_fullmodel = distance'*distance; %squared distance
MSE_reducedmodel = MSE_fullmodel;


display_labels = labels;
for b = 1:length(display_labels)
    display_labels{b} = [display_labels{b} ' '];
end

if p_fullmodel < 0.05
    disp(' ')
    disp(['full model: significance p = ' num2str(p_fullmodel)])
    disp([display_labels{1:end}])
    disp(['MSE = ' num2str(MSE_fullmodel)])
    for b = 1:length(betas)
        disp(['beta ' display_labels{b} ' = ' num2str(betas(b))])
    end

    % significant full model
    % what's driving that?
    % test other model variants by chopping model up
    
    intercept_column = find(contains(labels,'intercept'));
    X_reduced = X(:,intercept_column);
    % keep intercept term (usually 1st in matrix)
    
    % grab the predictors (or models) that has betas which
    % don't include 0 in confidence interval
    % (they have the same sign, both positive or both negative)
    tokeep = find(sign(BINT(:,1))==sign(BINT(:,2)));
    tokeep(tokeep==intercept_column) = [];
    predictors = [intercept_column];
    
    if isempty(tokeep)
        % the intercept is currently the only predictive term,
        % cycle over the others to see if tacking them on adds anything
        
        for pred = 2:length(labels)
            
            X_reduced = [X(:,1) X(:,pred)];
            [betas,BINT,~,~,stats] = regress(Y,X_reduced);
            p_reducedmodel = stats(3);
            
            if p_reducedmodel < 0.05
                
                % calculates mean squared error between regression predicted outputs and
                % true Y values
                predicted = X_reduced*betas;
                distance = predicted-Y;
                MSE_reducedmodel = distance'*distance; %squared distance
                
                disp(' ')
                disp(['Reduced model p = ' num2str(p_reducedmodel)])
                disp(['Model includes intercept and ' display_labels{pred}])
                disp(['Reduced model MSE = ' num2str(MSE_reducedmodel)])
                disp(['Beta ' num2str(betas(2))])
                disp(['Beta interval ' num2str(BINT(2,:))])
                predictors = [predictors pred];
            end

        end % of cycling over possible predictors
        
    elseif length(predictors) > 1 || ~isempty(tokeep)
        % are there MULTIPLE good predictors other than the intercept term?
        % then, make medium mixed model
        
        if length(predictors) == 1
            predictors = [predictors tokeep'];
        end

        X_medium_complexity = X(:,predictors);
        [betas,BINT,~,~,stats] = regress(Y,X_medium_complexity);
        p_medmodel = stats(3);
        predicted = X_medium_complexity*betas;
        distance = predicted-Y;
        MSE_medmodel = distance'*distance; %squared distance

        if (p_medmodel < 0.05)
            disp(' ')
            disp('Medium model wins!')
            disp(['MSE = ' num2str(MSE_medmodel)])
            disp(['p = ' num2str(p_medmodel)])
            disp(['Significant predictors are ' display_labels{predictors}])
            disp(['With betas of ' num2str(betas')])
            disp(['Confidence intervals at ' ])
            for bii = 1:size(BINT,1)
                disp([num2str(BINT(bii,:)) ','])
            end
            disp(' ')
        end
        
    end %of assessing "good" predictors & how they combine
    
else % full model NOT significant
    
    disp('No significant predictors of this variable.')
    
end % of assessing the significance in the full model

if wordy_flag
    % plot significant relationships
    if p_medmodel < 0.05 | p_reducedmodel < 0.05
        figure
        scatter(X(:,predictors(2)),Y)
        lsline
        xlabel(display_labels{predictors(2)})
        ylabel('Y')
    end
end


end

