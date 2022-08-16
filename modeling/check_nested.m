function [nestedflag] = checkNested(model1name,model2name)
%checkNested Check whether one model is nested in another
%   Because they're named in coc_createModels you can just compare the
%   names to see whether one contains the other 
%   Works best for one-param differences
%   Checks whether model 1 contains model 2

nestedflag = contains(model1name,model2name);

end

