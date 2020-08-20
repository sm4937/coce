function [groups] = tertileSplit(measurement)
%tertileSplit Takes in a measurement and returns tertile split groups
%   groups is the output, with 1 for lowest group membership, 2 for
%   second lowest, 3 for highest. Splits by bins.


matrix = [(1:length(measurement))' measurement NaN(length(measurement),1)];
%takes in column vector and appends "subject numbers"
%blank final variable for group membership

purescores = unique(measurement);
purescores = purescores(~isnan(purescores));
thirds = [ceil(length(purescores)/3) ceil(2*length(purescores)/3) length(purescores)]; %indices
lowbound = purescores(thirds(1)); highbound = purescores(thirds(2));
idx1 = matrix(:,2)<lowbound; idx2 = matrix(:,2)>=lowbound & matrix(:,2)<highbound; idx3 = matrix(:,2)>=highbound;
matrix(idx1,3) = 1; matrix(idx2,3) = 2; matrix(idx3,3) = 3;

groups = matrix(:,3);

%sortedm = sortrows(matrix,2,'ascend'); %sort by measurement
%thirds = [1 length(matrix)/3 2*length(matrix)/3 length(matrix)]; %indices
% sortedm(thirds(1):thirds(2),3) = 1;
% sortedm(thirds(2):thirds(3),3) = 2;
% sortedm(thirds(3):thirds(4),3) = 3;
% sortedm = sortrows(sortedm,1,'ascend'); %go back to first order
% groups = sortedm(:,3);

end

