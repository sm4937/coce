function [datawithcounts] = countEntries(data)
%countEntries Counts the number of unique entries in a matrix, where each
%entry is a row and each column contains distinct information
%   Adds an extra column which indexes how many of that entry there are in
%   the list so far
%   Helpful when trying to plot scatter plots with differing size dots,
%   which are meant to display how many data points are there in that
%   cluster by the size of the dot.

w = size(data,2);
h = size(data,1);

for i = 1:h %go through each entry
    row = data(i,:);
    count = sum(sum(data==row,2)==w); %number of exact matches
    datawithcounts(i,:) = [row count];
end

if isempty(data)
    datawithcounts = [NaN NaN NaN];
end

end

