function [scores] = gradeDifficultyRatings(raw_data)
%gradeDifficultyRatings Pull out subject difficulty ratings for each task
%   Take raw data, find relevant column and convert to 1-5 rating subjects
%   gave

ratings = raw_data.response;
ratings = ratings(raw_data.task==categorical(cellstr('difficultyrating')));
scores = [1+(ratings/25)]'; 

if isempty(ratings)
    scores = NaN(1,4);
end

end

