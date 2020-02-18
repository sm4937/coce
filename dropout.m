%%Analyze dropouts
%subject by subject, print out a new data table for failed people
function [failed_data] = dropout(raw_data)
    failed_data = table;
    subjnum = unique(raw_data.subjnum);
    quiz_long = unique(double(string(raw_data.BDMquizgrade)));
    quiz = quiz_long(1);
    passed_quiz = quiz>67;
    passed_practice = ~isnan(quiz);
    passed_consents = height(raw_data)>2;
    tasknums = double(string(raw_data.tasknum));
    nums = unique(tasknums(~isnan(tasknums)));
    blocknums = max(nums);
    passed_halfway_mark = blocknums>11;
    all_the_way = blocknums>22;
    if ~passed_quiz
        passed_halfway_mark = false;
        all_the_way = false;
    end
    failed_data.values = [passed_consents passed_practice passed_quiz passed_halfway_mark all_the_way];
    failed_data.labels = {'consents','practice','quiz','halfway','finished'};
    failed_data.exp_version = unique(raw_data.exp_version);
end
