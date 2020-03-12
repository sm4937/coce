%%Analyze dropouts
%subject by subject, print out a new data table for failed people
function [failed_data] = dropout(raw_data)
    default_length = 10;
    failed_data = table;
    subjnum = unique(raw_data.subjnum);
    quiz_long = unique(double(string(raw_data.BDMquizgrade)));
    quiz = quiz_long(1);
    passed_quiz = quiz>67;
    passed_practice = ~isnan(quiz);
    passed_consents = height(raw_data)>4;
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
    
    if passed_consents
        % add in practice performance for people who made it past practice
        practices = double(string(raw_data.number_practice_hard));
        acc = double(string(raw_data.practice_accuracy));
        accuracy = acc(raw_data.task==categorical({'debrief'})); 
        practice_accuracy = NaN(1,default_length); 
        if max(practices) > 0
            accuracy = accuracy(end-max(practices):end);
            practice_accuracy(:,1:length(accuracy)) = accuracy;
        end
        failed_data.number_practices = max(practices)+1;
        failed_data.practice_accuracy = practice_accuracy;
    else
        failed_data.number_practices = 0;
        failed_data.practice_accuracy = NaN(1,default_length);
    end
    
end
