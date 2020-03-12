function [total,normalized] = grade_scales(raw_data)
%gradeNFC Grade the need for cognition scale from raw text data
%   Takes raw answers to need for cognition scale and grades according to 
% The Efficient Assessment of need for cognition (Cacioppo, Petty, 1984)
screener_answers = {'extremely_characteristic_of_me','somewhat_agree'};
reverse = [3, 4, 5, 7, 8, 9, 12, 16, 17];

%rows = 383:385; %hard-coded for testing but will be pull-out-able later with:
responses = string(raw_data.responses(end));
separated = split(responses,'___');
screeners = find(contains(separated,'screener'));
total = 0; answered = 0;
all_NFC_answers = find(contains(separated,'characteristic')|contains(separated,'uncertain'));
screenidx = find(ismember(all_NFC_answers,screeners(1)+1)); %trim screener from NFC responses
passed_screeners(1) = screener_answers{1} == separated(all_NFC_answers(screenidx));
all_NFC_answers(screenidx)=[]; %don't count screener

responses = {'extremely_uncharacteristic_of_me','somewhat_uncharacteristic_of_me','uncertain','somewhat_characteristic_of_me','extremely_characteristic_of_me'};

for q = 1:length(all_NFC_answers) %all questions, of 18        
    response = separated(all_NFC_answers(q));
    grade = find(ismember(responses,response));
    if ismember(q,reverse)
       grade = 6 - grade;
    end
    if isempty(grade)
        grade = 0;
    end
    total = total+grade;
    answered = answered+(grade~=0);
end

normalized = total/answered;%control for # questions answered

end

