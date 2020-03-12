function [normalized_NFC,normalized_SAPS] = grade_scales(raw_data)
%gradeNFC Grade the need for cognition scale from raw text data
%   Takes raw answers to need for cognition scale and grades according to 
% The Efficient Assessment of need for cognition (Cacioppo, Petty, 1984)
screener_answers = {'extremely_characteristic_of_me','somewhat_agree'};
reverse = [3, 4, 5, 7, 8, 9, 12, 16, 17];

%rows = 383:385; %hard-coded for testing but will be pull-out-able later with:
responses = string(raw_data.responses(end));
separated = split(responses,'___');
screeners = find(contains(separated,'screener'));
totalNFC = 0; answered = 0;
all_NFC_answers = find(contains(separated,'characteristic')|contains(separated,'uncertain'));
if ~isempty(screeners) %some subjects don't have screeners in their questionnaires, from early
    screenidx = find(ismember(all_NFC_answers,screeners(1)+1)); %trim screener from NFC responses
    passed_screeners(1) = screener_answers{1} == separated(all_NFC_answers(screenidx));
    all_NFC_answers(screenidx)=[]; %don't count screener
else
    passed_screeners = 1;
    if length(all_NFC_answers)~=18
        disp(['subj num ' num2str(raw_data.subjnum(1)) ': Number of questions wrong for NFC, check out'])
    end
end

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
    totalNFC = totalNFC+grade;
    answered = answered+(grade~=0);
end

normalized_NFC = totalNFC/answered;%control for # questions answered

%% grade SAPS now, I think it's just a sum??

totalSAPS = 0; answered = 0;
all_SAPS_answers = find(contains(separated,'agree')|contains(separated,'neutral'));
if ~isempty(screeners)
    screenidx = find(ismember(all_SAPS_answers,screeners(2)+1)); %trim screener from SAPS responses
    passed_screeners(2) = screener_answers{2} == separated(all_SAPS_answers(screenidx));
    all_SAPS_answers(screenidx)=[]; %don't count screener
else %no screeners included
    passed_screeners = 1;
    if length(all_NFC_answers)~=8
        disp('Number of questions wrong for SAPS, check out')
    end
end

responses = {'strongly_disagree', 'disagree', 'somewhat_disagree', 'neutral', 'somewhat_agree', 'agree', 'strongly_agree'};

for q = 1:length(all_SAPS_answers) %all questions, of 18        
    response = separated(all_SAPS_answers(q));
    grade = find(ismember(responses,response));
    if isempty(grade)
        grade = 0;
    end
    totalSAPS = totalSAPS+grade;
    answered = answered+(grade~=0);
end

normalized_SAPS = totalSAPS/answered; %control for # questions answered

if sum(passed_screeners)==0
    normalized_SAPS = NaN;
    normalized_NFC = NaN;
    totalSAPS = NaN;
    totalNFC = NaN;
end

end

