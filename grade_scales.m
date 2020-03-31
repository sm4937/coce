function [normalized_NFC,normalized_SAPS] = grade_scales(raw_data)
%gradeNFC Grade the need for cognition scale from raw text data
%   Takes raw answers to need for cognition scale and grades according to 
% The Efficient Assessment of need for cognition (Cacioppo, Petty, 1984)
screener_answers = {'extremelycharacteristicofme','somewhatagree'};

%rows = 383:385; %hard-coded for testing but will be pull-out-able later with:
responses = string(raw_data.responses(end));
separated = split(responses,'___');
corr_flag = 1; %mark with a flag whether the data printed correctly or not
if length(separated)<=1 %the columns didn't get condensed properly in jspsych
    responsecol = find(ismember(raw_data.Properties.VariableNames,'responses'));
    responsecols = [responsecol find(contains(raw_data.Properties.VariableNames,'VarName'))];
    acrosscols = raw_data(end-6:end,responsecols);
    acrosscols = [string(acrosscols.responses) table2array(acrosscols(:,2:end))]; %deal with disparate variable types
    for input = 1:numel(acrosscols)
        row = acrosscols(input);
        new_input = strrep(row,' ','');
        new_input = strrep(new_input,'_','');
        new_input = strrep(new_input,':','');
        cleanup(input) = new_input;
    end
    separated = cleanup;
    corr_flag = 0;
    normalized_NFC = NaN; normalized_SAPS = NaN; %for now, leave scores out from people whose data didn't save properly
    return
end

% clean up separated a little
for item = 1:length(separated)
    response = separated(item);
    response2 = split(response,'_'); %split by _ character to erase it
    response3 = '';
    for word = 1:length(response2)
        response3 = response3 + response2(word); %new string, mashed words of response together
    end
    separated(item) = response3; %re-assign
end

%% Grade screeners
screeners = find(contains(separated,'screener'));
all_SAPS_answers = find(contains(separated,'agree')|contains(separated,'neutral'));
all_NFC_answers = find(contains(separated,'characteristic')|contains(separated,'uncertain'));
if ~isempty(screeners) %exclude screeners from each questionnaire
    for i = 1:length(screeners)
        screenidx = find(ismember(all_NFC_answers,screeners(i)+(corr_flag))); %trim screener from SAPS responses
        if ~isempty(screenidx)
            passed_screeners(i) = screener_answers{1} == separated(all_NFC_answers(screenidx));
            all_NFC_answers(screenidx)=[]; %don't count screener
        end
        screenidx = find(ismember(all_SAPS_answers,screeners(i)+(corr_flag)));
        if ~isempty(screenidx)
            passed_screeners(i) = screener_answers{2}==separated(all_SAPS_answers(screenidx));
            all_SAPS_answers(screenidx)=[]; %don't count screener
        end
    end
else
    passed_screeners = 1;
    if length(all_NFC_answers)~=18
        disp(['subj num ' num2str(raw_data.subjnum(1)) ': Number of questions wrong for NFC, check out'])
    end
end

%% Grade NFC
responses = {'extremelyuncharacteristicofme','somewhatuncharacteristicofme','uncertain','somewhatcharacteristicofme','extremelycharacteristicofme'};
question_labels = separated(all_NFC_answers-1); 
reverse = [3, 4, 5, 7, 8, 9, 12, 16, 17];

for l = 1:length(question_labels)
    label = question_labels(l);
    label = strsplit(label,'C'); 
    if length(label)>1
        label = str2num(char(label(2)));
    else
        if length(all_NFC_answers) == 18 %if the number of questions is right,
            label = l; %label according to order
        else %if it's not, no NFC score
            label = NaN;
            break
        end
    end  
    question_labels_2(l) = label;
end

totalNFC = 0; answered = 0;
for q = 1:length(all_NFC_answers) %all questions, of 18        
    response = separated(all_NFC_answers(q));
    grade = find(contains(responses,response));
    label = question_labels_2(q);
    if ismember(label,reverse)
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

responses = {'stronglydisagree', 'disagree', 'somewhatdisagree', 'neutral', 'somewhatagree', 'agree', 'stronglyagree'};

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

