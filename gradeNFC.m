function [total,normalized] = gradeNFC(raw_data)
%gradeNFC Grade the need for cognition scale from raw text data
%   Takes raw answers to need for cognition scale and grades according to 
% The Efficient Assessment of need for cognition (Cacioppo, Petty, 1984)

reverse = [3, 4, 5, 7, 8, 9, 12, 16, 17];
responses = {'extremely uncharacteristic','somewhat uncharacteristic','uncertain','somewhat characteristic','extremely characteristic'};

%rows = 383:385; %hard-coded for testing but will be pull-out-able later with:
rows = raw_data.task == categorical({'NFC'});
indx1 = find(contains(raw_data.Properties.VariableNames,'responses'));
NFC_pages = raw_data(rows,indx1:end);
expression = '\Q';
all_answers = [];
total = 0; answered = 0;
for i=1:length(NFC_pages)
    one_page = regexp(NFC_pages(i),expression,'split');
    page = one_page(2:end);
    all_answers = [all_answers page];
end

for q = 1:length(all_answers) %all questions, of 18        
    response = all_answers(q);
    for z = 1:length(responses)
        match = contains(response,responses(z));
        if match %kind of a convoluted way to do this, lazy about removing noise from questionnaire answers
            grade = z;
            if ismember(q,reverse)
                grade = 6 - grade;
            end
        else
            grade = 0;
        end
        total = total+grade;
        answered = answered+match;
        if match
            break
        end
    end
end

normalized = total/answered;%control for # questions answered
end

