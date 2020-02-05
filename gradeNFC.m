function [score] = gradeNFC(raw_data)
%gradeNFC Grade the need for cognition scale from raw text data
%   Takes raw answers to need for cognition scale and grades according to 
% The Efficient Assessment of need for cognition (Cacioppo, Petty, 1984)

reverse = [3, 4, 5, 7, 8, 9, 12, 16, 17];
responses = {'extremely uncharacteristic','somewhat uncharacteristic','uncertain','somewhat characteristic','extremely characteristic'};

rows = 382:384; %hard-coded for testing but will be pull-out-able later with:
%raw_data.task == categorical({'NFC'});

NFC_pages = raw_data.responses(rows);


end

