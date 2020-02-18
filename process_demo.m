function [demo_vars] = process_demo(raw_data)
%process_demo Pull out basic demographic information, like sex, age, and
%education levels
%   Pull out questionnaire responses

%rows = 386:387; %hard-coded for testing but will be pull-out-able later with:
rows = find(raw_data.task == categorical({'demographics'}));
rows = [rows rows+1]; %two relevant rows
%questions = {'age','education','sex','handedness'};
edulevels = ["Some elementary or middle school","Some high school","Graduated high school","Some college","Associate's degree","Bachelor's degree","Started higher education (MA, PhD, MD, etc)","Completed higher education (MA, PhD, MD, etc.)"];
sexes = ["Male","Female"];
hands = ["Right-handed","Left-handed","Ambidextrous"];

responses = string(raw_data.responses(rows));
responses = [responses string(table2array(raw_data(rows,end-3:end-1)))]; %this may need to change in the future
line1 = split(responses(1),'"');
line2 = split(responses(2),'"');

ageidx = find(contains(line2,'age'))+2;
educationidx = find(ismember(line1,'education'))+2; %this will all need to change
sexidx = find(ismember(line1,'sex'))+2; 
%handidx = find(ismember(line1,'handedness'))+2;

age = str2num(char(line2(ageidx)));
edu = line1(educationidx);
edulevel = find(edu==edulevels);
sex = find(line1(sexidx)==sexes);
%handedness = find(line1(handidx)==hands);
if isempty(sex)
    sex = NaN;
end
if isempty(age)
    age = NaN;
end
if isempty(edu)
    edu = NaN;
end

demo_vars.age = age;
demo_vars.edu = edu;
demo_vars.sex = sex;
%demo_vars.handedness = handedness;

end

