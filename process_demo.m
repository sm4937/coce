function [demo_vars] = process_demo(raw_data)
%process_demo Pull out basic demographic information, like sex, age, and
%education levels
%   Pull out questionnaire responses

rows = 386:387; %hard-coded for testing but will be pull-out-able later with:
%raw_data.task == categorical({'demographics'});
%questions = {'age','education','sex','handedness'};
edulevels = ["Some elementary or middle school","Some high school","Graduated high school","Some college","Associate's degree","Bachelor's degree","Started higher education (MA, PhD, MD, etc)","Completed higher education (MA, PhD, MD, etc.)"];
sexes = ["Male","Female"];
hands = ["Right-handed","Left-handed","Ambidextrous"];

responses = raw_data.responses(rows);
line1 = split(responses(1),'"');
line2 = split(responses(2),'"');

ageidx = find(ismember(line2,'age'))+2;
educationidx = find(ismember(line1,'education'))+2;
sexidx = find(ismember(line1,'sex'))+2;
handidx = find(ismember(line1,'handedness'))+2;

age = str2num(char(line2(ageidx)));
edu = line1(educationidx);
edulevel = find(edu==edulevels);
sex = find(line1(sexidx)==sexes);
handedness = find(line1(handidx)==hands);

demo_vars.age = age;
demo_vars.edu = edu;
demo_vars.sex = sex;
demo_vars.handedness = handedness;

end

