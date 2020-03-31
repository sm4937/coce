function [demo_vars] = process_demo(raw_data)
%process_demo Pull out basic demographic information, like sex, age, and
%education levels
%   Pull out questionnaire responses

%questions = {'age','education','sex','handedness'};
edulevels = ["Some elementary or middle school","Some high school","Graduated high school","Some college","Associate's degree","Bachelor's degree","Started higher education (MA, PhD, MD, etc)","Completed higher education (MA, PhD, MD, etc.)"];
sexes = ["Male","Female"];
hands = ["Right-handed","Left-handed","Ambidextrous"];

responses = string(raw_data.responses(end));
separated = split(responses,'___');
if responses=='NULL' %no demo data for this person
    demo_vars.age = NaN;
    demo_vars.edu = NaN;
    demo_vars.sex = NaN;
    demo_vars.hand = NaN;
    return 
end

if length(separated)>1
    ageidx = find(contains(separated,'age'))+1;
    sexidx = find(contains(separated,'male','IgnoreCase',true));
    handidx = find(contains(separated,'handed'));
    age = str2num(char(separated(ageidx(1))));
    edu = separated(sexidx-2);
    sex = find(separated(sexidx)==sexes);
else %the columns didn't get condensed properly in jspsych
    responsecol = find(ismember(raw_data.Properties.VariableNames,'responses'));
    responsecols = [responsecol find(contains(raw_data.Properties.VariableNames,'VarName'))];
    acrosscols = raw_data(end-8:end,responsecols);
    acrosscols = [string(acrosscols.responses) table2array(acrosscols(:,2:end))]; %deal with disparate variable types
    
    ageidx = find(contains(acrosscols,'age')); ageidx = ageidx(1); %second one is the 'age' in 'fair wage'
    sexidx = find(contains(acrosscols,'male','IgnoreCase',true));
    handidx = find(contains(acrosscols,'handed'));
    edu = NaN; %sexidx-1;
    
    age = acrosscols(ageidx); age = split(age,'":"'); number = split(age(2),'"'); age = str2num(char(number(1)));
    sex = split(acrosscols(sexidx),':"'); sex = char(sex(2)); sex = split(sex,'"'); sex = sex(1);
    sex = find(sexes==sex);
    separated = acrosscols;
end

handedness = separated(handidx);

expression = '\_handed';
splitStr = regexp(handedness,expression,'split');
if ~isempty(splitStr)
    hand = splitStr(1);
else
    if contains(separated,'Ambidextrous')
        hand = 'Ambidextrous';
    else
        hand = NaN;
    end
end

if isempty(sex)
    sex = NaN;
end
if isempty(age)
    age = NaN;
end
if isempty(edu)
    edu = NaN;
end
if isempty(hand)
    hand = NaN;
end

demo_vars.age = age;
demo_vars.edu = edu;
demo_vars.sex = sex;
demo_vars.hand = hand;

end

