%% Analysis of cost of control task data 
% meant to be posted on MTurk
% task written in JsPsych, data printed in CSV form
%% Process raw jspsych text file data
clear all
subjnum = 23893644;
filename = '12.02.2020_trial';

cutoff = 80; 
long_format = import_mturk_data(['./data/' filename '.csv']);
single = long_format(long_format.subjnum==subjnum,:);
data = make_data_table(single); 
task_progression = data.task_progression;
TOT = data.TOT;
perf_by_block = data.perf;
values = data.values;
BDM_rt = data.BDMrt/1000;

%% Plot some stuff
subjectID = [' ' num2str(subjnum) ' '];

if data.version == 1
    hardacc = data.nbackacc;
    easyacc = data.detectacc;
    tasklabels = {'detection','combine'};
    easyrts = data.detectrts;
    hardrts =  data.nbackrts;
    tasks = [categorical({'detection'});categorical({'combine'})];
    hardswitchcosts = data.combineswitchcosts{1};
elseif data.version == 2
    hardacc = data.nswitch9acc;
    easyacc = data.nswitch1acc;
    tasklabels = {'p(switch)=0.1','p(switch)=0.9'};
    easyrts = data.nswitch1rts;
    hardrts = data.nswitch9rts;
    tasks = [categorical({'pswitch0.1'}), categorical({'pswitch0.9'})];
    hardswitchcosts = data.allnswitchcosts{1};
end

figure
subplot(1,2,1)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
bar([nanmean(easyacc) nanmean(hardacc)])
hold on
errorbar(1:2,[nanmean(easyacc) nanmean(hardacc)],[nanstd(easyacc)/length(easyacc) nanstd(hardacc)/length(hardacc)],'k*','LineWidth',2)
xticklabels(tasklabels)
ylim([50 100])
title(['Subj' subjectID 'accuracy by task'])

subplot(1,2,2)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
bar([nanmean(easyrts) nanmean(hardrts)])
hold on
errorbar(1:2,[nanmean(easyrts) nanmean(hardrts)],[nanstd(easyrts)/length(easyrts) nanstd(hardrts)/length(hardrts)],'k*','LineWidth',2)
xticklabels(tasklabels)
title(['Subj' subjectID 'mean RT by task'])

figure
subplot(1,2,1)
ax = gca; fig = gcf;
passingnback = (perf_by_block>=cutoff&task_progression==tasks(2)); passingnback = [false passingnback(1:end-1)]; %offset by 1
passingdetect = (perf_by_block>=cutoff&task_progression==tasks(1)); passingdetect = [false passingdetect(1:end-1)]; %offset by 1
bluevalues = values(passingnback);
purplevalues = values(passingdetect);
redvalues = values(~passingnback&~passingdetect);
plot(find(passingnback),bluevalues,'o','MarkerFaceColor','b')
hold on
plot(find(passingdetect),purplevalues,'o','MarkerFaceColor','m')
plot(find(~passingnback&~passingdetect),redvalues,'o','MarkerFaceColor','r')
fig.Color = 'w';
ax.FontSize = 12;
title(['Subj' subjectID 'points by BDM round'])
xlabel('Block #')
ylabel('Points requested')
legend({'Passed hard','Passed easy','Did not Pass'},'Location','Best')

subplot(1,2,2)
ax = gca; fig = gcf;
plot(hardacc,'o')
lsline
title(['Subj' subjectID 'acc in hard task'])
fig.Color = 'w';
ax.FontSize = 12;
xlabel('Block #')
ylabel('Accuracy')

figure
bluevalues = BDM_rt(passingnback);
purplevalues = BDM_rt(passingdetect);
redvalues = BDM_rt(~passingnback&~passingdetect);
plot(find(passingnback),bluevalues,'o','MarkerFaceColor','b')
hold on
plot(find(~passingnback&~passingdetect),redvalues,'o','MarkerFaceColor','r')
plot(find(passingdetect),purplevalues,'o','MarkerFaceColor','m')
xlabel('Round')
ylabel('Reaction time (seconds)')
title(['Subj' subjectID 'RT by BDM round'])
legend({'Passed hard','Did not Pass','Passed easy'},'Location','Best')
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;

figure
next_trial = find(task_progression==tasks(2))+1;
y = hardacc(~isnan(hardacc))';
y(next_trial==length(task_progression)+1) = []; %prune lists the same way
next_trial(next_trial == length(task_progression)+1) = []; %outside the index, 13th block of 12
matrix = sortrows([y,values(next_trial)'],1);
plot(matrix(:,1),matrix(:,2),'o')
lsline
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
xlim([50 100])
title(['Subj' subjectID 'BDM value by prev. accuracy on pictured task'])
ylabel('BDM points')
xlabel('Accuracy')

%plot last BDM offer by subject request on next trial
%subject 10 tried to game this
figure
y = data.offers';
x = [NaN data.values(2:end)]';
matrix = sortrows([y,x],1);
plot(matrix(:,1),matrix(:,2),'o')
lsline
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title(['Subj' subjectID 'BDM value by prev. computer offer'])
ylabel('BDM points')
xlabel('Last offer')

%plot points at stake vs. following performance
%subject 9 felt more motivated by higher offers
figure
init = NaN(12,1);
init(data.offers>data.values) = data.offers(data.offers>data.values);
init(data.values>data.offers) = 1;
y = init;
x = [NaN perf_by_block(1:end-1)]';
matrix = sortrows([y,x],1);
plot(matrix(:,1),matrix(:,2),'o')
lsline
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;
title(['Subj' subjectID 'perf by points at stake'])
xlabel('Points')
ylabel('Accuracy')

figure
y = TOT(task_progression==tasks(1));
x = TOT(task_progression==tasks(2));
scatter(ones(length(y),1),y)
hold on
scatter(2*ones(length(x),1),x)
errorbar(1:2,[mean(y) mean(x)],[std(y)/length(y) std(x)/length(x)],'ko')
xlim([0.75 2.25])
ylabel('Time on task (seconds)')
xticks([1 2])
xticklabels(tasklabels)
ax = gca; fig = gcf;
fig.Color = 'w';
ax.FontSize = 12;

figure
% plot switch costs in hard switch task, easy switch task
plot(ones(sum(task_progression==tasks(1)),1),hardswitchcosts(task_progression==tasks(1)),'k*')
hold on
plot(2*ones(sum(task_progression==tasks(2)),1),hardswitchcosts(task_progression==tasks(2)),'k*')
ax = gca; fig = gcf;
title('switch costs by task')
xticks([1 2])
xlim([0.5 2.5])
xticklabels(tasklabels);
xtickangle(45);

