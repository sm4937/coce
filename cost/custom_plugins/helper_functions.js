// build functions to use across tasks

function pointsFeedback(){
  lasttrial = jsPsych.data.getLastTimelineData();
  tasknum = lasttrial.select("tasknum").values;
  tasknum = tasknum[0];
  points = points_list[tasknum];
  if(overall < cutoff){
    var cutoff_message = 'You were not accurate enough. You earned 0 points.';
    } else {
    var cutoff_message = 'Great job! You earned ' + String(points) + ' point(s) from this task.';
    };
  return cutoff_message;
};  

function generateOffer(){
  var computer = Math.floor(Math.random()*100);
  computer = 1 + computer*0.04;
  computer = Number.parseFloat(computer).toFixed(2);
  return computer;
}

function runBDM(points,offer,task_ID){
  var BDM_message = "<p style='font-size':25px> You want <strong> " + String(points) + " </strong> points. </p> <p> The computer offers <strong>" + String(offer) + "</strong> points. </p>"
  timeoutflag = points==null;
  if(points <= offer){
    BDM_message += "<p> You will be given " + String(offer) + " points if you achieve at least " + String(cutoff_percent)  + "% accuracy on this task. </p>" + fractals[task_ID] + " location='center'></img>";
    points = offer; 
  } else {
    BDM_message += "<p> You will do the following task instead to earn 1 point. </p> " + fractals[task_ID+1] + " location='center'></img>";
    points = 1; 
  }
  if(timeoutflag){
    BDM_message = "<p> You did not respond in time. You will do the following task to earn 1 point. </p>" + fractals[task_ID+1] + "location='center'></img>";
    points = 1;
  }
  BDM_message += "<p>[Press the space bar to continue.]</p>"
  points_list.push(points); offer_list.push(offer);
  data = {BDM_message: BDM_message, points_list: points_list, offer_list: offer_list}
  return data;
};

function getBDMresponse(){
  offer = generateOffer();
  offer = Number.parseFloat(offer).toFixed(2)
  var lasttrialdata = jsPsych.data.getLastTrialData();
  points = lasttrialdata.select("response").values;
  points = 1+points*0.04;
  points = Number.parseFloat(points).toFixed(2);
  rt = lasttrialdata.select('rt').values;
  if(rt[0]==null){
    points = null;
  }
  value_list.push(points);
  data = {offer: offer, value_list: value_list, points: points};
  return data
  };

function turnLogicalIntoAccuracy(input){
  count = 0;
  for(index = 0; index < input.length; index++){
    if(input[index] == true){
      count += 1;
    }
  }
  accuracy = count/input.length;
  return accuracy;
};

function keepTrackofPerformance(){
  sum = 0;
  for(i = 0; i < performance_list.length; i++){
    sum += parseInt(performance_list[i]);
  }
  overall = sum/performance_list.length;
  return overall;
};

function makeHTML(input){
  return "<p>" + String(input) + "</p>"
}

function sumPoints(input){
  sum = 0;
  for(i = 0; i < input.length; i++){
    point = Math.ceil(input[i]);
    if(performance_list[i]>=cutoff){
      sum += point;
    }
  }
  return sum;
}

var range = function(start, end, step) {
    var range = [];
    var typeofStart = typeof start;
    var typeofEnd = typeof end;
    if (step === 0) {
        throw TypeError("Step cannot be zero.");
    }
    if (typeofStart == "undefined" || typeofEnd == "undefined") {
        throw TypeError("Must pass start and end arguments.");
    } else if (typeofStart != typeofEnd) {
        throw TypeError("Start and end arguments must be of same type.");
    }
    typeof step == "undefined" && (step = 1);
    if (end < start) {
        step = -step;
    }
    if (typeofStart == "number") {
        while (step > 0 ? end >= start : end <= start) {
            range.push(start);
            start += step;
        }
    } else if (typeofStart == "string") {
        if (start.length != 1 || end.length != 1) {
            throw TypeError("Only strings with one character are supported.");
        }
        start = start.charCodeAt(0);
        end = end.charCodeAt(0);
        while (step > 0 ? end >= start : end <= start) {
            range.push(String.fromCharCode(start));
            start += step;
        }
    } else {
        throw TypeError("Only string and number types are supported");
    }
    return range;
}

function create_color_change_timeline(test_stimuli,key_presses){
  test = []; // empty test variable for practice, build more dynamic feedback by hijacking the plugins
    for(var i = 0; i < test_stimuli.length; i++){
      var stimulus_screen = {
        type: "categorize-html",
        stimulus: '<p style="color: ' + test_stimuli[i]['color'] + '" >' + test_stimuli[i]['stimulus'] + '</p>',
        data: test_stimuli[i]['data'],
        key_answer: test_stimuli[i]['correct_key'],
        prompt: key_presses,
        stimulus_duration: stim_timing,
        trial_duration: stim_timing,
        show_feedback_on_timeout: false,
        timeout_message: ' ',
        feedback_duration: 0,
        correct_text: ' ',
        show_stim_with_feedback: false,
        incorrect_text: ' '
      }
      //console.log(stimulus_screen)

      var feedback_screen = {
        type: "html-keyboard-response",
        stimulus: '<p style="color: ' + test_stimuli[i]['response_color'] + '">' + test_stimuli[i]['stimulus'] + '</p>',
        trial_duration: function(){
          grab = jsPsych.data.getLastTrialData();
          RT = grab.select('rt').values;
          if(RT[0]!=null){
            duration = stim_timing-RT[0];
          }else{
            duration = 0;
          }
          return duration
        },
        //stimulus_duration: feedback_timing,
        response_ends_trial: false, //allow no response, just feedback
        data: {task: 'feedback'}
      }

      var ITI = {
        type: "html-keyboard-response",
        stimulus: ' ',
        trial_duration: ITI_timing,
        response_ends_trial: false,
        data: {task: 'ITI'}
      }

      test.push(stimulus_screen);
      test.push(feedback_screen);
      test.push(ITI)
      // create your own dynamic feedback
    }; // end of test set up loop for trials 0-9
return test
}

function compareToCutoff(accuracy){
  if(!debug){
    if(accuracy<cutoff){
      return false; //cut off
    }
    if(accuracy>=cutoff){
      return true; //go right ahead
    }
    if(isNaN(accuracy)){
      return false; //past thing was skipped, skip rest of it
    }
  }else{
    return true; //debug mode, just do everything
  }
}

function run_task_by_ID(task_ID,loopi){
  //task list 
  //0 is combine
  //1 is detection
  //2 is n-back
  //3 is n-switch hard
  //4 is n-switch easy
  if(task_ID==0){
    var temp = setup_nBack(loopi);
    return temp
  }
  if(task_ID==1){
    var temp = setup_detection(false,loopi);
    return temp
  }
  if(task_ID==2){
    var temp = setup_practice_nback(3);
    return temp
  }
  if(task_ID==3){
    var temp = setup_nswitch(false,loopi,1);
    return temp
  }
  if(task_ID==4){
    var temp = setup_nswitch(false,loopi,2); //easy
    return temp
  }
}

function randomizeList(list){
  var oldlist = list;
  var newlist = [];
  var nend = oldlist.length;
  for(var i = 0; i < nend; i++){
    rand = Math.floor(Math.random()*oldlist.length);
    //console.log(rand)
    value = oldlist[rand];
    //console.log(idx)
    newlist.push(value);
    //console.log(trial_type)
    oldlist.splice(rand,1);//randomly index into indices, use that to reference trial_type, then delete the value
  }
  return newlist
}

function getColNames(){
  //sarah-written
  var xhr = new XMLHttpRequest();
  xhr.open('GET', 'php/get_data_needed_from_db.php');
  xhr.responseType = 'text';
  xhr.send();
  xhr.onload = function() {
    if(xhr.status == 200){
      //var response = JSON.parse(xhr.responseText);
      var response = xhr.responseText;
      console.log(response);
      console.log('status 200');
    }
  };
  col_names = readline();
  return col_names;
}

function saveDataToCsv(name, data){
  var xhr = new XMLHttpRequest();
  xhr.open('POST', 'php/write_data_to_csv.php'); // 'write_data.php' is the path to the php file described above.
  xhr.setRequestHeader('Content-Type', 'application/json');
  xhr.send(JSON.stringify({filename: name, filedata: data}));
}

function setDataColumns(){
  //this will be automatic some day, for now it's manual
  turkInfo = jsPsych.turk.turkInfo();
  turk = !turkInfo.outsideTurk;
  if(turk){
    var urlVar = jsPsych.data.urlVariables();
    var worker = urlVar.workerId;
    var assignment = urlVar.assignmentId;
    var hitid = urlVar.hitId;
  }
  if(!turk){
    var worker = null;
    var assignment = null;
    var hitid = null;
  }
  dict = {
    worker_ID: worker,
    assignment_ID: assignment,
    hid_ID: hitid,
    //task: "instructions",
    subjnum: subjnum*1,
    exp_version: exp_version
  }
  jsPsych.data.addProperties(dict)
  /*dict = {rt: null,
    responses:null,
    question_order: null, 
    trial_type: null, 
    trial_index: null, 
    time_elapsed: null, 
    internal_node_id: null,
    TOT: null,
    overall: null,
    performance_by_block: null,
    points_list: null,
    value_list: null,
    offer_list: null,
    exp_version: null,
    n_switch_rule_sides: null,
    subjnum: subjnum,
    worker_ID: worker,
    assignment_ID: assignment,
    hit_ID: hitid,
    success: null,
    stimulus: null,
    key_press: null,
    task: 'instructions',
    correct: null,
    detect: null,
    correct_key: null,
    tasknum: null,
    practice_accuracy: null,
    practice: null,
    nback: null,
    n: null,
    number_practice_hard: null,
    response: null,
    BQMquizgrade: null,
    task_displayed: null,
    slider_start: null,
    stimnum: null,
    interaction_data: null
  }*/ //end data dict, defining all fields from the beginning
  return dict
}