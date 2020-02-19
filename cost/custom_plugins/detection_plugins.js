// functions and such just for detection creation 
// set up the task parameters and stimulus list
function setup_detection(practice,loopi){
  var p_target = 0.2; //set probability of target stimulus
  detection_target = "T";

  var stimulus_list = []; var stimnum_list = []; 
  var data_list = []; var key_list = []; //initialize empty lists for task parameters

  var response_color = 'red'; //black to red upon response in n-back, detection
  
  if(!practice){ //regular task
    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function() {
          var initial_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[1] + '</strong> every time you see the target letter <strong>' + detection_target + '</strong>. </p><p>You must respond while the letter is on screen or your answer will not count.</p><p> Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this task.</p> <p>Press space to begin.</p>';
           //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
          return initial_instructs;
        },
        data: {task: 'instructions'},
        trial_duration: instructs_timing //30 seconds to respond
      }]
    };

    var presentation_screen = {
      type: "html-keyboard-response",
      stimulus: '+',
      trial_duration: 100
    };

    var debrief = {
      type: "html-keyboard-response",
      stimulus: function() {
        overall = accuracyDetection(); cutoff_message = pointsFeedback(); performance_list.push(overall);
        return "<p>" + cutoff_message + "</p> </p><p>Press the space bar to continue. </p>";
      },
      data: function(){
        dict = {task: 'debrief', perf: accuracyDetection()};
        return dict;
      }, 
      trial_duration: instructs_timing //30 seconds to respond
    };

  ntrials = 15;

  trial_type = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  indices = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14];
  trial_list = [];

  for(var trial = 0; trial < ntrials; trial++){
    rand = Math.floor(Math.random()*indices.length);
    //console.log(rand)
    idx = indices[rand];
    //console.log(idx)
    trial_list.push(trial_type[idx]);
    //console.log(trial_type)
    indices.splice(rand,1);//randomly index into indices, use that to reference trial_type, then delete the value
  }

  for(var trial = 0; trial < ntrials; trial++){
    stimnum = Math.floor(Math.random()*num_stim);
    //make it flexible
    if(trial_list[trial]==1){ //force it to be a detection trial
      stimnum = num_stim;
      answer_key = answer_keys[2];
    }
    if(trial_list[trial]==0){
      answer_key = answer_keys[1]; //placeholder key
    }
    stimnum_list.push(stimnum);
    data_list.push(stimnum==num_stim);
    key_list.push(answer_key);
  }; //end of task building

  var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
      test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], response_color: response_color, correct_key: key_list[trial], data: {detect: data_list[trial], nback: false, n: 0, stimnum: stimnum_list[trial], task: 'detection', correct_key: key_list[trial], tasknum: loopi}})
    }; 

  var detect = {
    timeline: create_color_change_timeline(test_stimuli,"")
  }

  } else { // make a practice task

    var presentation_screen = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: fractals[1],
        prompt: '<p>Practice task.</p><p>This picture will always be associated with the following task, like a picture label.</p><p>You will now learn the rules of this task and get a chance to practice.</p><p>Press the space bar to continue.</p>'
          }]
      };

    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
          stimulus: function() {
            //var initial_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[1] + '</strong> every time you see the target letter <strong>' + detection_target + '</strong>. </p><p>You must respond while the letter is on screen or your answer will not count.</p><p> Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this round.</p> <p>Press space to begin.</p>';
            var initial_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[1] + '</strong> every time you see the target letter <strong>' + detection_target + '</strong>. </p><p>You must respond while the letter is on screen or your answer will not count.</p><p>Press space to begin.</p>';            
             //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
            return initial_instructs;
              },
            data: {task: 'instructions'}
            }]
          };    

      var debrief = {
       type: "html-keyboard-response",
       stimulus: function() {
         overall = accuracyDetection();
         return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate to progress to the main experiment.</p><p>In the main experiment, you will not receive the same feedback as during practice.</p><p>Press the space bar to continue. </p>";
        },
        data: function(){
          accuracy = accuracyDetection();
          return {practice_accuracy: accuracy, practice: true, task: 'debrief'}
        }
      }

    ntrials = 10;

    trial_type = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0];
    indices = [0,1,2,3,4,5,6,7,8,9];
    trial_list = [];

    for(var trial = 0; trial < ntrials; trial++){
      rand = Math.floor(Math.random()*indices.length);
      //console.log(rand)
      idx = indices[rand];
      //console.log(idx)
      trial_list.push(trial_type[idx]);
      //console.log(trial_type)
      indices.splice(rand,1);//randomly index into indices, use that to reference trial_type, then delete the value
    }

    for(var trial = 0; trial < ntrials; trial++){
      stimnum = Math.floor(Math.random()*num_stim);
      //make it flexible
      if(trial_list[trial]==1){ //force it to be a detection trial
        stimnum = num_stim;
        answer_key = answer_keys[2];
      }
      if(trial_list[trial]==0){
        answer_key = answer_keys[1]; //placeholder key
      }
      stimnum_list.push(stimnum);
      data_list.push(stimnum==num_stim);
      key_list.push(answer_key);
    }; //end of task building

    var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
        test_stimuli.push({stimulus: makeHTML(stimuli[stimnum_list[trial]]), correct_key: key_list[trial], data: {detect: data_list[trial], correct_key: key_list[trial], task: 'detection', tasknum: -1}})
    };

    test = []; // empty test variable for practice, build more dynamic feedback by hijacking the plugins
    for(var i = 0; i < test_stimuli.length; i++){
      var stimulus_screen = {
        type: "categorize-html",
        stimulus: test_stimuli[i]['stimulus'],
        data: test_stimuli[i]['data'],
        key_answer: test_stimuli[i]['correct_key'],
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
        stimulus: function(){
          correct = jsPsych.data.getLastTrialData().select("correct").values; 
          correct = correct[0];
          response = jsPsych.data.getLastTrialData().select("key_press").values;
          response = response[0]; // remove list container
          answer = jsPsych.data.getLastTrialData().select("correct_key").values;
          answer = answer[0];
          detect = jsPsych.data.getLastTrialData().select("detect").values;
          detect = detect[0];
          if(correct){
            return "Correct!"
          }else{
            if(answer==answer_keys[1]&response==null){
              return " "
            } //not responding is automatically incorrect in categorize-html plugin but is right in this case
            //categorize-html is supposed to take care of this but okay
            if(response==answer){
              return "Correct!"
            }
            if(!answer_keys.includes(response)){ // didn't press a valid key
              return "Please press " + answer_key_names[1] + "."
            }
            if(detect){
              return "Press <strong>" + answer_key_names[1] + "</strong> when you see a T."
            }
            if(!detect){ //not an n-back trial
              return "Not a match! Press nothing."
            }
          } 
        },
        stimulus_duration: feedback_timing,
        trial_duration: feedback_timing,
        response_ends_trial: false, //allow no response, just feedback
        data: {task: 'feedback'}
      }
      test.push(stimulus_screen);
      test.push(feedback_screen);
      // create your own dynamic feedback
    }; // end of test set up loop for trials 0-9

    var detect = {
      timeline: test
    } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  }; //end of practice task creation

  timeline = [presentation_screen,instructs,detect,debrief];
  return timeline
};

function accuracyDetection(data){ //calculate accuracy for detection task (0-back)
  var correct_num = 0; var detect_num = 0;
  lasttrialdata = jsPsych.data.getLastTimelineData().filter({task: "detection"});
  detect = lasttrialdata.select("detect").values;
  var buttons = lasttrialdata.select("key_press").values; 
  ntrials = buttons.length;
  //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
  for(var trial = 0; trial < ntrials; trial++){
    if(detect[trial]){//yes it's a detection
      detect_num += 1;
      if(buttons[trial] == answer_keys[2]){
       correct_num+=4
      }
    } else { //not a detection trial
        if(buttons[trial]==null){
         correct_num+=1}
      }
    }
  ntrials = ntrials + detect_num*3; //weighted average (i.e. detection trials count for more)
  overall = (correct_num/ntrials)*100; 
  overall = Number.parseFloat(overall).toFixed(2);
  return overall;
 };