// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function setup_detection_n(practice_flag,loopi,n){ //called n-back but it's both that and the combine

  // initialize empty lists for variable storage

  taskID = 7; 

  var stamp = '<div style="position: absolute; left: 5px; top: 5px">';
  stamp += fractals[taskID]+' width="'+stamp_size+'" height ="'+stamp_size+'"';
  stamp += '</div>';

  var prob_oddball = 0.20; //0.20; //kool botvinick 2014 overlapping detection task
  var task_label = 'ndetection';

  var response_color = 'red'; //black to red upon response in n-back, detection

  // set up initial global variables
  var answer_key_names = [universal_key,'H'];
  var answer_keys = [jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[0]),jsPsych.pluginAPI.convertKeyCharacterToKeyCode('z'),jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[1])];

  var n = 2; //set how many back they're looking for
  //this should just get fed in at setup_detection_n, but for some reason it's not registering all the time
  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        //var n_back_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + String(n) + ' letters ago.</p> <p>If the image on screen was not the one you saw ' + String(n) + ' letters ago, do not press anything.</p><p>If you see a T, press <strong>'+ answer_key_names[1] + '</strong> as fast as you can.</p><p>Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this task.</p> <p>Press space to begin.</p>';
        var instructs_text = '<p>In this task, you will press a button every time you see the same letter ' + String(n+1) + ' times in a row (we call this a match).</p><p>When you see a match, you need to press <strong>' + answer_key_names[0] + '</strong>.</p><p>If the letter on screen is not the same one you saw ' + String(n-1) + ' and ' + String(n) + ' letters ago, do not press anything.</p><p>There may be multiple targets in a row.</p><p>Press any key to begin.</p>'; 
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return instructs_text;
      },
      trial_duration: instructs_timing //30 seconds to respond
    }]
  };

  var presentation_screen = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: '+',
      trial_duration: 100
    }]
  };

  //length of each block
  ntrials = 15;

  //nmatches = 3;
  rand = Math.random();
  if(rand<1){ //some blocks with more matches, randomly, approx 1/3 divided
    nmatches = 2;
  }
  if(rand<0.66){
    nmatches = 3;
  }
  if(rand<0.33){
    nmatches = 4;
  }

  if(practice_flag){
    ntrials = 10;
    nmatches = 2;
    var presentation_screen = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: '<p> Task: </p> ' + fractals[taskID]+' width="'+fractal_size+'" height ="'+fractal_size+'"',
        prompt: '<p>This picture will always be associated with the following task, like a picture label.</p><p>You will now learn the rules of this task and get a chance to practice.</p><p>Press any key to continue.</p>'
      }]
    };
  }

  var notyet = true;
  while(notyet){
    var stimnum_list = []; var data_list = []; var key_list = []; var detect_list = [];
    var trial_type = new Array(ntrials-n);
    trial_type.fill(0);
    trial_type.fill(1,0,nmatches);
    //trial_type = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; //add 2's in to make this the combine
    var temp = new Array(n);
    temp.fill(0);
    var trial_list = temp.concat(randomizeList(trial_type)); //after n trials, begin random assignment of trial identity
    for(var trial = 0; trial < ntrials; trial++){
      stimnum = Math.floor(Math.random()*num_stim);
      //make it flexible
      stimolds = [stimnum_list[trial-n], stimnum_list[trial-1]]; //this is fixed for now but it shouldn't be
      if(trial_list[trial]==1){ //force it to be an n-back trial
        stimnum = stimolds[0];
        stimnum_list[trial-1] = stimolds[0]; //reference first one in the 3 stimuli, make those retroactively the same
        answer_key = answer_keys[0]; //this is written specifically for n = 2
      }
      if(trial_list[trial]==0){
        answer_key = answer_keys[1]; //placeholder key
      }
      // n-backs may happen organically, so make sure this is accounted for
      if((stimnum==stimolds[0]&stimnum==stimolds[1])==1){
        if(trial_list[trial]==1){
          answer_key = answer_keys[0];
        }
        if(trial_list[trial]==0){
          answer_key = answer_keys[1]; //placeholder key
          var stimnums = range(0,num_stim-1,1)
          var idx = stimnums.indexOf(stimolds[0]); //remove this stimulus from the list of possible stimuli
          stimnums.splice(idx,1);
          rand = Math.floor(Math.random()*(stimnums.length));
          stimnum = stimnums[rand];
        } // change the stimnum so it's not a match, thereby standardizing the number of matches each round
      }
      stimnum_list.push(stimnum);
      detect_list.push(trial_list[trial]==1);
      data_list.push(false);
      key_list.push(answer_key);
    }; //end of task building

    //catch trials where the backwards stimulus assignment procedure creates extra matches, re-counterbalance
    //if nmatches == nactualmatches, then move forward
    var matchcount = 0;
    for(var trial = 0; trial < ntrials; trial++){
      stims = [stimnum_list[trial-n], stimnum_list[trial-1],stimnum_list[trial]]; 
      if((stims[2]==stims[0]&stims[2]==stims[1])==1){
        matchcount+=1;
      }
    }//end of matchcounting

    if(matchcount==nmatches){
      notyet = false;
    } //discrepancy between goal matches and actual matches

  }; // end of while loop, move forward in task creation

  var test_stimuli = [];
  for(var trial = 0; trial < ntrials; trial++){
    test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], response_color: response_color, correct_key: key_list[trial], data: {nback: data_list[trial], stimnum: stimnum_list[trial], correct_key: key_list[trial], detect: detect_list[trial], n: n, nmatches: nmatches, task: task_label, tasknum: loopi}})
  };

  var detection = {
      timeline: create_color_change_timeline(test_stimuli,"",stamp)
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  var debrief = {
    type: "html-keyboard-response",
    stimulus: function() {
     overall = accuracyNdetect(); cutoff_message = pointsFeedback(); performance_list.push(overall);
     return "<p style='font-size:25px'>" + cutoff_message + "</p> </p><p style='font-size:25px'>Press any key to continue. </p>";
    },
    trial_duration: instructs_timing, //30 seconds to respond
    data: function(){
      dict = {task: "debrief", perf: accuracyNdetect()}
      return dict;
    }
  };

  if(practice_flag){ //write practice task n = 10 trials
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
          detect = detect[0]==1;
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
            if(detect){ //n-back match but missed
              return "<p>You missed a match!</p><p>Press <strong>" + answer_key_names[0] + "</strong> when you see a match.</p>"
            }
            if(!detect){ //not an n-back trial
              return "Not a match! Press nothing."
            }
            if(!answer_keys.includes(response)&response!=null){ // didn't press a valid key
              return "Please press " + answer_key_names[0] + "."
            }
          } 
        },
        stimulus_duration: feedback_timing,
        trial_duration: feedback_timing,
        response_ends_trial: false, //allow no response, just feedback
        data: {task: 'feedback', practice: true}
      }
      
      test.push(stimulus_screen);
      test.push(feedback_screen);
      // create your own dynamic feedback
    }; // end of test set up loop for trials 0-9

    var detection = {
      timeline: test
    } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

    var debrief = {
      type: "html-keyboard-response",
      data: function(){
        accuracy = accuracyNdetect();
        return {practice_accuracy: accuracy, practice: true, task: 'debrief'}
      },
      stimulus: function(){
        accuracy = accuracyNdetect(); 
        //return "<p>Your accuracy was <strong>" + accuracy + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate during the study to earn points.</p><p>Please remember that you will not receive feedback like this during the real experiment.</p><p>Press the space bar to continue.</p>";
        return "<p>Your accuracy was <strong>" + accuracy + "%</strong>.</p><p>Please remember that you will not receive feedback like this during the real experiment.</p><p>Press any key to continue.</p>";    
      }
    };

  } // end of practice flag

  timeline = [presentation_screen,instructs,detection,debrief];
  return timeline;

}; // end of set up function 

function accuracyNdetect(data){ //calculate accuracy for Nback
  var correct_num = 0;
  var multiplier = 3; var detect_num = 0;
  lasttrialdata = jsPsych.data.getLastTimelineData().filter({task: 'ndetection'});
  var buttons = lasttrialdata.select('key_press').values; 
 //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
  var trialtype = lasttrialdata.select('detect').values;
  ntrials = buttons.length;
  for(var trial = 0; trial < ntrials; trial++){
    if(trialtype[trial]==1){//yes it's an n-back trial
      detect_num += 1;
      if(buttons[trial] == answer_keys[0]){
       correct_num+=(1*multiplier);
     }
    }
    if(trialtype[trial]==0){ //not an n-back trial
      if(buttons[trial]==null){
       correct_num+=1}
      }
    }
  overall = Math.round((correct_num/(ntrials+(multiplier-1)*detect_num))*100);
  overall = Number.parseFloat(overall).toFixed(2);
  return overall;
}
