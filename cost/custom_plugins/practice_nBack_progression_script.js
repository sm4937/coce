// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function setup_practice_nBack(n){
  var stimulus_list = []; var stimnum_list = []; var feedback_list = [];
  var data_list = []; var key_list = []; var detect_list = [];
  // initialize empty lists for variable storage

  var prob_n_back = 0.20; //kool botvinick 2014

   var presentation_screen = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: fractals[2],
        prompt: '<p>Practice game.</p><p>This picture will always be associated with the following task, like a picture label.</p><p>You will now learn the rules of this task and get a chance to practice.</p><p>Press the space bar to continue.</p>'
      }]
    };

  var back_rule = '<strong>' + String(n) + '</strong> letters ago';
  if(n==1){
    back_rule = '<strong>' + String(n) + '</strong> letters ago';
  }

  var back_diagram = "./static/img/"+String(n)+"-back.png";
  all_images.push(back_diagram)

  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function(){
        var n_back_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + back_rule + '.</p> <p>If the letter on screen was not the one you saw ' + back_rule + ', do not press anything.</p><p>Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate by the last practice game, you will not progress to the real experiment.</p>'
        n_back_instructs += '<img src='+back_diagram+'>';
        n_back_instructs += '<p>Press space to begin.</p>';
        return n_back_instructs;
      }
    }]
  };
  
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
    stimold = stimnum_list[trial-n]; //keep track of n trials back, starting with trial 0
    if(trial_list[trial]==1&trial>=n){ //force it to be an n-back trial
      stimnum = stimold;
      answer_key = answer_keys[0];
    }
    if(trial_list[trial]==0 || trial<n){
      answer_key = answer_keys[1]; //placeholder key
    }
    // n-backs may happen organically, so make sure this is accounted for
    if(stimnum==stimold){
      answer_key = answer_keys[0];
    }
    stimnum_list.push(stimnum);
    data_list.push(stimnum==stimold);
    key_list.push(answer_key);
  }; //end of task building

    //word_instructs = instructs;
    //instructs = [presentation_screen, instructs]; 
  var test_stimuli = [];
  for(var trial = 0; trial < ntrials; trial++){
     test_stimuli.push({stimulus: makeHTML(stimuli[stimnum_list[trial]]), correct_key: key_list[trial], data: {nback: data_list[trial], detect: detect_list[trial], correct_key: key_list[trial], task: 'n-back', tasknum: -1, n: n, practice: true}})
  };      
  var debrief = {
    type: "html-keyboard-response",
    data: function(){
      accuracy = accuracyNback();
      return {practice_accuracy: accuracy, practice: true}
    },
    stimulus: function(){
      overall = accuracyNback(); 
      return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate during the study to earn points.</p><p>Please remember that you will not receive feedback like this during the real experiment.</p><p>Press the space bar to continue.</p>";
    }
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
        nback = jsPsych.data.getLastTrialData().select("nback").values;
        nback = nback[0]==1;
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
            return "Please press " + answer_key_names[0] + "."
          }
          if(nback){ //n-back match but missed
            return "Press <strong>" + answer_key_names[0] + "</strong> when you see a match."
          }
          if(!nback){ //not an n-back trial
            return "Not a match! Press nothing."
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

  var nback = {
    timeline: test
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  timeline = [presentation_screen, instructs, nback, debrief];
  return timeline;

};//end of set up function


// create the practice detection/n-back task
function setup_practice_combo(){
  var stimulus_list = []; var stimnum_list = []; 
  var data_list = []; var key_list = []; var feedback_list = []; var detect_list = [];
  // initialize empty lists for variable storage

  var prob_n_back = 0.20; //kool botvinick 2014
  var prob_oddball = 0.20;

  var n = 3; //set how many back they're looking for
  // it's a 3-back in Kool Botvinick 2014

  var presentation_screen = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: fractals[0],
      prompt: '<p>Practice game.</p><p>This picture will always be associated with this task, like a picture label.</p><p>Press the space bar to continue.</p>'
    }]
  };

  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        var n_back_instructs = '<p>This task is a combination of the two tasks you just practiced.</p><p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + String(n) + ' letters ago (we call this a "3-back match").</p> <p>If the letter on screen was not the one you saw ' + String(n) + ' trials ago (not a "3-back match"), do not press anything.</p><p>If you see a T, press <strong>'+ answer_key_names[1] + '</strong> as fast as you can.</p><p>T\'s count in the sequence of letters for a 3-back match.</p><p>A T might be a 3-back match, but ignore that. Just press ' + answer_key_names[1] + ' when you see a T.</p><p>You will now be given six chances to practice this task until you are at least ' + String(cutoff_percent) + '% accurate on it.</p><p>If you don\'t reach ' +String(cutoff_percent)+ '% accuracy in six attempts, you will not move on to the real experiment.</p><p>Press space to begin.</p>';
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return n_back_instructs;
      }
    }]
  };
  
  //length of each block
  ntrials = 15;

  trial_type = [1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
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
    if(trial_list[trial]==2){
      stimnum = num_stim;
      answer_key = answer_keys[2]; //detection task
    }
    //make it flexible
    stimold = stimnum_list[trial-n]; //keep track of n trials back, starting with trial 0
    if(trial_list[trial]==1&trial>=n){ //force it to be an n-back trial
      stimnum = stimold;
      answer_key = answer_keys[0];
    }
    if(trial_list[trial]==0 || trial<n){
      answer_key = answer_keys[1]; //placeholder key
    }
    // n-backs may happen organically, so make sure this is accounted for
    if(stimnum==stimold){
      answer_key = answer_keys[0];
    }
    stimnum_list.push(stimnum);
    detect_list.push(stimnum==num_stim);
    data_list.push((stimnum==stimold)&(stimnum!=num_stim));
    key_list.push(answer_key);
  }; //end of task building

  var test_stimuli = [];
  for(var trial = 0; trial < ntrials; trial++){
       test_stimuli.push({stimulus: makeHTML(stimuli[stimnum_list[trial]]), correct_key: key_list[trial], fb: feedback_list[trial], data: {nback: data_list[trial], practice: true, correct_key: key_list[trial], detect: detect_list[trial], n: 3, task: 'combine', tasknum: -1}})
  };
  //console.log(test_stimuli)
  //console.log(test_stimuli[0]['stimulus'])

  var debrief = {
    timeline: [{
    type: "html-keyboard-response",
    data: function(){
      accuracy = accuracyCombine();
      return {practice_accuracy: accuracy, practice: true}
    },
    stimulus: function() {
      overall = accuracyCombine();
      if(overall>=cutoff){
        return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate during the study to earn points.</p><p>Please remember that you will not receive feedback like this during the real experiment.</p><p>Press the space bar to continue.</p>";
      } else {
        return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate during the study to earn points.</p><p>Please remember that you will not receive feedback like this during the real experiment.</p><p>Press the space bar to practice this again.</p>";
      } // end of n == 3 conditional
    }
  }]
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
        nback = jsPsych.data.getLastTrialData().select("nback").values;
        nback = nback[0]==1;
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
            return "Please press " + answer_key_names[0] + " or " + answer_key_names[1] + "."
          }
          if(detect&response!=answer){ //why do I have to do this?
            return "Press <strong>" + answer_key_names[1] + "</strong> when you see a T."
          }
          if(detect&response==answer){
            return "Correct!"
          }
          if(nback){ //n-back match but missed
            return "Press <strong>" + answer_key_names[0] + "</strong> when you see a match."
          }
          if(!nback&!detect&response!=null){ //not an n-back trial, but they pressed something
            return "Not a match! Press nothing."
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

  var combine = {
    timeline: test
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  timeline = [presentation_screen, instructs, combine];//, debrief];
  return timeline;

};//end of set up function