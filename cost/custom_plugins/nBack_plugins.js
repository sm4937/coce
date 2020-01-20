// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function setup_nBack(practice){
  var stimulus_list = []; var stimnum_list = []; 
  var data_list = []; var key_list = []; var feedback_list = []; var detect_list = [];
  // initialize empty lists for variable storage

  var prob_oddball = 0.20; //kool botvinick 2014

  //also, in kool botvinick 2014 it's a 20% chance of an n-back match trial
  //soft code this 

  var n = 3; //set how many back they're looking for
  // it's a 3-back in Kool Botvinick 2014
  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        var n_back_instructs = '<p>In this task, you will need to press ' + answer_key_names[0] + ' every time you see the same letter you saw ' + String(n) + ' trial(s) ago.</p> <p>If the image on screen was not the one you saw ' + String(n) + ' trials ago, do not press anything.</p><p>If you see a T, ignore the other rules and just press '+ answer_key_names[1] + ' as fast as you can.</p><p>Pay attention! If you are not ' + String(cutoff) + '% accurate, you will not earn points for completing this task.</p> <p>Press space to begin.</p>';
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return n_back_instructs;
          }
    }]
  };

  if(!practice){ //so regular task, not practice parameters
    ntrials = 15; 
    for (var trial = 0; trial < ntrials; trial++) {
      stimnum = Math.floor(Math.random()*num_stim);
      if(Math.random() < prob_oddball){
        stimnum = num_stim+1;
      }
      stimnum_list.push(stimnum);
      //stimulus_list.push(makeHTML(stimuli[stimnum]));
      //feedback_list.push(fb_img_paths[stimnum]);
      detect_list.push(stimnum == num_stim+1); //last stimulus
      if(trial < n){
      key_list.push(null);
      data_list.push(false)}
      //make it flexible
      if(trial >= n){
        stimold = stimnum_list[trial-n]; //keep track of n trials back
        data_list.push(stimnum == stimold);
        answer_key = null; //withheld response
        if(stimnum==stimold){
          answer_key = answer_keys[0];
        } 
        if(stimnum == num_stim+1){
           answer_key = answer_keys[2];
         }
        key_list.push(answer_key);
      }
    };

    var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
       test_stimuli.push({stimulus: makeHTML(stimuli[stimnum_list[trial]]), correct_key: key_list[trial], fb: feedback_list[trial], data: {nback: data_list[trial], correct_key: key_list[trial], detect: detect_list[trial], task: 'n-back', tasknum: loopi}})
    };      

  /*var test = {
    timeline: [{
      type: "categorize-html",
        //choices: [answer_keys],
        stimulus: jsPsych.timelineVariable('stimulus'),
        data: jsPsych.timelineVariable('data'),
        key_answer: jsPsych.timelineVariable('correct_key'),
        stimulus_duration: 1500,
        trial_duration: 1500,
        show_feedback_on_timeout: false,
        feedback_duration: 250,
        timeout_message: ' ',
        correct_text: jsPsych.timelineVariable('fb'),
        incorrect_text: jsPsych.timelineVariable('fb'),
        show_stim_with_feedback: false
      }],
      timeline_variables: test_stimuli
    };*/ 
    //categorize HTML can't keep the stimulus on screen when the trial is over, so let's see if we can make it work with 
    //html keyboard response since we code for accuracy ourselves anyway

    var test = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: jsPsych.timelineVariable('stimulus'),
        data: jsPsych.timelineVariable('data'),
        trial_duration: 1750,
        stimulus_duration: 1500,
        show_feedback_on_timeout: false,
        response_ends_trial: false
      }],
      timeline_variables: test_stimuli
    };

    var debrief = {
       type: "html-keyboard-response",
       stimulus: function() {
         overall = accuracyNback(); cutoff_message = pointsFeedback();
         return "<p style='font-size:25px'>Your accuracy was <strong>" + overall + "%</strong>.</p><p style='font-size:25px'>" + cutoff_message + "</p> </p><p style='font-size:25px'>Press the space bar to continue. </p>";
        }
      };

    var presentation_screen = {
      type: "html-keyboard-response",
      stimulus: '+',
      trial_duration: 100
      };

  }else{ // put in practice paradigm
    ntrials = 10; 

    //var data_list = [Array(n).fill(false)];//first trial always false w/ n-back
    //var key_list = [Array(n).fill(answer_keys[1])];
    for (var trial = 0; trial < ntrials; trial++) {
      stimnum = Math.round(Math.random());
      stimnum_list.push(stimnum);
      stimulus_list.push(makeHTML(stimuli[stimnum]));
      if(trial < n){
      key_list.push(answer_keys[1]);
      data_list.push(false)}
      //make it flexible
      if(trial >= n){
        stimold = stimnum_list[trial-n]; //keep track of n trials back, starting with trial 0
        data_list.push(stimnum==stimold);
        detect_list.push(stimnum==2);
          if(stimnum==stimold){
            answer_key = answer_keys[0];
              } else {
            answer_key = answer_keys[1];
               }
            key_list.push(answer_key);
          }
      };

    var presentation_screen = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: fractals[0],
        prompt: '<p>Practice game.</p><p>Press the space bar to continue.</p>'
          }]
      };

    //word_instructs = instructs;
    //instructs = [presentation_screen, instructs]; 
    var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
       test_stimuli.push({stimulus: stimulus_list[trial], correct_key: key_list[trial], data: {nback: data_list[trial], detect: detect_list[trial], correct_key: key_list[trial], task: 'n-back', tasknum: -1}})
    };      
    var debrief = {
       type: "html-keyboard-response",
       stimulus: function() {
         overall = accuracyNback(); performance_list.push(overall);
         return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be 80% accurate during the study to earn points.</p><p>Please remember that you will not receive the same amount of feedback during the real experiment.</p><p>Press the space bar to continue. </p>";
        }
    };

    var test = {
    timeline: [{
      type: "categorize-html",
      choices: answer_keys,
      stimulus: jsPsych.timelineVariable('stimulus'),
      data: jsPsych.timelineVariable('data'),
      key_answer: jsPsych.timelineVariable('correct_key'),
      stimulus_duration: 1500,
      trial_duration: 1500,
      show_feedback_on_timeout: false,
      timeout_message: ' ',
      feedback_duration: 750,
      correct_text: 'Correct',
      show_stim_with_feedback: false,
      incorrect_text: function(){
        response = jsPsych.data.getLastTrialData().select("key_press").values;
        answer = jsPsych.data.getLastTrialData().select("correct_key").values;
        type = jsPsych.data.getLastTrialData().select("nback").values;
        detect = jsPsych.data.getLastTrialData().select("detect").values;
        if(type[0]){
            return "Wrong response key!"
          }
        if(!type[0]){ //not an n-back trial
            return "Not a match!"
          } 
        if(detect[0]) {
            return "Wrong response key!"
          }
        },
      }],
      timeline_variables: test_stimuli
    //sample: {type: 'fixed-repetitions', size: 2},
    //randomization: true
    };

  };//end of practice set up

  timeline = [presentation_screen,instructs,test,debrief];
  return timeline;

}; // end of set up function 

  function accuracyNback(data){ //calculate accuracy for Nback
    var correct_num = 0;
    lasttrialdata = jsPsych.data.getLastTimelineData();
    var buttons = lasttrialdata.select('key_press').values; 
    //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
    var trialtype = lasttrialdata.select('nback').values;
    for(var trial = 0; trial < ntrials; trial++){
      if(trialtype[trial]){//yes it's an n-back trial
        if(buttons[trial] == answer_keys[0]){
         correct_num+=1
       }
      } 
      if(!trialtype[trial]) { //not an n-back trial
          if(buttons[trial]==null){
           correct_num+=1}
        }
      }
   overall = Math.round((correct_num/ntrials)*100);
   overall = Number.parseFloat(overall).toFixed(2);
   return overall;
   }
