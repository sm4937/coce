// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function setup_nBack(loopi){ //called n-back but it's the combine
  var stimulus_list = []; var stimnum_list = []; var feedback_list = [];
  var data_list = []; var key_list = []; var detect_list = [];
  // initialize empty lists for variable storage

  var prob_oddball = 0.20; //kool botvinick 2014
  var prob_n_back = 0.20;
  var cdf = [prob_n_back, prob_n_back+prob_oddball,1];

  var response_color = 'red'; //black to red upon response in n-back, detection

  var n = 3; //set how many back they're looking for
  // it's a 3-back in Kool Botvinick 2014
  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        var n_back_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + String(n) + ' letters ago.</p> <p>If the image on screen was not the one you saw ' + String(n) + ' letters ago, do not press anything.</p><p>If you see a T, press <strong>'+ answer_key_names[1] + '</strong> as fast as you can.</p><p>Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this task.</p> <p>Press space to begin.</p>';
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return n_back_instructs;
      },
      trial_duration: instructs_timing //30 seconds to respond
    }]
  };
  
  //length of each block
  ntrials = 15;

  trial_type = [1, 1, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0];
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

  /*for(var trial = 0; trial < ntrials; trial++){
    rand = Math.random();
    trial_type = 0; //standard is non-event
    if(rand<cdf[1]){
      trial_type = 2; //oddball
    }
    if((rand<cdf[0])&(trial>n)){
      trial_type = 1; //n-back
    }
    trial_list.push(trial_type);
  }*/

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
    if(trial_list[trial]==0){
      answer_key = answer_keys[1]; //placeholder key
    }
    // n-backs may happen organically, so make sure this is accounted for
    if(stimnum==stimold){
      answer_key = answer_keys[0];
    }
    stimnum_list.push(stimnum);
    detect_list.push(stimnum == (num_stim));
    data_list.push((stimnum==stimold)&(stimnum!=num_stim));
    key_list.push(answer_key);
  }; //end of task building

  var test_stimuli = [];
  for(var trial = 0; trial < ntrials; trial++){
    test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], response_color: response_color, correct_key: key_list[trial], fb: feedback_list[trial], data: {nback: data_list[trial], correct_key: key_list[trial], detect: detect_list[trial], n: 3, task: 'combine', tasknum: loopi}})
  };

  var combine = {
      timeline: create_color_change_timeline(test_stimuli,"")
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  var debrief = {
      type: "html-keyboard-response",
      stimulus: function() {
       overall = accuracyCombine(); cutoff_message = pointsFeedback(); performance_list.push(overall);
       return "<p style='font-size:25px'>" + cutoff_message + "</p> </p><p style='font-size:25px'>Press the space bar to continue. </p>";
      },
      trial_duration: instructs_timing //30 seconds to respond
  };

  var presentation_screen = {
    type: "html-keyboard-response",
    stimulus: '+',
    trial_duration: 100
  };

  timeline = [presentation_screen,instructs,combine,debrief];
  return timeline;

}; // end of set up function 

function accuracyNback(data){ //calculate accuracy for Nback
  var correct_num = 0;
  lasttrialdata = jsPsych.data.getLastTimelineData().filter({task: 'n-back'});
  var buttons = lasttrialdata.select('key_press').values; 
 //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
  var trialtype = lasttrialdata.select('nback').values;
  ntrials = trialtype.length;
  for(var trial = 0; trial < ntrials; trial++){
    if(trialtype[trial]){//yes it's an n-back trial
      if(buttons[trial] == answer_keys[0]){
       correct_num+=1
     }
    }
    if(!trialtype[trial]){ //not an n-back trial
      if(buttons[trial]==null){
       correct_num+=1}
      }
    }
  overall = Math.round((correct_num/ntrials)*100);
  overall = Number.parseFloat(overall).toFixed(2);
  return overall;
}

function accuracyfinalNback(data){//calculate practice n-back accuracy to determine whether the experiment continues
    var correct_num = 0;
    lasttrialdata = jsPsych.data.get().filter({n: 3}).filter({practice: true}).filter({task: 'n-back'});
    //lasttrialdata = lasttrialdata.slice(-10);
    var buttons = lasttrialdata.select('key_press').values.slice(-10); 
    //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
    var trialtype = lasttrialdata.select('nback').values.slice(-10);
    ntrials = 10;
    for(var trial = 0; trial < ntrials; trial++){
      if(trialtype[trial]){//yes it's an n-back trial
        if(buttons[trial] == answer_keys[0]){
         correct_num+=1
       }
      } 
      if(!trialtype[trial]){ //not an n-back trial
          if(buttons[trial]==null){
           correct_num+=1}
        }
      }
   overall = Math.round((correct_num/ntrials)*100);
   overall = Number.parseFloat(overall).toFixed(2);
   return overall;    
}


  function accuracyFinalCombine(){//calculate practice n-back accuracy to determine whether the experiment continues
    var correct_num = 0;
    lasttrialdata = jsPsych.data.get().filter({n: 3}).filter({practice: true}).filter({task: 'combine'});
    //lasttrialdata = lasttrialdata.slice(-10);
    ntrials = 15;
    var buttons = lasttrialdata.select('key_press').values.slice(-ntrials); 
    //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
    var trialtype = lasttrialdata.select('nback').values.slice(-ntrials);
    var detect = lasttrialdata.select('detect').values.slice(-ntrials);
    for(var trial = 0; trial < ntrials; trial++){
      if(trialtype[trial]==1){//yes it's an n-back trial
        if(buttons[trial] == answer_keys[0]){
         correct_num+=1
       }
      }
      if(detect[trial]){
        if(buttons[trial]==answer_keys[2]){
          correct_num+=1;
        }
      } 
      if(trialtype[trial]==0&!detect[trial]){ //not an n-back or detect trial
        if(buttons[trial]==null){
         correct_num+=1}
      }
    }
    overall = Math.round((correct_num/ntrials)*100);
    overall = Number.parseFloat(overall).toFixed(2);
    return overall;    
  }

function accuracyCombine(data){ //calculate accuracy for Nback
  var correct_num = 0;
  lasttrialdata = jsPsych.data.getLastTimelineData();
  var buttons = lasttrialdata.filter({task: 'combine'}).select('key_press').values; 
 //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
  var trialtype = lasttrialdata.filter({task: 'combine'}).select('nback').values;
  var detect = lasttrialdata.filter({task: 'combine'}).select('detect').values;
  ntrials = trialtype.length;
  for(var trial = 0; trial < ntrials; trial++){
    if(trialtype[trial]==1){//yes it's an n-back trial
      if(buttons[trial] == answer_keys[0]){
       correct_num+=1
     }
    }
    if(detect[trial]){
      if(buttons[trial]==answer_keys[2]){
        correct_num+=1;
      }
    } 
    if(trialtype[trial]==0&!detect[trial]){ //not an n-back or detect trial
      if(buttons[trial]==null){
       correct_num+=1}
    }
  }
  overall = Math.round((correct_num/ntrials)*100);
  overall = Number.parseFloat(overall).toFixed(2);
  return overall;
}

//code graveyard

// this first thing is the code I wrote to have the stimulus change color on response, 
// which only didn't work because the feedback was smaller than the stimuli themselves
// could still work in theory
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