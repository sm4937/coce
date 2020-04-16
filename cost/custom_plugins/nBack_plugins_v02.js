// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function setup_nBack(loopi,n){ //called n-back but it's both that and the combine

  var stimulus_list = []; var stimnum_list = []; var feedback_list = [];
  var data_list = []; var key_list = []; var detect_list = [];
  // initialize empty lists for variable storage

  taskID = n;
  
  var stamp = '<div style="position: absolute; left: 5px; top: 5px">';
  stamp += fractals[taskID]+' width="'+stamp_size+'" height ="'+stamp_size+'"';
  stamp += '</div>';

  var prob_oddball = 0; //0.20; //kool botvinick 2014 overlapping detection task
  var prob_n_back = 0.20;
  var cdf = [prob_n_back, prob_n_back+prob_oddball,1];
  var task_label = 'combine';
  if(prob_oddball==0){
    var task_label = 'n-back';
  }

  var response_color = 'red'; //black to red upon response in n-back, detection

  // set up initial global variables
  var answer_key_names = [universal_key,'H'];
  var answer_keys = [jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[0]),jsPsych.pluginAPI.convertKeyCharacterToKeyCode('z'),jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[1])];

  //var n = 3; //set how many back they're looking for
  // it's a 3-back in Kool Botvinick 2014
  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        //var n_back_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + String(n) + ' letters ago.</p> <p>If the image on screen was not the one you saw ' + String(n) + ' letters ago, do not press anything.</p><p>If you see a T, press <strong>'+ answer_key_names[1] + '</strong> as fast as you can.</p><p>Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this task.</p> <p>Press space to begin.</p>';
        var n_back_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[0] + '</strong> every time you see the same letter you saw ' + String(n) + ' letters ago.</p> <p>If the image on screen was not the one you saw ' + String(n) + ' letters ago, do not press anything.</p><p>There might be multiple matches in a row. You must respond to every match you see.</p><p>Press any key to begin.</p>'; 
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return n_back_instructs;
      },
      trial_duration: instructs_timing //30 seconds to respond
    }]
  };
  
  //length of each block
  ntrials = 15;

  //nmatches = 3;
  rand = Math.random();
  if(rand<1){ //some blocks with more matches, randomly, approx 1/3 divided
    nmatches = 5;
  }
  if(rand<0.66){
    nmatches = 4;
  }
  if(rand<0.33){
    nmatches = 3;
  }
  var trial_type = new Array(ntrials-n);
  trial_type.fill(0);
  trial_type.fill(1,0,nmatches);
  //trial_type = [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; //add 2's in to make this the combine
  var temp = new Array(n);
  temp.fill(0);
  var trial_list = temp.concat(randomizeList(trial_type));

  for(var trial = 0; trial < ntrials; trial++){
    stimnum = Math.floor(Math.random()*num_stim);
    if(trial_list[trial]==2){
      stimnum = num_stim;
      answer_key = answer_keys[2]; //detection task
    }
    //make it flexible
    stimold = stimnum_list[trial-n]; //keep track of n trials back, starting with trial 0
    if(trial_list[trial]==1){ //force it to be an n-back trial
      stimnum = stimold;
      answer_key = answer_keys[0];
    }
    if(trial_list[trial]==0){
      answer_key = answer_keys[1]; //placeholder key
    }
    // n-backs may happen organically, so make sure this is accounted for
    if(stimnum==stimold){
      if(trial_list[trial]==1){
        answer_key = answer_keys[0];
      }
      if(trial_list[trial]==0){
        answer_key = answer_keys[1]; //placeholder key
        var stimnums = range(0,num_stim-1,1)
        var idx = stimnums.indexOf(stimold); //remove this stimulus from the list of possible stimuli
        stimnums.splice(idx,1);
        rand = Math.floor(Math.random()*(stimnums.length));
        stimnum = stimnums[rand];
      } // change the stimnum so it's not a match, thereby standardizing the number of matches each round
    }
    stimnum_list.push(stimnum);
    detect_list.push(stimnum == (num_stim));
    data_list.push((stimnum==stimold)&(stimnum!=num_stim));
    key_list.push(answer_key);
  }; //end of task building

  var test_stimuli = [];
  for(var trial = 0; trial < ntrials; trial++){
    test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], response_color: response_color, correct_key: key_list[trial], fb: feedback_list[trial], data: {nback: data_list[trial], stimnum: stimnum_list[trial], correct_key: key_list[trial], detect: detect_list[trial], n: n, nmatches: nmatches, task: task_label, tasknum: loopi}})
  };

  var combine = {
      timeline: create_color_change_timeline(test_stimuli,"",stamp)
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  var debrief = {
      type: "html-keyboard-response",
      stimulus: function() {
       overall = accuracyNback(); cutoff_message = pointsFeedback(); performance_list.push(overall);
       return "<p style='font-size:25px'>" + cutoff_message + "</p> </p><p style='font-size:25px'>Press any key to continue. </p>";
      },
      trial_duration: instructs_timing, //30 seconds to respond
      data: function(){
        dict = {task: "debrief", perf: accuracyNback()}
        return dict;
      }
  };

  timeline = [instructs,combine,debrief];
  return timeline;

}; // end of set up function 

function accuracyNback(data){ //calculate accuracy for Nback
  var correct_num = 0;
  var multiplier = 3; var nback_num = 0;
  lasttrialdata = jsPsych.data.getLastTimelineData().filter({task: 'n-back'});
  var buttons = lasttrialdata.select('key_press').values; 
 //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
  var trialtype = lasttrialdata.select('nback').values;
  ntrials = buttons.length;
  for(var trial = 0; trial < ntrials; trial++){
    if(trialtype[trial]==1){//yes it's an n-back trial
      nback_num += 1;
      if(buttons[trial] == answer_keys[0]){
       correct_num+=(1*multiplier);
     }
    }
    if(trialtype[trial]==0){ //not an n-back trial
      if(buttons[trial]==null){
       correct_num+=1}
      }
    }
  overall = Math.round((correct_num/(ntrials+(multiplier-1)*nback_num))*100);
  overall = Number.parseFloat(overall).toFixed(2);
  return overall;
}

function accuracyfinalNback(data){//calculate practice n-back accuracy to determine whether the experiment continues
    var correct_num = 0;
    var multiplier = 3; var nback_num = 0;
    lasttrialdata = jsPsych.data.get().filter({practice: true}).filter({task: 'n-back'});
    //lasttrialdata = lasttrialdata.slice(-10);
    var buttons = lasttrialdata.select('key_press').values.slice(-10); 
    //when you add other tasks on top of the n-back, you'll need to make this specific to n-back trials by filter for task:nback
    var trialtype = lasttrialdata.select('nback').values.slice(-10);
    ntrials = buttons.length;
    for(var trial = 0; trial < ntrials; trial++){
      if(trialtype[trial]==1){//yes it's an n-back trial
        nback_num += 1;
        if(buttons[trial] == answer_keys[0]){
         correct_num+=(1*multiplier);
       }
      } 
      if(trialtype[trial]==0){ //not an n-back trial
          if(buttons[trial]==null){
           correct_num+=1}
        }
      }
   overall = Math.round((correct_num/(ntrials+(multiplier-1)*nback_num))*100);
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