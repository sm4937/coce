// functions and such just for detection creation 
// set up the task parameters and stimulus list

function setup_nswitch(practice,loopi,condition){
  var stimulus_list = []; var stimnum_list = []; color_list = [];
  var data_list = []; var key_list = []; rule_list = []; //initialize empty lists for task parameters

  var stimuli = ["1","2","3","4","6","7","8","9"];
  var num_stim = stimuli.length;
  var response_color = 'black';
  var colors = n_switch_colors;
 
  var rules = [1,2];
  if(Math.random()<0.5){
    rules = [rules[1], rules[0]] // counterbalance which rule starts out
  }
  
  // condition = 1 if hard, 2 if easy
  if(condition==2){
    p_switch = 0.1;
  }else{ //condition is 1
    p_switch = 0.9;
  }

  fractal = fractals[2+condition]; //fractal 3 for hard, 4 for easy

  //answer keys are 
  var answer_key_names = ['S','F','H','K']
  var response_names = n_switch_response_names;
  var answer_keys = [jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[0]),jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[1]),jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[2]),
      jsPsych.pluginAPI.convertKeyCharacterToKeyCode(answer_key_names[3])];
  labels = answer_key_names;

  //create div tags for displaying rules with correct buttons
  var key_presses = '<div>';
  var left_offsets = [-100, -50, 50, 100]
    for(var j=0; j < labels.length; j++){
      var width = 100/(labels.length);
      var left_offset = left_offsets[j]//(j * (100 /(labels.length))) - (width*2);
      key_presses += '<div style="display: inline-block; position: relative; left:'+left_offset+'%; text-align: center; width: '+width+'%;">';
      key_presses += '<span style="text-align: center; font-size: 25px;">'+labels[j]+'</span>';
      key_presses += '</div>'
    }
  key_presses += '</div>';
    //some css to make spaced out labels work, instead of building new plugin entirely
    // now layering key names with rules
  key_presses += '<div>';
  for(var j=0; j < labels.length; j++){
    var width = 100/(labels.length);
    var left_offset = left_offsets[j]//(j * (100 /(labels.length))) - (width*2);
    key_presses += '<div style="display: inline-block; position: relative; left:'+left_offset+'%; text-align: center; width: '+width+'%;">';
    key_presses += '<span style="text-align: center; font-size: 25px;">'+response_names[j]+'</span>';
    key_presses += '</div>'
  }
  key_presses += '</div>';
     
  space_bar_message = "<p> [Press the space bar to continue.] </p>";

  if(!practice){ //regular task
    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function() {
          var initial_instructs = '<p>In this task, you will need to follow different rules depending on the color of the numbers on screen.</p>';
          initial_instructs += '<p style="color: '+ colors[0] + '">If the number is ' + colors[0] + ', you will need to respond to its ' + rule_names[0] + '.</p>';
          initial_instructs += '<p style="color: '+ colors[1] + '">If the number is ' + colors[1] + ', you will need to respond to its ' + rule_names[1] + '.</p>';
          initial_instructs += '<p> As a reminder, the response keys are: </p>'
          initial_instructs += '<p><strong> S   F   H   K </strong></p>';
          initial_instructs += '<p>Please place your fingers on the response keys now.</p>'
          initial_instructs += '<p>Press any button to begin.</p>';
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
        overall = accuracySwitch(); cutoff_message = pointsFeedback(); performance_list.push(overall);
        return "<p>" + cutoff_message + "</p> </p><p>Press the space bar to continue. </p>";
      },
      data: function(){
        dict = {task: 'debrief', perf: accuracySwitch()};
        return dict;
      },
      trial_duration: instructs_timing //30 seconds to respond
    };

  ntrials = 15;
  practice_flag = false;

  } else { // make a practice task

    var presentation_screen = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: fractal,
        prompt: '<p>Practice task.</p><p>This picture will always be associated with the following task, like a picture label.</p><p>You will now learn the rules of this task and get a chance to practice.</p><p>Press the space bar to continue.</p>',
        trial_duration: instructs_timing
          }]
      };

    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
          stimulus: function() {
            var initial_instructs = '<p>Now, let\'s practice responding to both rules.</p>'
            initial_instructs += '<p>You will be responding to both rules, one at a time, depending on the color of the number on screen.</p>';
            initial_instructs += '<p>You will use four different keys, <strong>S, F, H, and K</strong> to respond.</p>'
            initial_instructs += '<p style="color: '+ colors[0] + '">If the number is ' + colors[0] + ', you will need to respond to its ' + rule_names[0] + '.</p>';
            initial_instructs += '<p style="color: '+ colors[1] + '">If the number is ' + colors[1] + ', you will need to respond to its ' + rule_names[1] + '.</p>';
            initial_instructs += '<p>Press the space bar to practice this task.</p>';
            return initial_instructs;
              },
            data: {task: 'instructions'},
            trial_duration: instructs_timing
            }]
          };    

      var debrief = {
       type: "html-keyboard-response",
       stimulus: function() {
         overall = accuracySwitch();
         return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate to earn points in the real experiment.</p><p>In the real experiment, you will not receive the same feedback as during practice.</p><p>Press the space bar to continue. </p>";
        },
        data: function(){
          accuracy = accuracySwitch();
          return {practice_accuracy: accuracy, practice: true, task: 'debrief'}
        },
        trial_duration: instructs_timing
      }

    ntrials = 15;
    practice_flag = true;

  } // end practice task creation

  trial_list = [];

  for(var trial = 0; trial < ntrials; trial++){
    rand = Math.random();
    //console.log(rand)
    sw = rand<p_switch;
    //console.log(sw)
    trial_list.push(sw);
    //console.log(trial_type)
  }

  /*if(practice){
    trial_list = [false, false, false, false, false, false, false, false, false, false];
    for(var j = 0; j <= 3; j++){ //up to 3 random switch trials in practice block to practice that
      idx = Math.floor(Math.random()*ntrials);
      trial_list[idx] = true; //add one random switch trial in there
      practice_flag = true;
    }
  }*/

  rule = rules[0];
  color = colors[rule-1];
  for(var trial = 0; trial < ntrials; trial++){
    stimnum = Math.floor(Math.random()*num_stim);
    stim = stimuli[stimnum]*1;
    //make it flexible
    if(trial_list[trial]){ //it's a switch trial
      rules = [rules[1], rules[0]];
      rule = rules[0];
      color = colors[rule-1];
    }
    if(rule==1){ //magnitude
      if(stim<=4){
        answer_key = answer_keys[n_switch_response_names.indexOf("smaller")];
      }if(stim>5){
        answer_key = answer_keys[n_switch_response_names.indexOf("larger")];
      }
    }
    if(rule==2){ //parity
      if(stim%2==0){
        answer_key = answer_keys[n_switch_response_names.indexOf("even")];
      }if(stim%2>0){
        answer_key = answer_keys[n_switch_response_names.indexOf("odd")];
      }
    }
    stimnum_list.push(stimnum);
    data_list.push(trial_list[trial]);
    key_list.push(answer_key);
    rule_list.push(rule);
    color_list.push(color);
  }; //end of task building

  var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
      test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], color: color_list[trial], response_color: response_color, correct_key: key_list[trial], data: {switch: data_list[trial], rule: rule_list[trial], stimnum: stimnum_list[trial], task: 'n-switch', prob_switch: p_switch, correct_key: key_list[trial], tasknum: loopi, practice: practice_flag}})
  }; 

  if(!practice){

  var nswitch = {
    timeline: create_color_change_timeline(test_stimuli,key_presses)   
  }

  }else{
    test = []; // empty test variable for practice, build more dynamic feedback by hijacking the plugins
    for(var i = 0; i < test_stimuli.length; i++){
      var stimulus_screen = {
        type: "categorize-html",
        stimulus: '<p style="color: ' + test_stimuli[i]['color'] + '">' + test_stimuli[i]['stimulus'] + '</p>',
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
        stimulus: function(){
          correct = jsPsych.data.getLastTrialData().select("correct").values; 
          correct = correct[0];
          response = jsPsych.data.getLastTrialData().select("key_press").values;
          response = response[0]; // remove list container
          answer = jsPsych.data.getLastTrialData().select("correct_key").values;
          answer = answer[0];
          sw = jsPsych.data.getLastTrialData().select("switch").values;
          sw = sw[0]==1;
          if(correct){
            return "Correct!"
          }else{
            if(response==answer){
              return "Correct!"
            }
            if(response == null){
              return "Please respond in time!"
            }
            if(!answer_keys.includes(response)){ // didn't press a valid key
              return "Please press S, F, H or K."
            }
            if(response!=answer){
              return "Incorrect!"
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

    var nswitch = {
      timeline: test
    } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  }; //end of practice task creation

  timeline = [presentation_screen,instructs,nswitch,debrief];
  return timeline
};

function accuracySwitch(data){ //calculate accuracy for detection task (0-back)
  lasttrialdata = jsPsych.data.getLastTimelineData().filter({task: "n-switch"});
  correct_list = lasttrialdata.select('correct').values;
  accuracy = turnLogicalIntoAccuracy(correct_list);
  return Number.parseFloat(accuracy*100).toFixed(2);
 };

function accuracyfinalSwitch(){
  lasttrialdata = jsPsych.data.get().filter({practice: true}).filter({task: 'n-switch'});
  ntrials = 15;
  correct_list = lasttrialdata.select('correct').values.slice(-ntrials);
  accuracy = turnLogicalIntoAccuracy(correct_list);
  return Number.parseFloat(accuracy*100).toFixed(2);
}

function practice1rule(rule){
  var stimulus_list = []; var stimnum_list = []; color_list = [];
  var data_list = []; var key_list = []; rule_list = []; //initialize empty lists for task parameters

  var stimuli = ["1","2","3","4","6","7","8","9"];
  var num_stim = stimuli.length;
  var response_color = 'black';
  var colors = n_switch_colors;
  var blank_placeholder = '<p style="color:white"> _ </p>';
  var answer_key_names = ["S","F","H","K"];
  var answer_key_labels = [blank_placeholder,blank_placeholder,blank_placeholder,blank_placeholder];
  var response_names = [blank_placeholder,blank_placeholder,blank_placeholder,blank_placeholder];

  //fix indices appropriately for each rule
  if(rule==1){ //magnitude
    var indices = [n_switch_response_names.indexOf("smaller"),n_switch_response_names.indexOf("larger")];
    answer_key_labels[indices[0]] = answer_key_names[indices[0]];
    answer_key_labels[indices[1]] = answer_key_names[indices[1]];
    response_names[indices[0]] = 'smaller';
    response_names[indices[1]] = 'larger';
  }
  if(rule==2){ //parity
    var indices = [n_switch_response_names.indexOf("even"),n_switch_response_names.indexOf("odd")];
    answer_key_labels[indices[0]] = answer_key_names[indices[0]];
    answer_key_labels[indices[1]] = answer_key_names[indices[1]];
    response_names[indices[0]] = 'even';
    response_names[indices[1]] = 'odd';
  }
  var answer_keys = [jsPsych.pluginAPI.convertKeyCharacterToKeyCode('S'),jsPsych.pluginAPI.convertKeyCharacterToKeyCode('F'),jsPsych.pluginAPI.convertKeyCharacterToKeyCode('H'),
      jsPsych.pluginAPI.convertKeyCharacterToKeyCode('K')];
  labels = answer_key_labels;

  //create div tags for displaying rules with correct buttons
  var key_presses = '<div>';
  var left_offsets = [-100, -50, 50, 100]
    for(var j=0; j < labels.length; j++){
      var width = 100/(labels.length);
      var left_offset = left_offsets[j]//(j * (100 /(labels.length))) - (width*2);
      key_presses += '<div style="display: inline-block; position: relative; left:'+left_offset+'%; text-align: center; width: '+width+'%;">';
      key_presses += '<span style="text-align: center; font-size: 25px;">'+labels[j]+'</span>';
      key_presses += '</div>'
    }
  key_presses += '</div>';
    //some css to make spaced out labels work, instead of building new plugin entirely
    // now layering key names with rules
  key_presses += '<div>';
  for(var j=0; j < labels.length; j++){
    var width = 100/(labels.length);
    var left_offset = left_offsets[j]//(j * (100 /(labels.length))) - (width*2);
    key_presses += '<div style="display: inline-block; position: relative; left:'+left_offset+'%; text-align: center; width: '+width+'%;">';
    key_presses += '<span style="text-align: center; font-size: 25px;">'+response_names[j]+'</span>';
    key_presses += '</div>'
  }
  key_presses += '</div>';
     
  space_bar_message = "<p> [Press the space bar to continue.] </p>";

  var instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        initial_instructs = '<p>Let\'s practice following the '  + rule_names[rule-1] + ' rule.</p>'
        initial_instructs += '<p>On the next screen, please respond to the ' + rule_names[rule-1] + ' of the numbers on screen.</p>'
        initial_instructs += '<p>Press <strong>' + answer_key_names[indices[0]] + '</strong> if the number is ' + response_names[indices[0]] + '.</p>'
        initial_instructs += '<p>Press <strong>' + answer_key_names[indices[1]] + '</strong> if the number is ' + response_names[indices[1]] + '.</p>'
        initial_instructs += '<p>You will have ' + String(stim_timing) + ' milliseconds to respond to each number. Try to respond as quickly as you can.</p>'
        initial_instructs += space_bar_message;
        return initial_instructs;
      },
      data: {task: 'instructions'},
      trial_duration: instructs_timing //30 seconds to respond
      }]
    };

  var debrief = {
    type: "html-keyboard-response",
    stimulus: function() {
      overall = accuracySwitch();
      return "<p>You were <strong>" + overall + "% accurate</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate to earn points in the real experiment.</p><p>In the real experiment, you will not receive the same feedback as during practice.</p><p>Press the space bar to continue. </p>";
    },
    data: function(){
      accuracy = accuracySwitch();
      return {practice_accuracy: accuracy, practice: true, task: 'debrief'}
    },
    trial_duration: instructs_timing
  }

  ntrials = 8;

  var data_list = new Array(ntrials); var rule_list = new Array(ntrials);
  data_list.fill(false);
  rule_list.fill(rule);

  color = colors[rule-1];
  for(var trial = 0; trial < ntrials; trial++){
    stimnum = Math.floor(Math.random()*num_stim);
    stim = stimuli[stimnum]*1;
    //make it flexible
    if(rule==1){ //magnitude
      if(stim<=4){
        answer_key = answer_keys[indices[0]];
      }if(stim>5){
        answer_key = answer_keys[indices[1]];
      }
    }
    if(rule==2){ //parity
      if(stim%2==0){
        answer_key = answer_keys[indices[0]];
      }if(stim%2>0){
        answer_key = answer_keys[indices[1]];
      }
    }
    stimnum_list.push(stimnum);
    key_list.push(answer_key);
    color_list.push(color);
  }; //end of task building

  var test_stimuli = [];
    for(var trial = 0; trial < ntrials; trial++){
      test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], color: color_list[trial], response_color: response_color, correct_key: key_list[trial], data: {switch: data_list[trial], rule: rule_list[trial], stimnum: stimnum_list[trial], task: 'n-switch', prob_switch: 0, correct_key: key_list[trial], tasknum: -1, practice: true}})
  }; 

  test = []; // empty test variable for practice, build more dynamic feedback by hijacking the plugins
  for(var i = 0; i < test_stimuli.length; i++){
    var stimulus_screen = {
      type: "categorize-html",
      stimulus: '<p style="color: ' + test_stimuli[i]['color'] + '">' + test_stimuli[i]['stimulus'] + '</p>',
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

    var feedback_screen = {
      type: "html-keyboard-response",
      stimulus: function(){
        correct = jsPsych.data.getLastTrialData().select("correct").values; 
        correct = correct[0];
        response = jsPsych.data.getLastTrialData().select("key_press").values;
        response = response[0]; // remove list container
        answer = jsPsych.data.getLastTrialData().select("correct_key").values;
        answer = answer[0];
        sw = jsPsych.data.getLastTrialData().select("switch").values;
        sw = sw[0]==1;
        if(correct){
          return "Correct!"
        }else{
          if(response==answer){
            return "Correct!"
          }
          if(response == null){
            return "Please respond in time!"
          }
          if(!answer_keys.includes(response)){ // didn't press a valid key
            return "Please press " + answer_key_names[indices[0]] + " or " + answer_key_names[indices[1]] + "."
          }
          if(response!=answer){
            return "Incorrect!"
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

  var nswitch = {
     timeline: test
  } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  timeline = [instructs,nswitch,debrief];
  return timeline
};