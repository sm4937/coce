// functions and such just for detection creation 
// set up the task parameters and stimulus list
function setup_test_detection(loopi){
  var p_target = 0.2; //set probability of target stimulus
  detection_target = "T";

  var stimulus_list = []; var stimnum_list = []; 
  var data_list = []; var key_list = []; //initialize empty lists for task parameters
  
  var presentation_screen = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: fractals[1],
      prompt: '<p>Practice task.</p><p>This picture will always be associated with this task, like a picture label.</p><p>Press the space bar to continue.</p>'
    }]
  };

    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
          stimulus: function() {
            var initial_instructs = '<p>In this task, you will need to press <strong>' + answer_key_names[1] + '</strong> every time you see the target letter <strong>' + detection_target + '</strong>. </p><p>You must respond while the letter is on screen or your answer will not count.</p><p> Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing this round.</p> <p>Press space to begin.</p>';
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
         return "<p>Your accuracy was <strong>" + overall + "%</strong>.</p><p>You will need to be at least " + String(cutoff_percent) + "% accurate to earn points in the real experiment.</p><p>In the real experiment, you will not receive the same feedback as during practice.</p><p>Press the space bar to continue. </p>";
        },
        data: function(){
          accuracy = accuracyDetection();
          return {training_performance: accuracy, task: 'debrief'}
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
        test_stimuli.push({stimulus: stimuli[stimnum_list[trial]], correct_key: key_list[trial], data: {detect: data_list[trial], correct_key: key_list[trial], task: 'detection', tasknum: -1}})
    };

    test = create_color_change_timeline(test_stimuli);

    var detect = {
      timeline: test
    } // do not TOUCH THIS PIECE OF CODE HOLY mother of god!

  timeline = [presentation_screen,instructs,detect,debrief];
  return timeline

};