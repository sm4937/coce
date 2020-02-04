// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function create_experiment_instructions(){

  space_bar_message = "<p>[Press the space bar to continue.]</p>";

  instructs = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        experiment_instructs = '<p>In this experiment, you\'ll complete up to ' + String(n_repetitions) + ' short rounds of two tasks.</p>'
        experiment_instructs +='<p>If you are ' + String(cutoff_percent) + '% correct or better during one round of a task, then you pass that round and earn bonus points for it. If you are less than ' + String(cutoff_percent) + '% correct, you will not earn any points for that round.</p>'
        experiment_instructs +='<p>Keep in mind that you should be about ' + String(cutoff_percent) + '% accurate for up to ' + String(n_repetitions) + ' rounds.</p>'        
        experiment_instructs +='<p>You can earn 1 to 5 points each round.</p>';
        experiment_instructs +='<p>At the end of the experiment, the points you earned will be added up and turned into a bonus payment. $$ The more points you earn, the more money you will be paid. $$</p>'
        experiment_instructs +='<p>We can tell if you\'re not paying attention; if you\'re not, then this experiment will end early and you will earn a smaller payment. $$</p>'
        experiment_instructs +='<p>[Press space to continue.]</p>';
        //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
        return experiment_instructs;
      },
      data: {task: 'instructions'}
    }]
  };

  var full_combine_practice = {
    timeline: create_practice_tasks(exp_version)
  }

  BDM_quiz = {
    timeline: bdmQuiz(),
    conditional_function: function(){
      if(exp_version==1){
        accuracy = accuracyFinalCombine()
      }
      if(exp_version==2){
        accuracy = accuracyfinalSwitch()
      }
      return compareToCutoff(accuracy);
    }
  };


  debrief1 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        experiment_instructs = '<p>The real experiment will now start. </p>'
        experiment_instructs += '<p>For each round of the experiment, you will be asked how many points you would like to earn to complete a task.</p>'
        experiment_instructs += '<p>If you ask for more points than a randomly generated number, you will not complete that task. </p>'
        experiment_instructs += '<p>Instead, you will complete a different task and earn 1 point.</p>'
        experiment_instructs += '<p>At the end of the experiment, the points you earn will be counted and turned into a bonus payment. $$</p>'
        experiment_instructs += space_bar_message;
           //return "<p style='font-size:25px'>" + n_back_instruct
        return experiment_instructs;
        },
      data: {task: 'instructions'},
      trial_duration: instructs_timing //30 seconds to respond
    }],
    conditional_function: function(){
      return compareToCutoff(gradeBDMquiz())
    }
  };

  debrief2 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function() {
        experiment_instructs = '<p>Pay attention! If you are not at least ' + String(cutoff_percent) + '% accurate, you will not earn points for completing each task.</p>'
        experiment_instructs += '<p>You will not be getting the same amount of feedback on your performance in the real experiment.</p>'
        //experiment_instructs += '<p>If you press a button, the stimulus on screen will turn ' + response_color + ' to show you that the computer registered your response.</p>'
        experiment_instructs += '<p>If you press a button, the stimulus on screen will change color to show you that the computer registered your response.</p>'
        experiment_instructs += '<p>At the end of the experiment, we\'ll ask you to complete some optional survey questions.</p>'
        experiment_instructs += '<p>If we can tell that you are not paying attention (not accurate enough or not responding quickly enough), then this experiment will end early and you will earn a smaller payment. $$</p>'
        experiment_instructs += '<p>[Press the space bar to begin the real experiment.]</p>'
        return experiment_instructs;
        },
      data: {task: 'instructions'},
      trial_duration: instructs_timing //30 seconds to respond
    }],
    conditional_function: function(){
      return compareToCutoff(gradeBDMquiz())
    }
  };

if(debug){
  return timeline = [debrief1, debrief2]
}

//return timeline = [instructs, practice_detection, practice_1_back, practice_2_back, practice_3_back, check_accuracy, repeat_3_back, practice_combine, check_final_accuracy, repeat_combine, reach_asymptote, repeat_combine, BDM_quiz, debrief]//setup_practice_nBack(1),setup_practice_nBack(2),setup_practice_nBack(3), BDM_quiz, debrief];
//return timeline = [full_combine_practice, BDM_quiz, debrief]//setup_practice_nBack(1),setup_practice_nBack(2),setup_practice_nBack(3), BDM_quiz, debrief];
//return timeline = [instructs, practice_detection, practice_1_back, practice_2_back, practice_3_back, check_accuracy, repeat_3_back, full_combine_practice, BDM_quiz, debrief1, debrief2]//setup_practice_nBack(1),setup_practice_nBack(2),setup_practice_nBack(3), BDM_quiz, debrief];
return timeline = [instructs, full_combine_practice, BDM_quiz, debrief1, debrief2]//setup_practice_nBack(1),setup_practice_nBack(2),setup_practice_nBack(3), BDM_quiz, debrief];

};



//create dynamic practice for different versions of the experiment 

function create_practice_tasks(exp_version){
  if(exp_version==1){
    practice_flag = true;

    practice_detection = {
        timeline: setup_detection(practice_flag)
    };

    var practice_1_back = {
      timeline: setup_practice_nBack(1)
    };

    var practice_2_back = {
      timeline: setup_practice_nBack(2)
    };

    var practice_3_back = {
      timeline: setup_practice_nBack(3)
    };

    var check_accuracy = {
      timeline: [{
        type: "html-keyboard-response",
        data: {practice_accuracy: accuracyfinalNback(),task: 'instructions'},
        stimulus: function(){
          accuracy = accuracyfinalNback();
          if(accuracy>=cutoff){
            return "<p>Good job! You have reached " + String(cutoff_percent) + "% accuracy on this practice task. We will now move on.</p><p>[Press the space bar to continue.]</p>"
          } else {
            return "<p>Let's try this practice task once more. You have to be " + String(cutoff_percent) + "% accurate to progress to the real experiment.</p><p>[Press the space bar to continue.]</p>"
          }
        }
      }]
    }

    var repeat_3_back = {
      timeline: setup_practice_nBack(3),
      conditional_function: function(){
        accuracy = accuracyfinalNback();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    };

    var practice_combine = {
      timeline: setup_practice_combo()
    }  

    var reach_asymptote = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: "<p>You're not at " + String(cutoff_percent) + "% accuracy on this task.</p><p>Let's try one more time. Try your best! If we can't reach " + String(cutoff_percent) + "% accuracy, then this experment will end.</p><p>[Press the space bar to continue.]</p>",
        data: {practice_accuracy: accuracyFinalCombine(), task: 'instructions'}
      }],
      conditional_function: function(){
        accuracy = accuracyFinalCombine();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    };

    // flexible framework for repeating practice combine as many times as a subject might need it
    max_repeats = 5;
    combine_practice_timeline = [practice_combine];
    for(var reps=0;reps<max_repeats;reps++){
    starter_message = '<p>You are not at ' + String(cutoff_percent) + '% accuracy yet.</p><p>Let\'s try this task again.</p>';
    var notice = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function(){
          accuracy = accuracyFinalCombine();
          message = '<p> You were ' + String(accuracy) + '% accurate.</p>';
          message += starter_message;
          message += space_bar_message;
          return message;
        },
        data: function(){
          accuracy = accuracyFinalCombine();
          return {practice_accuracy: accuracy}
        }
      }],
      conditional_function: function(){
        accuracy = accuracyFinalCombine();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }

    var log = {
      timeline: [{
        type: "html-keyboard-response",
        data: {number_practice_combines: reps+1},
        trial_duration: instructs_timing,
        stimulus: '<p>You will get ' + String(max_repeats-reps) + ' more tries to get to ' + String(cutoff_percent) + '% accuracy on this task.</p><p>Don\'t worry, most subjects need a few rounds of practice on this task.</p><p>If the feedback is more distracting than helpful, try your best to just ignore it.</p><p>[Press space bar to try again.]</p>',
        response_ends_trial: true
      }],
      conditional_function: function(){
        accuracy = accuracyFinalCombine();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }

    var repeat_combine = {
      timeline: setup_practice_combo(),
      conditional_function: function(){
        accuracy = accuracyFinalCombine();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }
    
    combine_practice_timeline.push(notice);
    combine_practice_timeline.push(log);
    combine_practice_timeline.push(repeat_combine);
      //combine_practice_timeline.push(reach_asymptote)
    }

    var wrap_it_up = {
      timeline: [{
        type: "html-keyboard-response",
        data: {task: 'debrief'},
        stimulus: function(){
          if(accuracy>=cutoff){
            return "<p>Good job! You have reached " + String(cutoff_percent) + "% accuracy on this practice task. We will now explain how the payment process for this experiment will work.</p><p>[Press the space bar to continue.]</p>"
          }else{
            return "<p>Sorry, but you haven't reached " + String(cutoff_percent) + "% accuracy on this practice task.</p><p>You will not be moving forward into the real experiment.</p><p>[Press the space bar to continue.]</p>"
          }
        }   
      }]
    }

    combine_practice_timeline.push(wrap_it_up);

    var full_combine_practice = {
      timeline: combine_practice_timeline
    }

    return temp = [practice_detection, practice_1_back, practice_2_back, practice_3_back, check_accuracy, repeat_3_back, full_combine_practice];
  }
  // n-switch version of the instructions
  if(exp_version==2){
    var instructs = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function() {
          var initial_instructs = '<p>In this task, you will need to follow different rules depending on the color of the numbers on screen.</p>';
          initial_instructs += '<p style="color: '+ n_switch_colors[0] + '">If the number is ' + n_switch_colors[0] + ', you will need to respond to its ' + rule_names[0] + '.</p>';
          initial_instructs += '<p style="color: '+ n_switch_colors[1] + '">If the number is ' + n_switch_colors[1] + ', you will need to respond to its ' + rule_names[1] + '.</p>';
          initial_instructs += '<p>Let\'s practice responding to these rules, one at a time.</p>'
          initial_instructs += '<p>[Press the space bar to continue.]</p>'
          //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
          return initial_instructs;
        },
        data: {task: 'instructions'}
      }]
    };
    var practice_timeline = [];
    practice_timeline.push(instructs);
    var rules = [1,2];
    if(Math.random()<0.5){
      rules = [rules[1], rules[0]] // counterbalance which rule starts out
    }
    var magnitude = {
      timeline: practice1rule(rules[0])}
    var parity = {
      timeline:practice1rule(rules[1])} 
    practice_timeline.push(magnitude);
    practice_timeline.push(parity);
    var together = {
      timeline: setup_nswitch(true,-1,2) //easy block
    }
    var practice_hard = {
      timeline: setup_nswitch(true,-1,1) //hard block
    }
    //temp.push(together);

    max_repeats = 3;
    practice_timeline.push(together)
    practice_timeline.push(practice_hard)
    for(var reps=0;reps<max_repeats;reps++){
    starter_message = '<p>You are not at ' + String(cutoff_percent) + '% accuracy yet.</p><p>Let\'s try this task again.</p>';
    var notice = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function(){
          accuracy = accuracySwitch();
          message = '<p> You were ' + String(accuracy) + '% accurate.</p>';
          message += starter_message;
          message += space_bar_message;
          return message;
        },
        data: function(){
          accuracy = accuracySwitch();
          return {practice_accuracy: accuracy}
        }
      }],
      conditional_function: function(){
        accuracy = accuracyfinalSwitch();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }

    var log = {
      timeline: [{
        type: "html-keyboard-response",
        data: {number_practice_combines: reps+1},
        trial_duration: instructs_timing,
        stimulus: '<p>You will get ' + String(max_repeats-reps) + ' more tries to get to ' + String(cutoff_percent) + '% accuracy on this task.</p><p>Don\'t worry, most subjects need a few rounds of practice on this task.</p><p>If the feedback is more distracting than helpful, try your best to just ignore it.</p><p>[Press space bar to try again.]</p>',
        response_ends_trial: true
      }],
      conditional_function: function(){
        accuracy = accuracyfinalSwitch();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }

    var repeat_task = {
      timeline: setup_nswitch(true,-1,1),
      conditional_function: function(){
        accuracy = accuracyfinalSwitch();
        if(accuracy<cutoff){
          return true;
        }else{
          return false;
        }
      }
    }
    
    practice_timeline.push(notice);
    practice_timeline.push(log);
    practice_timeline.push(repeat_task);
      //combine_practice_timeline.push(reach_asymptote)
    }

    var wrap_it_up = {
      timeline: [{
        type: "html-keyboard-response",
        data: {task: 'debrief'},
        stimulus: function(){
          if(accuracy>=cutoff){
            return "<p>Good job! You have reached " + String(cutoff_percent) + "% accuracy on this practice task. We will now explain how the payment process for this experiment will work.</p><p>[Press the space bar to continue.]</p>"
          }else{
            return "<p>Sorry, but you haven't reached " + String(cutoff_percent) + "% accuracy on this practice task.</p><p>You will not be moving forward into the real experiment.</p><p>[Press the space bar to continue.]</p>"
          }
        }   
      }]
    }

    practice_timeline.push(wrap_it_up);

    return practice_timeline;
  }
}

// code graveyard - counterbalance practice games for random presentation order
  /*counterb = 1 + Math.round(Math.random());
  practice_nback1 = {
    timeline: setup_nBack(practice_flag),
    conditional_function: function(){
      if(counterb == 2){
        return true // run n-back first
      } else { //ie condition = 1
        return false //skip n-back first
      }
    }
  };
  console.log('practice n-back created')
  practice_detection = {
      timeline: setup_detection(practice_flag)
    };

  practice_nback2 = {
    timeline: setup_nBack(practice_flag),
    conditional_function: function(){
      if(counterb == 2){
        return false //skip second n-back
      } else {
        return true // run n-back second
      }
    }
  };*/

  /*var check_final_accuracy = {
    timeline: [{
      type: "html-keyboard-response",
      data: {practice_accuracy: accuracyFinalCombine(), task: 'instructions'},
      stimulus: function(){
        return "<p>Good job! You have reached " + String(cutoff_percent) + "% accuracy on this practice task. We will now move on.</p><p>[Press the space bar to continue.]</p>"
      }    
    }],
    conditional_function: function(){
      accuracy = accuracyFinalCombine();
      if(accuracy>=cutoff){
        return true;
      }
      if(accuracy<cutoff){
        return false;
      }
    }
  }*/