// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function bdmQuiz(){
  space_bar_message = "<p>[Press the space bar to continue.]</p>";

  screen1 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>In this experiment, you will be asked to repeatedly complete two of the tasks you just practiced.</p>'
          experiment_instructs +='<p>The next few screens provide detailed instructions about the way we will determine your final payment.</p>'
          experiment_instructs +='<p>The way we will determine your final payment is a bit unusual, so we are going to try and explain it well.</p>'
          experiment_instructs +='<p>It is very important that you understand the explanation, since your bonus payment $$ will depend on your ability to make good decisions.</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
            }
          }]
        };

  /*screen2 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>In this task, you\'ll complete up to ' + String(n_repetitions) + ' short rounds of two tasks.</p>'
          experiment_instructs +='<p>If you are ' + String(cutoff_percent) + '% correct or better during one round of a task, then you will pass that round and earn bonus points for it. If you are less than ' + String(cutoff_percent) + '% correct, you will not earn any points for that round.</p>'
          experiment_instructs +='<p>Keep in mind that you have to maintain this level of performance or better for up to ' + String(n_repetitions) + ' rounds.</p>'
          experiment_instructs +='<p>You can earn 1 to 5 points each round.</p>';
          experiment_instructs +='<p>At the end of the experiment, the points you earned will be added up and turned into a bonus payment. $$ The more points you earn, the more money you will be paid. $$</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
            }
          }]
        };*/

  screen3 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = "<p>Please read these instructions carefully, as you will be tested on your understanding of them.</p>"
          experiment_instructs +="<p>Before each round begins, you will be asked how many points is a fair wage to be at least " + String(cutoff_percent) + "% accurate at one of the tasks you just practiced.</p>"
          experiment_instructs +='<p>We want you to tell us truthfully how many points is a fair wage, and so we are using a payment process where <strong>telling the truth is the best strategy.</strong></p>'
          //experiment_instructs +='<p>In other words, you will be <strong>selling</strong> your task performance to us. </p>'
          experiment_instructs +='<p>This means that <strong>you have some control over what you will be paid</strong> at the end of the experiment.</p>'
          experiment_instructs +='<p>Whenever given a choice, <strong> the best thing you can do is to always ask for what you consider to be a fair wage for at least '+ String(cutoff_percent) + '% accuracy on the task.</strong></p>'
          experiment_instructs += space_bar_message;
           //return "<p style='font-size:25px'>" + n_back_instructs + " </p>";
          return experiment_instructs;
            }
          }]
        };

  screen4 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>What do we mean by "fair wage"?</p>'
          experiment_instructs +='<p>When you get hired to complete some work, you may ask for more or less payment based on how much you would enjoy doing that work.</p>'
          experiment_instructs +='<p>For example, you might dislike mowing the lawn more than watering the garden, and may ask for more in return for mowing the lawn.</p>'
          experiment_instructs +='<p>A fair wage for you is the amount of points you think it would be fair for you to receive in exchange for your completion of the task.</p>'
          experiment_instructs += '<p>While you won\'t get paid exactly what you ask for, the best strategy is to ask for what you truly believe is a fair wage for you.</p>'
          experiment_instructs += '<p>We\'re not going to get into the exact reasons why until the end of the experiment.</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
            }
          }]
        };      

  /*screen5 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>Why is the best strategy to ask for <strong> what you truly believe is a fair wage for you</strong>?</p>'
          experiment_instructs +='<p>Because the auction has one special rule.</p>'
          experiment_instructs +='<p>Its an unusual rule, but its implications are simple.</p>'
          experiment_instructs +='<p>There is no way of gaming the auction.</p>'
          experiment_instructs +='<p>The <strong> best </strong> thing you can do in every trial is <strong> ask yourself how many points you want to be paid if you achieve at least ' + String(cutoff_percent) + '% accuracy that round, and then ask for that amount.</strong></p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
            }
          }]
        };*/

  screen5 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function(){
        experiment_instructs ='<p>In short, once you report your fair wage, the computer will present you with a <strong> random offer.</strong></p>'
        experiment_instructs +='<p>If you ask for fewer points than the computer offers you, then you will complete the task and be given the points the computer offered.</p>'
        experiment_instructs +='<p>Say you ask for 3 points and the computer offers 4, then if you pass the task with at least ' + String(cutoff_percent) + '% accuracy, you will be given 4 points.</p>'
        experiment_instructs +='<p>If you ask for more points than the computer offers you, then you will do a different task on that round. The alternate task is worth <strong>1 point.</strong></p>'
        experiment_instructs +='<p>For example, if you ask for 5 points and the computer offers 4 points, then you will do something else and be given 1 point.</p>'
        experiment_instructs +='<p>There is no way to game the payment system. You will either get the computer\'s random point offer, or 1 point for an alternate task.</p>'
        experiment_instructs += space_bar_message;
        return experiment_instructs;
      }
    }]
  };

  screen6 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>What is the rule of the auction? Once you tell us what you think is a fair wage for you, the computer will present you with a <strong> random counter-offer. </strong></p>'
          experiment_instructs +='<p>If you ask for fewer points than the computer offers you, then you will complete the task and be given the points the computer offered.</p>'
          experiment_instructs +='<p>Say you ask for 3 points and the computer offers 4, then if you complete the task with at least '+ String(cutoff_percent) + '% accuracy, you will be given 4 points.</p>'
          experiment_instructs +='<p>If you ask for more points than the computer offers you, then you will do a different task on that round. The alternate task is worth <strong> 1 point. </strong></p>'
          experiment_instructs +='<p>For example, if you ask for 5 points and the computer offers 4 points, then you will do something else and be given 1 point.</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
        }
    }]
  };

  screen7 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>Why is it in your best interest to ask for a truly fair wage?</p>'
          experiment_instructs +='<p>You might think your best strategy is to always ask for a lot of points. This is <strong>incorrect.</strong></p>'
          experiment_instructs +='<p>The points you earn are determined by the computer\'s random offer, and <strong> not </strong> just what you ask for.</p>'
          experiment_instructs +='<p>Thus, if you ask for more than your truly fair wage, you will not affect your final point payment. Instead, you might lose the opportunity to do that task and earn more than 1 point.</p>'
          experiment_instructs += '<p>This happens if the computer\'s random number lies between your true fair wage and the value you reported.</p>'
          experiment_instructs +='<p>It follows that by asking for a fair wage for you, you will earn that amount or more if the computer is willing to pay as much as your fair wage.</p>'
          experiment_instructs +='<p>What happens if you ask for less than a truly fair wage? Then you will earn fewer points than you believe you deserve!</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
        }
    }]
  };

  exampleBDM = {
    timeline: exampleBDM()
  };

  screen8 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p><strong> Your best strategy is to always ask for what you believe is a fair wage for hitting at least ' + String(cutoff_percent) + '% accuracy on 1 round of a task.</strong></p>'
          experiment_instructs +='<p>We will now quiz you on the instructions you\'ve just read. If you fail this quiz, you cannot participate in this experiment.</p>'
          experiment_instructs +='<p>Press R if you\'d like to read these instructions again.</p>'
          experiment_instructs +='<p>Press the space bar if you\'d like to move on to the quiz.</p>'
          return experiment_instructs;
            }
          }]
        };

  screen9 = {
    timeline: [screen1, screen3, screen4, screen5, screen8],
    conditional_function: function(){
      lasttrialdata = jsPsych.data.getLastTrialData();
      key = lasttrialdata.select("key_press").values;
      if(key == jsPsych.pluginAPI.convertKeyCharacterToKeyCode('r')){
        return true;
      }
      if(key == jsPsych.pluginAPI.convertKeyCharacterToKeyCode(' ')){
        return false;
      }
    }
  }

  questions = [];
  response_prompt = "<p'font-size':25px>[Press A or B to indicate your answer.]</p>";
  quiz_questions = ["<p'font-size':25px>Imagine you think a true fair wage is 2 points. Which amount should you ask for?</p><p 'font-size':25px> A. 4</p><p 'font-size':25px> B. 2</p>", 
  "<p><img src='./static/img/offer_big.png'></p><p'font-size':25px>What happens when you ask for <strong>more</strong> points than the computer offers?</p><p>A. You will start the pictured task, and get the points the computer offered if you pass that round.</p><p 'font-size':25px> B. You will start a different task and be given 1 point if you pass that round.</p>",
  "<p'font-size':25px>Under which of the following two conditions will you do the pictured task?</p><p>A. When you ask for 2 points, and the computer randomly offers 3 points.</p><p 'font-size':25px> B. When you ask for 2 points, and the computer randomly offers 1 point.</p>",
  "<p><img src='./static/img/offer_small.png'></p><p'font-size':25px>What happens when you ask for <strong>fewer</strong> points than the computer offers?</p><p>A. You will start the task, and get the points the computer offered if you pass that round.</p><p 'font-size':25px> B. You will start a different task and be given 1 point if you pass that round.</p>",
  "<p'font-size':25px>Under which of the following two conditions will you <strong>not</strong> do the pictured task?</p><p>A. When you ask for 3 points, and the computer randomly offers 2 points.</p><p 'font-size':25px> B. When you ask for 3 points, and the computer randomly offers 4 points.</p>",
  "<p'font-size':25px>Which is the correct strategy for asking for points?</p><p>A. Ask for as many points as possible.</p><p 'font-size':25px> B. Tell the truth about what is a fair wage for you.</p>"];
  for(i=0;i<quiz_questions.length;i++){
    questions[i] = "<p'font-size':15px>[Question " + String(i+1) + " of " + String(quiz_questions.length) + "]</p>" + response_prompt + quiz_questions[i];
  }
  quiz_variables = [{stimulus: questions[0], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('b'),data:{task: "quiz"}},
  {stimulus: questions[1], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('b'),data:{task: "quiz"}},
  {stimulus: questions[2], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('a'),data:{task: "quiz"}},
  {stimulus: questions[3], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('a'),data:{task: "quiz"}},
  {stimulus: questions[4], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('a'),data:{task: "quiz"}},
  {stimulus: questions[5], key_answer: jsPsych.pluginAPI.convertKeyCharacterToKeyCode('b'),data:{task: "quiz"}}];
  all_images.push(["./static/img/offer_small.png","./static/img/offer_big.png"]);

  screen10 = {
    timeline: [{
      type: "categorize-html",
      stimulus: jsPsych.timelineVariable("stimulus"),
      key_answer: jsPsych.timelineVariable("key_answer"),
      data: jsPsych.timelineVariable("data"),
      show_stim_with_feedback: false,
      choices: [jsPsych.pluginAPI.convertKeyCharacterToKeyCode('a'), jsPsych.pluginAPI.convertKeyCharacterToKeyCode('b')]
    }],
    timeline_variables: quiz_variables,
    //randomize_order: true,
    data: {task: 'BDMquiz'}
  };

  screen11 = { 
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function(){
        accuracy = gradeBDMquiz();
        if(accuracy<=bdmcutoff){
          return "<p>Sorry, but you did not pass the instructions quiz. This experiment will now end and you will be compensated for the time you spent on it.</p>" + space_bar_message;
        } else {
          return "<p>Great job! We will now begin the real experiment.</p>" + space_bar_message;
        };
      },
      data: function(){
        accuracy = gradeBDMquiz();
        return {BDMquizgrade: accuracy, task: 'BDMquiz'}
      }
    }]
  };


return timeline = [screen3, screen4, screen5, exampleBDM, screen8, screen9, screen10, screen11];
//return timeline = [screen8, screen9];
};

function gradeBDMquiz(){
  BDM_quiz_grade = jsPsych.data.get().filter({task: "quiz"}).select("correct").values;
  accuracy = turnLogicalIntoAccuracy(BDM_quiz_grade)*100;
  return accuracy;
};


function moreBDMexplanation(){
  screen5 = {
  timeline: [{
    type: "html-keyboard-response",
      stimulus: function() {
      experiment_instructs = '<p>Why is the best strategy to ask for <strong> what you truly believe is a fair wage for you</strong>?</p>'
        experiment_instructs +='<p>Because the payment process of this experiment has one special rule.</p>'
        experiment_instructs +='<p>Its an unusual rule, but its implications are simple.</p>'
        experiment_instructs +='<p>There is no way of gaming the process.</p>'
        experiment_instructs +='<p>The <strong> best </strong> thing you can do in every trial is <strong> ask yourself how many points you want to be paid if you achieve at least ' + String(cutoff_percent) + '% accuracy that round, and then ask for that amount.</strong></p>'
        experiment_instructs += space_bar_message;
        return experiment_instructs;
      }
    }]
  };

  screen6 = {
  timeline: [{
    type: "html-keyboard-response",
      stimulus: function() {
        experiment_instructs = '<p>What is the rule? Once you tell us what you think is a fair wage for you, the computer will present you with a <strong> random counter-offer. </strong></p>'
        experiment_instructs +='<p>If you ask for fewer points than the computer offers you, then you will complete the task and be given the points the computer offered.</p>'
        experiment_instructs +='<p>For example, say you ask for 3 points and the computer offers 4, then if you complete the task with at least '+ String(cutoff_percent) + '% accuracy, you will be given 4 points.</p>'
        experiment_instructs +='<p>If you ask for more points than the computer offers you, then you will do a different task on that round. The alternate task is worth <strong> 1 point. </strong></p>'
        experiment_instructs +='<p>For example, if you ask for 5 points and the computer offers 4 points, then you will do something else and be given 1 point.</p>'
        experiment_instructs += space_bar_message;
        return experiment_instructs;
      }
    }]
  };

  screen7 = {
    timeline: [{
      type: "html-keyboard-response",
        stimulus: function() {
          experiment_instructs = '<p>Why is it in your best interest to ask for a truly fair wage?</p>'
          experiment_instructs +='<p>You might think your best strategy is to always ask for a lot of points. This is <strong>not the case.</strong></p>'
          experiment_instructs +='<p>The points you earn are determined by the computer\'s random offer, and <strong> not </strong> just what you ask for.</p>'
          experiment_instructs +='<p>Thus, if you ask for more than your truly fair wage, you will not affect your final point payment. Instead, you might lose the opportunity to do that task and earn more than 1 point.</p>'
          experiment_instructs +='<p>This happens if the computer\'s random number lies between your true fair wage and the value you reported.</p>'
          experiment_instructs +='<p>It follows that by asking for a fair wage for you, you will earn that amount or more if the computer offers to pay as much as your fair wage.</p>'
          experiment_instructs +='<p>What happens if you ask for less than a truly fair wage? Then you will earn fewer points than you believe you deserve!</p>'
          experiment_instructs += space_bar_message;
          return experiment_instructs;
            }
          }]
        };

  screen12 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: "<p>We hope this explanation has cleared things up.</p><p>Press the space bar to continue.</p>"
    }]
  }

  return timeline = [screen5, screen6, screen7, screen12];
}

function exampleBDM(){
  var BDM_trial = [{stimulus: fractals[2], data: {task: 'instructs'}}];//initialize BDM trial parameters
  var BDM1 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: jsPsych.timelineVariable('stimulus'),
      data: jsPsych.timelineVariable('data'),
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: 0,
      step: 4,
      trial_duration: instructs_timing,
      prompt: "<p>This is how the payment process will work.</p><p>You will report your fair wage with the slider on the right.</p><p>You have thirty seconds to respond.</p><p>Try moving the slider up just a little bit and see what happens.</p>"
      }],
    timeline_variables: BDM_trial
  };

  var BDM_outcome1 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function(data){
        offer = 4.88;
        var lasttrialdata = jsPsych.data.getLastTrialData();
        points = lasttrialdata.select("response").values;
        points = 1+points*0.04;
        points = Number.parseFloat(points).toFixed(2);
        message = "<p>You asked for " + String(points) + " points.</p><p>The computer offered " + String(offer) + " points.</p><p>Because the computer's offer is larger, you would move onto the pictured task and, if you achieved " + String(cutoff_percent) + "% accuracy on this task, you would be given " + String(offer) + " points at the end.</p><p>Press the space bar to continue.</p>"
        return message;      
      },
      trial_duration: instructs_timing
    }]
  }

  var BDM2 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: jsPsych.timelineVariable('stimulus'),
      data: jsPsych.timelineVariable('data'),
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: 0,
      step: 4,
      trial_duration: instructs_timing,
      prompt: "<p>Now, try moving the slider up a lot and see what happens.</p><p>You have thirty seconds to respond.</p>"
      }],
    timeline_variables: BDM_trial
  };

  var BDM_outcome2 = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: function(data){
        offer = 1.24;
        var lasttrialdata = jsPsych.data.getLastTrialData();
        points = lasttrialdata.select("response").values;
        points = 1+points*0.04;
        points = Number.parseFloat(points).toFixed(2);
        message = "<p>You asked for " + String(points) + " points.</p><p>The computer offered " + String(offer) + " points.</p><p>Because the computer's offer is smaller, you would skip the pictured task and, if you achieved " + String(cutoff_percent) + "% accuracy on an alternate task, you would be given 1 point at the end.</p><p>Press the space bar to continue.</p>"
        return message;      
      },
      trial_duration: instructs_timing
    }]
  }

return timeline = [BDM1,BDM_outcome1,BDM2,BDM_outcome2];

}