// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function create_difficulty_ratings(){
  var debrief = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: "<p>The next few screens have some survey questions on them.</p><p>Press respond to every question you are comfortable answering.</p><p>Press any key to continue.</p>"
    }]
  }

  var rating = []; var rating_timeline = [];
  uniqueValues.push(easy_task_ID); //list tasks used in this version of the exp

  rating_timeline.push(debrief);

  var srt_point = Math.random()*100;
  //var task_ID = uniqueValues[i];
  var rating1 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: stim = fractals[0]+'" width="'+fractal_size+'" height ="'+fractal_size+'"',
      data:{task: 'difficultyrating', task_displayed: 0, slider_start: srt_point},
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: srt_point,
      step: 4,
      prompt: "<p>How difficult was this task for you (where 1 is very easy, and 5 is very hard)?</p>",
      trial_duration: 60000 //60 seconds to respond
    }]
  };
  var srt_point = Math.random()*100;
  var rating2 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: stim = fractals[1]+'" width="'+fractal_size+'" height ="'+fractal_size+'"',
      data:{task: 'difficultyrating', task_displayed: 1, slider_start: srt_point},
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: srt_point,
      step: 4,
      prompt: "<p>How difficult was this task for you (where 1 is very easy, and 5 is very hard)?</p>",
      trial_duration: 60000 //60 seconds to respond
    }]
  };
  var srt_point = Math.random()*100;
  var rating3 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: stim = fractals[2]+'" width="'+fractal_size+'" height ="'+fractal_size+'"',
      data:{task: 'difficultyrating', task_displayed: 2, slider_start: srt_point},
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: srt_point,
      step: 4,
      prompt: "<p>How difficult was this task for you (where 1 is very easy, and 5 is very hard)?</p>",
      trial_duration: 60000 //60 seconds to respond
    }] 
  };
  var rating4 = {
    timeline: [{
      type: "html-slider-response",
      stimulus: stim = fractals[7]+'" width="'+fractal_size+'" height ="'+fractal_size+'"',
      data:{task: 'difficultyrating', task_displayed: 7, slider_start: srt_point},
      require_movement: true,
      labels: ['5','4','3','2','1'],
      start: srt_point,
      step: 4,
      prompt: "<p>How difficult was this task for you (where 1 is very easy, and 5 is very hard)?</p>",
      trial_duration: 60000 //60 seconds to respond
    }] 
  };

  rating_timeline.push(rating1);
  rating_timeline.push(rating2);
  rating_timeline.push(rating3);
  rating_timeline.push(rating4);

  return rating_timeline;
} //end of function

function create_need_for_cognition(){

  space_bar_message = "<p>[Press the space bar to continue.]</p>";

  NFC_instructs = "<p><strong>For each of the statements below, please indicate whether or not the statement is characteristic of you or of what you believe.</strong></p>"

  var options = ["extremely uncharacteristic of me", "somewhat uncharacteristic of me", "uncertain", "somewhat characteristic of me", "extremely characteristic of me"];

  //var standards = {options: options, required:false, horizontal:true};
  var question_list = ["1. I prefer complex to simple problems.",
    "2. I like to have the responsibility of handling a situation that requires a lot of thinking.",
    "3. Thinking is not my idea of fun.",
    "4. I would rather do something that requires little thought than something that is sure to challenge my thinking abilities.",
    "5. I try to anticipate and avoid situations where there is a likely chance I will have to think in depth about something.",
    "6. I find satisfaction in deliberating hard and for long hours.",
    "7. I only think as hard as I have to.",
    "8. I prefer to think about small daily projects to long term ones.",
    "9. I like tasks that require little thought once I've learned them.",
    "10. The idea of relying on thought to make my way to the top appeals to me.",
    "11. I really enjoy a task that involves coming up with new solutions to problems.",
    "12. Learning new ways to think doesn't excite me very much.",
    "13. I prefer my life to be filled with puzzles I must solve.",
    "14. The notion of thinking abstractly is appealing to me.",
    "15. I would prefer a task that is intellectual, difficult, and important to one that is somewhat important but does not require much thought.",
    "16. I feel relief rather than satisfaction after completing a task that requires a lot of mental effort.",
    "17. It's enough for me that something gets the job done; I don't care how or why it works.",
    "18. I usually end up deliberating about issues even when they do not affect me personally."]

    // make a list of dictionaries [{},{},{}] with standard settings

  npages = 3;  
  ends = [question_list.length/npages, 2*(question_list.length/npages), question_list.length];
  formattedqs1 = [];
  for(var q=0; q<ends[0]; q++){
    formattedqs1.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true, name:"NFC"+question_list[q][0]+question_list[q][1]});
  }
  //add a Screener question
  formattedqs1.push({prompt: "<strong>7. Please select 'extremely characteristic of me' if you've read this question.</strong>", options: options, required:false, horizontal:true,name:'screener'})
  formattedqs2 = [];
  for(var q=ends[0]; q<ends[1]; q++){
    formattedqs2.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true, name:"NFC"+question_list[q][0]+question_list[q][1]});
  }
  formattedqs3 = [];
  for(var q=ends[1]; q<ends[2]; q++){
    formattedqs3.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true, name:"NFC"+question_list[q][0]+question_list[q][1]});
  }
    // standards['prompt'] = "This is a question?"

  var NFC1 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs1,
      preamble: NFC_instructs,
      data: {task: 'NFC'}
    }]
  }
  var NFC2 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs2,
      preamble: NFC_instructs,
      data: {task: 'NFC'}
    }]
  }
  var NFC3 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs3,
      preamble: NFC_instructs,
      data: {task: 'NFC'}
    }]
  }
    //NFC.push(multi_choice_block_horizontal);

  return timeline = [NFC1, NFC2, NFC3];

};

function create_demographics(){

  space_bar_message = "<p>[Press the space bar to continue.]</p>";

  //var standards = {options: options, required:false, horizontal:true};

  var question_list = [{prompt: "How much education do you have?", options: ["Some elementary or middle school","Some high school","Graduated high school","Some college","Associate's degree","Bachelor's degree","Started higher education (MA, PhD, MD, etc)","Completed higher education (MA, PhD, MD, etc.)"],name: "education"},
  {prompt: "What sex are you?", options: ["Male","Female","Other","Prefer not to state"],name: "sex"},
  {prompt: "Are you colorblind?", options: ["Yes","No","I don\'t know","Prefer not to state"], name: "colorblindness"},
  {prompt: "How hard was this task for you?", options: ["Very Easy", "Easy","Okay","Hard","Very Hard"], name: "difficulty"},
  {prompt: "Are you right- or left-handed?", options: ["Right-handed","Left-handed","Ambidextrous"],name: "handedness"}];

  var all_questions = [];
  all_questions = [];
  for(var q=0; q<question_list.length; q++){
    all_questions.push({prompt: "<strong>" + question_list[q]['prompt'] + "</strong>", options: question_list[q]['options'], required:false, horizontal:true});
  }
    // standards['prompt'] = "This is a question?"

  var demo = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: all_questions,
      data: {task: 'demographics'}
    }]
  }

  free_response_questions = [{prompt: "How old are you?",rows: 1,name: "age"},
  {prompt: "What was your strategy for determining your fair wage for each task? Did it change over time?", rows: 5, name: "fair wage"},
  {prompt: "Do you have any other comments about this task? For example, was it fun to do? Did you encounter any problems with it?", rows: 5, name: "other"}]

  var textboxes = {
    timeline: [{
      type: 'survey-text',
      questions: free_response_questions,
      data: {task: 'free_response'}
    }]
  }

  /*var full_timeline = {
    timeline: [demo, textboxes]
  }*/
  return timeline = [demo, textboxes];

};

function create_SAPS(){
  var SAPS_instructions = '<p><strong>The following items are designed to measure certain attitudes people have toward themselves, their performance, and toward others. It is important that your answers be true and accurate for you. For each statement, please select an answer from strongly disagree to strongly agree to describe your degree of agreement with each item.</p></strong>';

  var options = ["strongly disagree", "disagree", "somewhat disagree", "neutral", "somewhat agree", "agree", "strongly agree"];

  //var standards = {options: options, required:false, horizontal:true};
  var question_list = ["1. I have high expectations for myself.",
    "2. Doing my best never seems to be enough.",
    "3. I set very high standards for myself.",
    "4. I often feel disappointment after completing a task because I know I could have done better.",
    "5. I have a strong need to strive for excellence.",
    //screener
    "6. My performance rarely measures up to my standards.",
    "7. I expect the best from myself.",
    "8. I am hardly ever satisfied with my performance."]

    // make a list of dictionaries [{},{},{}] with standard settings

  ends = [5, 8];
  formattedqs1 = [];
  for(var q=0; q<ends[0]; q++){
    formattedqs1.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true, name:"SAPS"+question_list[q][0]+question_list[q][1]});
  }
  formattedqs2 = [];
  formattedqs2.push({prompt: "<strong>6. Please select 'somewhat agree' if you've read this question.</strong>", options: options, required:false, horizontal:true, name:'screener'})
  for(var q=ends[0]; q<ends[1]; q++){
    formattedqs2.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true, name:"SAPS"+question_list[q][0]+question_list[q][1]});
  }
  // standards['prompt'] = "This is a question?"

  var SAPS1 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs1,
      preamble: SAPS_instructions,
      data: {task: 'SAPS'}
    }]
  }
  var SAPS2 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs2,
      preamble: SAPS_instructions,
      data: {task: 'SAPS'}
    }]
  }
    //NFC.push(multi_choice_block_horizontal);

  return timeline = [SAPS1, SAPS2];

}