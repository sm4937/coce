// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function create_need_for_cognition(){

  space_bar_message = "<p>[Press the space bar to continue.]</p>";

  NFC_instructs = "<p><strong>For each of the statements below, please indicate whether or not the statement is characteristic of you or of what you believe.</strong></p>"

  var debrief = {
    timeline: [{
      type: "html-keyboard-response",
      stimulus: "<p>The next few screens have some survey questions on them.</p><p>Press respond to every question you are comfortable answering.</p>" + space_bar_message
    }]
  }

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
    formattedqs1.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true});
  }
  formattedqs2 = [];
  for(var q=ends[0]; q<ends[1]; q++){
    formattedqs2.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true});
  }
  formattedqs3 = [];
  for(var q=ends[1]; q<ends[2]; q++){
    formattedqs3.push({prompt: "<strong>" + question_list[q] + "</strong>", options: options, required:false, horizontal:true});
  }
    // standards['prompt'] = "This is a question?"

  var NFC1 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs1,
      preamble: NFC_instructs
    }]
  }
  var NFC2 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs2,
      preamble: NFC_instructs
    }]
  }
  var NFC3 = {
    timeline: [{
      type: 'survey-multi-choice',
      questions: formattedqs3,
      preamble: NFC_instructs
    }]
  }
    //NFC.push(multi_choice_block_horizontal);

  return timeline = [debrief, NFC1, NFC2, NFC3];

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
      questions: all_questions
    }]
  }

  free_response_questions = [{prompt: "How old are you?",rows: 1,name: "age"},
  {prompt: "What was your strategy for determining your fair wage for each task? Did it change over time?", rows: 5, name: "fair wage"},
  {prompt: "Do you have any other comments about this task? For example, was it fun to do? Did you encounter any problems with it?", rows: 5, name: "other"}]

  var textboxes = {
    timeline: [{
      type: 'survey-text',
      questions: free_response_questions
    }]
  }

  /*var full_timeline = {
    timeline: [demo, textboxes]
  }*/
  return timeline = [demo, textboxes];

};