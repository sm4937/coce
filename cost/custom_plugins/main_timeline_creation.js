// functions and such just for n-back creation 
// set up the n-back task parameters and stimulus list
function create_main_timeline(blocknums, task_list){

  var practice_flag = false; //set tasks to run normally again
  
  loop_timeline = []; //initialize empty loop

  for(var loopi=blocknums[0]; loopi <= blocknums[1]; loopi++){ // do the whole thing n times ??
    srt_point = Math.random()*100;
    points = 0; //updated every BDM loop
    task_ID = task_list[loopi];
    var BDM_trial = [{stimulus: fractals[task_ID], data: {task: 'BDM', task_displayed: task_ID, tasknum: loopi, slider_start: srt_point}}];//initialize BDM trial parameters
    var BDM = {
      timeline: [{
        type: "html-slider-response",
        stimulus: jsPsych.timelineVariable('stimulus'),
        data: jsPsych.timelineVariable('data'),
        require_movement: true,
        labels: ['5','4','3','2','1'],
        start: srt_point,
        step: 4,
        prompt: "<p>How many points is a fair wage for completing this task?</p>",
        trial_duration: instructs_timing //30 seconds to respond
        }],
      timeline_variables: BDM_trial
    };

    var BDM_outcome = {
      timeline: [{
        type: "html-keyboard-response",
        stimulus: function(data){
          getBDMresponse();
          BDM_message = runBDM(points,offer,task_ID);
          return BDM_message;
         },
         trial_duration: instructs_timing
       }]
    }
      
  // conditional function for hard vs. easy task timeline running
    var run_hard_task = {
      timeline: run_task_by_ID(task_ID,loopi),//setup_nBack(loopi),
      conditional_function: function(){
          if(points<=offer & points!=null){
            return true;
          } else {
            console.log('skipping hard task')
            return false;
          }
        }
      };

    var run_easy_task = {
      timeline: run_task_by_ID(task_ID+1,loopi),
      conditional_function: function(){
          if(points>offer || points==null){ //didn't do hard n-back
            points = 1;
            return true;
          } else {
            console.log('skipping easy task')
            return false;
          }
        }
      };

      loop_timeline.push(BDM);
      loop_timeline.push(BDM_outcome);
      loop_timeline.push(run_hard_task);
      loop_timeline.push(run_easy_task);

    } // outer trial (task) loop

  return loop_timeline;

};