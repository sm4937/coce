READ ME for Costs of Cognitive Effort Task Code

Within this folder is all the code for running versions 1 through 4 of the
Costs of Cognitive Effort task on Amazon Mechanical Turk. There is also now 
obsolete code for saving data on the Kyb Lab (AGPD) server at the Max Planck 
Institute for Biological Cybernetics. 

The code is written in Javascript, with some HTML & CSS elements. The code
leverages a package for building psychology experiments to run
online, called jsPsych. jsPsych is a timeline creation package, which
pre-loads all the events of an experiment before the experiment even 
begins. It also provides structure for creating the elements of this timeline
through jsPsych plugins. The plugins are essentially templates for what
each trial of a psych experiment could look like - plugins which display
images and catalog responses to those images, for example. These plugins
are the strength of jsPsych, in that someone who doesn't know anything about 
front-end website design doesn't have to mess around with CSS unless they
really want to. The weakness of jsPsych is the timeline, which makes the 
implementation of jsPsych experiments fairly rigid. In Index.html, 
main_timeline_creation.js, and other files in this folder, you'll see the 
workarounds to the rigidity of jsPsych that I implemented, including many
conditional_function calls. main_timeline_creation.js is the cleanest example
of this.

Index.html is the central hub of the task code, where all code is brought
together to run the task. At the beginning of Index.html is a line of
code where you can switch between task versions:

var exp_version = 4; 
      // 1 is the combine, detection version
      // 2 is the n-switch version
      // 3 is the n-back version
      // 4 is the n-back version with a 3-target detection task (ndetection)

Version 4 is the final version which we are attempting to publish the results 
of now. We piloted versions 1-4 and found that version 1 was too hard, and
version 2 was too easy. Version 3 is a simpler iteration of version 4, which
includes one more task, the 3-detection task.

There are four folders of code here:
- custom_plugins: Javascript scripts written by yours truly (Sarah Master),
which create timelines to display our questionnaires, run the different psych tasks,
run the BDM auction procedure, and create the practice blocks and quiz. This is
the experimental code, split up into many different files for readability, and 
assembled to run the experiment in Index.html.
- jsPsych: This is the jsPsych code base, not written by me. This is available
for anyone to download on the jsPsych website. Some of the base code
has been modified by me, including plugins/jspsych-html-slider-response.js, 
which I made aesthetic changes to in HTML. 
- php: This is the code I wrote (and some code borrowed from Franziska
Broeker) to save data to my SQL database. It's pretty much all obsolete,
given that data is saved to the server in a different way now, but it 
may be a useful template for someone to learn how to do data saving in
SQL using PHP scripts. Note that for this type of system to work it's
important you pre-create the data tables that you later call upon in db_config.php,
with all the right variable names & types.
- static: This is where I stored all the images used in the experiment, namely
the fractals used as task labels. These are not my fractal images and I'm probably
violating some laws somewhere by using them without the proper licensing but, oh
well. Storing the images in a separate "static" folder is a convention I picked
up when I learned to write MTurk tasks with Psiturk. You don't have to do it 
this way, but I found it easier to wrap my head around this structure.

The task code should work right out of the box. To run the entire experiment, just
open index.html on any browser (preferably Google Chrome, for formatting reasons).
To upload the code and run it online, move this entire folder onto the appropriate
server, or mount a Psiturk experiment inside it (or something similar), and it
should work. The only issue is data saving, which will need to be re-configured
to play nice with your server setup. jsPsych provides some template data-saving
codes, but there will be issues of read/write permissions, or SQL database
passwords, that you will need to figure out for your setup.

Have fun!


