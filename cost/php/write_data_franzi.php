<?php

//list every variable in my data, which I hate but okay
$rt = $_POST['rt'];
$responses = $_POST['responses'];
$question_order = $_POST['question_order'];
$trial_type = $_POST['trial_type'];
$trial_index = $_POST['trial_index'];
$time_elapsed = $_POST['time_elapsed'];
$internal_node_id = $_POST['internal_node_id'];
$TOT = $_POST['TOT'];
$overall =  $_POST['overall'];
$performance_by_block = $_POST['performance_by_block'];
$value_list = $_POST['value_list'];
$offer_list = $_POST['offer_list'];
$exp_version = $_POST['exp_version'];
$n_switch_rule_sides = $_POST['n_switch_rule_sides'];
$subjnum = $_POST['subjnum'];
$worker_ID = $_POST['worker_ID'];
$assignment_ID = $_POST['assignment_ID'];
$hit_ID = $_POST['hit_ID'];
$success = $_POST['success'];
$stimulus = $_POST['stimulus'];
$key_press = $_POST['key_press'];
$task = $_POST['task'];
$correct = $_POST['correct'];
$detect =  $_POST['detect'];
$correct_key = $_POST['correct_key'];
$tasknum = $_POST['tasknum'];
$practice_accuracy = $_POST['practice_accuracy'];
$practice = $_POST['practice'];
$nback = $_POST['nback'];
$n = $_POST['n'];
$number_practice_hard = $_POST['number_practice_hard'];
$response = $_POST['response'];
$BQMquizgrade = $_POST['BQMquizgrade'];
$task_displayed = $_POST['task_displayed'];
$slider_start = $_POST['slider_start'];
$stimnum = $_POST['stimnum'];
$interaction_data = $_POST['interaction_data'];
echo 'list created';

$path = "db_config.php";
require_once $path;

$config = new DatabaseConfiguration();
echo 'config created';
$conn = $config->createConnection();
// Insert trial entry into trial table
// Protect against injection using prepared statements and parameterized queries
$stmt = mysqli_prepare($conn, "INSERT INTO cost1 (rt, responses, question_order, trial_type, trial_index, time_elapsed, internal_node_id​, TOT, overall, performance_by_block, points_list, value_list, offer_list, exp_version, n_switch_rule_sides, subjnum, worker_ID, assignment_ID, hit_ID, success, stimulus, key_press, task, correct, detect, correct_key, tasknum, practice_accuracy, practice, nback, n, number_practice_hard, response, BQMquizgrade, task_displayed, slider_start, stimnum, interaction_data) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
mysqli_stmt_bind_param($stmt, 'ssssssssssssssssssssssssssssssssssssss', $rt, $responses, $question_order, $trial_type, $trial_index, $time_elapsed, $internal_node_id​, $TOT, $overall, $performance_by_block, $points_list, $value_list, $offer_list, $exp_version, $n_switch_rule_sides, $subjnum, $worker_ID, $assignment_ID, $hit_ID, $success, $stimulus, $key_press, $task, $correct, $detect, $correct_key, $tasknum, $practice_accuracy, $practice, $nback, $n, $number_practice_hard, $response, $BQMquizgrade, $task_displayed, $slider_start, $stimnum, $interaction_data);
$result = mysqli_stmt_execute($stmt);

if ($result === TRUE) {
    echo "Recorded new trial entry \n";
} else {
    echo "Error: " . $conn->error . "\n";
}
mysqli_stmt_close($stmt);
mysqli_close($conn);
?>

        function saveData_franzi() {
          $.post("php/write_data_franzi.php", dataToPost).fail( function () {
          console.log("Failed to create new assignment entry in database.");
          });
        }