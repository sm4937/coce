<?php
//Franzi wrote this
$MTurkWorkerId = $_POST['MTurkWorkerId'];
$MTurkAssignmentId = $_POST['MTurkAssignmentId'];
$MTurkHitId = $_POST['MTurkHitId'];
$status = $_POST['status'];
$code = $_POST['code'];
$browser = $_POST['browser'];
$platform = $_POST['platform'];
$browserLanguage = $_POST['browserLanguage'];
$isMobile = $_POST['isMobile'];
$beginHit = $_POST['beginHit'];
$codeVersion = $_POST['codeVersion'];

$path = "db_config.php";
require_once $path;

$config = new DatabaseConfiguration();
$conn = $config->createConnection();

// Insert assignment entry into assignment table
// Protect against injection using prepared statements and parameterized queries
$stmt = mysqli_prepare($conn, "INSERT INTO assignment (MTurkWorkerId, MTurkAssignmentId, MTurkHitId, status, code, browser, platform, browserLanguage, isMobile, beginHit, codeVersion) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
mysqli_stmt_bind_param($stmt, 'sssssssssss', $MTurkWorkerId, $MTurkAssignmentId, $MTurkHitId, $status, $code, $browser, $platform, $browserLanguage, $isMobile, $beginHit, $codeVersion);
$result = mysqli_stmt_execute($stmt);
$new_assignment_id = mysqli_insert_id($conn);

//if ($result === TRUE) {
//    echo "New assignment entry successfully created.";
//} else {
//    echo "Error: " . $conn->error;
//}

echo $new_assignment_id;
mysqli_stmt_close($stmt);
mysqli_close($conn);
?>