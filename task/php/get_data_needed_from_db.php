<?php
//josh de leeuw wrote this
// this path should point to your configuration file.
include('db_config.php');

$conn = new PDO("mysql:host=$servername;port=$port;dbname=$dbname", $username, $password);
$conn->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);
// First stage is to get all column names from the table and store
// them in $col_names array.
$stmt = $conn->prepare("SHOW COLUMNS FROM `$table`");
$stmt->execute();
$col_names = array();
while($row = $stmt->fetchColumn()) {
  $col_names[] = $row;
}

for($i = 0; $i < count($col_names); $i++){
    $name = $col_names[$i];
    echo $name . '/';
}

return $col_names;
?>