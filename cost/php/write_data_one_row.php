<?php
//josh de leeuw wrote this
// this path should point to your configuration file.
include('db_config_trial.php');

$data_array = json_decode(file_get_contents('php://input'), true);

try {
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
  //echo count($col_names);
  // Second stage is to create prepared SQL statement using the column
  // names as a guide to what values might be in the JSON.
  // If a value is missing from a particular trial, then NULL is inserted
  $sql = "INSERT INTO $table VALUES(";
  for($i = 0; $i < count($col_names); $i++){
    $name = $col_names[$i];
    $sql .= ":$name";
    if($i != count($col_names)-1){
      $sql .= ", ";
    }
  }
  $sql .= ");";
  $insertstmt = $conn->prepare($sql);
  //for($i=0; $i < count($data_array); $i++){
  //$i = (count($data_array))-1; //php bug here //added extra parentheses to see whether this helps
  $i = (!empty($data_array) ? count($data_array)-1 : 0);
    for($j = 0; $j < count($col_names); $j++){
      $colname = $col_names[$j];
      if(!isset($data_array[$i][$colname])){
        $insertstmt->bindValue(":$colname", null, PDO::PARAM_NULL);
      } else {
        $insertstmt->bindValue(":$colname", $data_array[$i][$colname]);
        //echo 'insert statement: ' . $insertstmt;
      }
    }
    $insertstmt->execute(); //josh
    //$result = mysqli_stmt_execute($insertstmt); //sarah
  //}
  //echo '{"success": true}';
} catch(PDOException $e) {
  echo '{"success": false, "message": ' . $e->getMessage() . '}. ';
  //echo ' spit out : ' . $data_array[50][$col_names[40]];
  echo 'i is' . $i ;
  echo ' count output is ' . count($data_array) ;
  // to try: go row by row
}
$conn = null;
?>