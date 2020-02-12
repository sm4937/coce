<?php
// josh de leeuw, not working
$post_data = json_decode(file_get_contents('php://input'), true);

$name = "data/".$post_data['filename'].".csv";
$data = $post_data['filedata'];
// tried changing direction of /'s
file_put_contents($name, $data);
?>