#!/bin/bash

nameData=Data

suffix=_label
iteration=10
num_genes=30000
use_time_as_label=False

data="$nameData".txt
model_name="$nameData"_tf"$suffix"_"$num_genes"_it"$iteration".pickle

echo $model_name, $data
python CSHMM_json.py input/$data model/"$model_name" input/tool/tfDNA_predicted_100.txt.update $num_genes 0.6 $use_time_as_label


