#!/bin/bash

nameData=Data

num_genes=30000
suffix=_label

nameInit="$nameData"

data="$nameData".txt
init_file="$nameInit"_init"$suffix".txt
model_name="$nameData"_tf"$suffix"_"$num_genes"

echo "Use $init_file as initialization file, use $data as data file, output $model_name ; using genes: $num_genes" 


python CSHMM_train_modified_20180403.py \
--data_file input/$data \
-tf input/tool/tfDNA_predicted_100.txt.update \
--structure_file output/init/$init_file \
--n_split 100 -ng "$num_genes" \
--n_iteration 10 \
--cross_validation 0 \
--random_seed 5 \
-k 10 \
-opt genlasso \
--model_name model/"$model_name" \
