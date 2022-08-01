#!/bin/bash

nameData=Data

useConfig=True
label2cluster=1
nPCA=1000
ancestor=0


if [ "$useConfig" == "True" ]
then 
	if [ "$label2cluster" -eq "1" ]
	then
		output_file="$nameData"_init_label.txt
	else
		output_file="$nameData"_init_"$nPCA".txt
	fi
else
	output_file="$nameData"_init_auto_"$nPCA".txt
fi

data="$nameData".txt

echo Running for "$nameData" and use config: "$useConfig", use known label for cluster: "$label2cluster", number of PCA: "$nPCA", larget type: True, ancestor: "$ancestor"


if [ "$useConfig" == "True" ]
then
	config="$nameData".config

	python scdiff_init.py -tf input/tool/tfDNA_predicted_100.txt.update -d input/$data -o output/$output_file -l 1 -c input/$config -lc $label2cluster -np "$nPCA" -a "$ancestor"
else
	python scdiff_init.py -tf input/tool/tfDNA_predicted_100.txt.update -d input/$data -o output/$output_file -l 1 -np "$nPCA" -a "$ancestor"
fi

