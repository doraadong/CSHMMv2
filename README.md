# Introduction

This is an modified CSHMM version, which is based on the CSHMM tool and the web visualization tool initially developed by [Chieh Lin](https://github.com/jessica1338) & [Jun Ding](https://github.com/phoenixding) and is associated with the published work [[1]](#1) (refered to as Lin & Ding version). 

# What changed compared to Lin & Ding version
1. Enable initialize the model using existing cell type labels 
2. In the web visualization, 
  * enable displaying the expression level of a TF's targets per path  
  * for each DE gene, enable displaying the TFs regulating it in the DE results 
  * change the coloring rules to color all cells with darker color indicating higher expression. 
  * enable displaying average expression per path 
  * enable downloading TFs and eTFs. 

# Installation
See the original CSHMM for setting up the environments. 

# Example usages
1. generate initialization file by using [run_init.q](run_init.q)
2. run CSHMM using [run_cshmm.q](run_cshmm.q). 
3. for visualization
  * generate .json files for visualization by running [run_json.q](run_json.q). 
  * copy .json files to the web_visualization folder and rename them as data.json and CellViz.json. 
  * generate .db (database) by runing [run_create_table.sh](web_visualization/utils/run_create_table.sh) and move it to ht_bin folder.
  * boot server by command: python simpleServer.py

# To-do list
1. Update the readme.pdf file for the web visualization tool. 

# References
<a id="1">[1]</a> 
Hurley, Killian, Jun Ding, Carlos Villacorta-Martin, Michael J. Herriges, Anjali Jacob, Marall Vedaie, Konstantinos D. Alysandratos, et al. 2020. “Reconstructed Single-Cell Fate Trajectories Define Lineage Plasticity Windows during Differentiation of Human PSC-Derived Distal Lung Progenitors.” Cell Stem Cell 26 (4): 593–608.e8.
