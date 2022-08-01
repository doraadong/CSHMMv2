# __Cosimo__-<ins>Co</ins>ntinuous <ins>S</ins>tate H<ins>i</ins>dden <ins>M</ins>arkov M<ins>o</ins>del
![images/cosimo.png](images/cosimo.png)    

Cosimo is an interactively visualized Continuous hidden Markov model to reconstruct the cell differentiation trajectories.  
With the model,  users can generate many visualizations based on their specific needs. They can also make queries and predictions with their own interests. 
The model (right) can be configured by multiple panels (left), which were explained in details below. 


# Cosimo model 
![images/model.jpg](images/model.jpg)
 __cell__: a dot in the model, which is color-coded by the capture time of the cell.  
 __P0-P10__: paths in the model. Cells are differentiating into different cell fates via different paths.   
 __N0-N11__: Nodes in the model, which represent specific cell states (e.g. divergence states or end states).   
 __TFs/eTFs__: The regulators which dictate the cell differentiation in the path ending at that specific node.  
 The TFs are inferred based on their target genes.  The target genes of the TFs are differentially expressed in the regulated paths. 
 eTFs are predicted based on the expression of those TFs. To be more specific, the expression of the eTFs in the regulated paths are significantly different compared with the expression in the parent and sibling paths. 

# Global Config Panel

Set background, Node color, Text color, Line color, and Cell Size of the model

  

# Cell Plots Panel
![images/cellplot.png](images/cellplot.png)  
Get the t-SNE (https://lvdmaaten.github.io/tsne/), PCA (https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) and Diffusion map (https://academic.oup.com/bioinformatics/article/31/18/2989/241305) for all the cells in the dataset.  
You can also choose how would you want to color the cells using the "Color the cells by:" drop-down menu (capture time/Model path of the cell)

For example, t-sne plot and color-coded by the Model path (cells at different paths will be colored differently). 
![images/tsne.png](images/tsne.png)

# CELL Config Panel
![images/cellconfig.png](images/cellconfig.png)  
__(1) ColorCellsByGene__:   
![images/colorcell.png](images/colorcell.png)   
Color all the cells based on the expression level of input genes.  
If > the level specified by the slider below => red    
Otherwise => gray    
The cell size can also be set by the Global config panel. 

![images/colorcellres.png](images/colorcellres.png)   
You can use the slider to customize the threshold. 
The coloring process usually takes a few secs but it varies based on the connection. 

__(2) Cell Group analysis__:   
![images/cellgroups.png](images/cellgroups.png)  
	__(2.1) Set Group1 as selected:__    
	While SHIFT pressed, use mouse to drag select the cells on the right. The selected cells will turn into gray.     
	![images/dragselect.png](images/dragselect.png)    
	Then, press the "Set Group 1 as selected" to appoint the selected cells as group 1.     
	The appointed cells will be colored as the specified group1 color.    
	![images/g1.png](images/g1.png)       
	__(2.2) Set Group2 as selected:__       
	Use similar method to select Group2 cells.   
	![images/g2.png](images/g2.png)    
	Or, you can use select all the cells in an edge as group1/group2 using the "/Select Group1(2) as edge" button.    
	__(2.3) Get different genes between selected group1 and group2__     
	You can get the DE genes between two selected groups using "DE genes of selected cells" button.   
	![images/de.png](images/de.png)  
	Using these button, you can get all the up-regulated DE genes between group 1 and group 2.  
	If you want to get down-regulated DE genes between two groups, please switch the order of group1 and group2.     
	At the bottom of the DE gene list, you can see a button "functional analysis", which can be used to perform GO enrichment analysis for those genes.   
	![images/dego.png](images/dego.png)   
	__(2.4) Map to SRPING result__    
	You can map the model to SPRING visualization using the "SPRING Visualization" button.   
	![images/spring.png](images/spring.png)   
	__(2.5) Expression in selected groups__  
	You can check the expression of input genes in selected cells.   
	![images/ges.png](images/ges.png)   
	__(2.6) Clear all selected__  
	You can unset all the groups by using this button. 
	
# Gene Config Panel  
![images/geneconfig.png](images/geneconfig.png)   
You can explore the expression for specific genes using this Panel.   
(1) TFs are predicted based on their targets. Basically, target genes of a specific TF are differentially expressed in a specific path.

This inference is relying on the TF-gene interaction annotation we have.  If we don't have good annotation for some TFs, we might miss them in the prediction. 

That's why we also need eTF as a complement.   eTF stands for expressed-based TF.  In other words, if the expression of TF in a specific path/edge is significantly different compared with both parents and siblings.

In other words, it's unique.  We also list those TFs as eTFs.  

For the eTFs listed on N7, basically, it means the TFs regulating the path ending at N7 (P6). 

To be more accurate, the expression of ATF4, CEBPD and ZNF503 are "different" in P6 compared with the parent (P3) and the sibling (P9). z

# TF Config Panel  
![images/tfconfig.png](images/tfconfig.png)   
You can show/hide TFs and eTFs using this Panel.   
(1) TFs:  are predicted based on their targets.     
Basically, target genes of a specific TF are differentially expressed in a specific path.
This inference is relying on the TF-gene interaction annotation we have.  
If we don't have good annotation for some TFs, we might miss them in the prediction. 
That's why we also need eTF as a complement.   

(2) eTF stands for expressed-based TF.      
If the expression of a TF in a specific path/edge is significantly different 
compared with both parents and siblings, we regard it as an eTF.   
For example the eTFs listed on N7 means TFs regulating the path ending at N7 (P6), which was inferred based on their expression.  
To be more specific, the expression of ATF4, CEBPD and ZNF503 (at N7)s are "different" in P6 compared with the parent (P3) and the sibling (P9). 
  
# Download Panel  
![images/downloadconfig.png](images/downloadconfig.png)   
You can generate/download a figure for current visualization using the "Generate Figure" Button.  
With the "Download Selected Cells", you will be able to download the cells you selected using the "Cell Config Panel" (Group1 and/or Group2 cells).   

