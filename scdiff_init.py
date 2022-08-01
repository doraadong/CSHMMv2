import scdiff.scdiff as S
import argparse
def run_scdiff_init(data_file,tfdna=None,output_file=None, config="auto",large=None,virtualAncestor=None,
                    label2cluster=False,nPCA=10):

    # 1) : read in gene expression
    AllCells=[]
    print("reading cells...")
    with open(data_file,'r') as f:
        line_ct=0
        for line in f:
            #print line
            if line_ct==0:
                S.GL=line.strip().split("\t")[3:]
                print("Number of genes: ", len(S.GL))
            else:
                line=line.strip().split("\t")
                iid=line[0]
                ti=float(line[1])
                li=line[2]
                ei=[round(float(item),2) for item in line[3:]]
                ci=S.Cell(iid,ti,ei,li,S.GL)
                if len(ei) != len(S.GL):
                    print("Different length, discarded: ",len(ei), ti, li)
                else:
                	AllCells.append(ci)
            line_ct+=1
            print('cell:'+str(line_ct))
            
    firstTime=min([float(item.T) for item in AllCells])
    
    print("config",config, "largeType", large, "label2cluster", label2cluster, "nPCA", nPCA, "ancestor", virtualAncestor)

    G1=S.Graph(AllCells,tfdna,config,largeType=large,dsync=None,virtualAncestor=virtualAncestor,
               label2cluster=label2cluster,nPCA=nPCA)

#     if large:
#         G1=S.Graph(AllCells,'auto','True')  #Cells: List of Cell instances 
#     else:
#         G1=S.Graph(AllCells,'auto',None)  #Cells: List of Cell instances 
    if output_file is None:
        out_file = open('init_cluster_'+data_file+'.txt','w')
    else:
        out_file = open(output_file, 'w')
    pairs=[]
    pairs2=[]
    print 'G1'
    for node in G1.Nodes:
        print 'ID: ', node.ID
        print 'index: ',G1.Nodes.index(node)
        ic=G1.Nodes.index(node)
        if node.P is not None:
            print 'P index: ',G1.Nodes.index(node.P)
            ip=G1.Nodes.index(node.P)
            pairs.append(str(ip)+' '+str(ic))
        else:
            print 'P index: ',None
        print 'ST: ',node.ST
        print 'T: ',node.T
        for cell in node.cells:
            print 'cell.ID: ', cell.ID
            pairs2.append(cell.ID+" "+str(ic))
    print pairs
    print pairs2
    out_file.write('\t'.join(pairs)+'\n')
    out_file.write('\t'.join(pairs2)+'\n')
    out_file.close()
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-tf',"--tf_dna_file", help="specify the tf_target file")
    parser.add_argument('-d',"--data_file", help="specify the data file")
    parser.add_argument('-o', "--output_file", required=False, default=None, help="specify the output file")
    parser.add_argument('-l',"--large_dataset", type=int, required=False, default=0, help="specify whether or not to speed up the initialization for large dataset (a different algorithm (PCA+k-means) will be applied)")
    parser.add_argument('-c', "--config_file", required=False, default=None, help="configuration file for known number of clusters at each time")
    parser.add_argument("-lc", '--labelAsCluster', required=False, default=0, type=int,
                        help="(1/0), optional, whether to use the input cell label as the cluster for initialization, default 0")
    parser.add_argument("-np", '--numberPCA', required=False, default=10, type=int,
                        help="integer, optional, number of principal components to use for large data for initialization, default 10")
    parser.add_argument("-a", '--virtualAncestor', required=False, default=None,
                        help="1 or None, optional, if using virtual ancestor or not")

    args=parser.parse_args()

    large=None
    if args.large_dataset>0:
        large="True"

    if args.labelAsCluster == 1:
        label2cluster = True
    else:
        label2cluster = False

    if args.data_file is not None:
        data_file=args.data_file
        output_file = args.output_file
        tfdna=args.tf_dna_file
        config_file = args.config_file
        nPCA = args.numberPCA
        ancestor = args.virtualAncestor

        if ancestor is not None:
            ancestor = str(ancestor)

        run_scdiff_init(data_file,tfdna,output_file=output_file, large=large, config=config_file,
                        label2cluster=label2cluster, nPCA=nPCA, virtualAncestor=ancestor)
        #if it takes too long, set large="True" for large dataset
    else:
        print "please use -d to specify the data file to run"
        print "please use -tf to specify the tf-target file to run"

