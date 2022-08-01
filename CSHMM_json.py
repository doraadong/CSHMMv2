"""
Updated 2-13-20: add mean expression to edges.
Updated 8-12-20: get RTFs (combined with getRTF.py); save a new model file with RTFs and output into data.json.

"""
import sys
import json
import pickle

import numpy as np
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
from scipy.stats import ttest_ind

# import CSHMM_analysis_modified_20181018 as ANA
import CSHMM_train_modified_20180403 as ML


def tellDifference(exps, exps_compare, j, fcut=0.6):
    """

    Modified based on SCDIFF's tellDifference function.
    Author: Jun Ding
    Project: https://github.com/phoenixding/scdiff

    """

    if len(exps_compare) == 0:
        return [0, 1, 0]

    # get expression for jth gene
    exp = exps[:, j]
    exp_compare = exps_compare[:, j]

    fc = np.mean(exp) - np.mean(exp_compare)
    pv = ttest_ind(exp, exp_compare)[1]

    pcut = 0.05
    if pv < pcut and fc > fcut:
        return [1, pv, fc]
    if pv < pcut and fc < -1*fcut:
        return [-1, pv, fc]
    return [0, pv, fc]


def getMeanPath(model, hid_var, cell_exps):
    cluster_names = list(np.unique(hid_var['cell_path']))
    print ('All paths: ', cluster_names)
    cluster_cellExps_dict = {}

    # get reference to the cluster of cells assigned to each cluster/path
    for cluster in cluster_names:
        cluster_cellExps_dict[cluster] = list(cell_exps[hid_var['cell_path'] == cluster])

    # initialize dictionary to store average gene expression per path
    mean_paths = {}

    # for each path
    for i, path in enumerate(model['path_info']):
        # get expression of cells assigned to the current path
        exps_cur = np.array(cluster_cellExps_dict[path['ID']])
        mean_paths[path['ID']] = exps_cur.mean(axis=0)

    return mean_paths


def get_sibling_path_idx(p_idx, path_info):
    """
    Get path p_idx's siblings.
    Modified from Chieh Lin's CSHMM_train.py.

    """
    sib_idxs = []

    for i, p in enumerate(path_info):
        if p_idx in p["child_path"]:
            sib_idxs = list(p["child_path"])
            sib_idxs.remove(p_idx)  # remove the path itself

    return np.array(sib_idxs)

def getRTF(model, hid_var, cell_exps, gene_names, model_file, fcut):
    """

    Compare with its siblings (branch from the same parent).

    """

    cluster_names = list(np.unique(hid_var['cell_path']))
    print ('All paths: ', cluster_names)
    cluster_cellExps_dict = {}

    # get reference to the cluster of cells assigned to each cluster/path
    for cluster in cluster_names:
        cluster_cellExps_dict[cluster] = list(cell_exps[hid_var['cell_path'] == cluster])

    # initialize eTF list to store all eTFs for different paths
    eTFs = []
    eTF_count = 0

    # for each path
    for i, path in enumerate(model['path_info']):
        # get expression of cells assigned to the current path
        exps_cur = np.array(cluster_cellExps_dict[path['ID']])

        # get expression of its parent' path
        try:
            exps_parent = np.array(cluster_cellExps_dict[path['parent_path']])
        except KeyError:
            exps_parent = []

        path_siblings = get_sibling_path_idx(path['ID'], model['path_info'])

        # get expression from siblings
        exps_siblings = [cluster_cellExps_dict[k] for k in path_siblings]
        if len(exps_siblings) > 0:
            exps_siblings = np.vstack(exps_siblings)

        peTFs = []
        # for each TF
        for j, tf in enumerate(gene_names):
            [jflag1, pvp, fcp] = tellDifference(exps_cur, exps_parent, j, fcut)
            [jflag2, pvs, fcs] = tellDifference(exps_cur, exps_siblings, j, fcut)

            # case 1: significantly lower than both parent/siblings;
            # case 2: significantly higher than both parent/siblings;
            # case 3: significantly higher or lower than parent;
            if (jflag1 * jflag2 > 0) or (jflag1 != 0 and len(exps_siblings) == 0):
                peTFs.append([pvp, j, fcp, fcs])
        # sort by p-value
        peTFs.sort()

        # update the model
        path['eTF_list'] = [gene_names[item[1]] for item in peTFs]
        path['eTF_method'] = np.repeat("pv_cutoff_0.05_FC_cutoff_" + str(fcut), len(peTFs))
        path['eTF_pv'] = [item[0] for item in peTFs]
        path['eTF_fcp'] = [item[2] for item in peTFs]
        path['eTF_fcs'] = [item[3] for item in peTFs]

        # add into eTFs list
        for j, item in enumerate(peTFs):
            eTFs += [{"etf_name": gene_names[item[1]], "fcp": item[2], "fcs": item[3], "id": eTF_count,
                      "location": path["ID"], "p-value": item[0]}]
            eTF_count += 1

    # save updated model with eTF
    filename = model_file.split('.')[0] + '_etf' + '.pickle'
    with open(filename, 'wb') as handle:
        pickle.dump(model, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return eTFs



def json_visualization( gene_names,cell_names,cell_exps,cell_day,model=None, hid_var=None, tf_file=None, model_file=None,
                        eTF_index=None, fcut=0.5, use_time_as_label=False):

    # generate visualization
    xmatrix = np.array(cell_exps)

    # pca
    pca = PCA(n_components=50)
    pca2 = PCA(n_components=2)
    xpca2_matrix = pca2.fit_transform(xmatrix)
    xpca_matrix = pca.fit_transform(xmatrix)
    # pdb.set_trace()

    # tsne
    tnse = TSNE(n_components=2)
    xtsne_matrix = tnse.fit_transform(xpca_matrix)

    # isomap
    isomap = Isomap(n_components=2)
    isomap_matrix = isomap.fit_transform(xmatrix)


    if model is None or hid_var is None:
        model,hid_var=ML.load_model(model_file)

    path_info=model["path_info"]
    g_param=model["g_param"]
    sigma_param=model['sigma_param']
    K_param=model['K_param']
    gene_start_time=model['gene_start_time']
    cell_path = hid_var["cell_path"]
    cell_time = hid_var["cell_time"]
    cell_labels = hid_var["cell_labels"]

    print(model_file)

    if use_time_as_label:
        cell_labels = cell_day

    edges=[]
    node_dict={}
    tfs=[]

    tf_count=0
    ML.dTD=ML.getdTD(tf_file)
    TF_names = np.array(map(lambda x: x.lower(),ML.dTD.keys()))

    # for each path, get mean expression value of a gene
    mean_paths = getMeanPath(model, hid_var, cell_exps)

    # for each path
    for i, path in enumerate(path_info):
        print(mean_paths[path["ID"]].shape)
        edges+=[{"id":path["ID"], "source_node":path["Sp_idx"], "end_node":path["Sc_idx"], \
                 "expression":mean_paths[path["ID"]].tolist()}]
        node_dict[path["Sp_idx"]]={"id": path["Sp_idx"], "time": path["level"], "expression":g_param[path["Sp_idx"],:].tolist()}
        node_dict[path["Sc_idx"]]={"id": path["Sc_idx"], "time": path["level"]+1, "expression":g_param[path["Sc_idx"],:].tolist()}
        if "TF_list" in path.keys():
            for tf,tfp in zip(path["TF_list"],path["TF_pvalue"]):
                tf_idx=TF_names.tolist().index(tf)
                tfs+=[{ "id":tf_count, "tf_name": tf, "location": (path["ID"], model["TF_start_time"][path["ID"],tf_idx]), "p-value":tfp}]
                tf_count+=1
    nodes=node_dict.values()

    cells = []
    vizCells = []

    for i,(cn,ce,cl,cp,ct) in enumerate(zip(cell_names,cell_exps,cell_labels,hid_var["cell_path"],hid_var["cell_time"])):
        p=cell_path[i]
        t=cell_time[i]
        Sp_idx=path_info[p]['Sp_idx']
        Sc_idx=path_info[p]['Sc_idx']
        g_a=g_param[Sp_idx]
        g_b=g_param[Sc_idx]
        mu_x_i=g_b+(g_a-g_b)*np.exp(-K_param[p]*np.maximum(t-gene_start_time[p],0))  # mean learned from the model
        diff = cell_exps[i]-mu_x_i
        ssq = np.sum(diff**2)
        sign = np.sign(np.sum(diff))

        cells += [{"id": cn, "type": cl, "location": (cp, ct, sign * ssq / len(gene_names))}]

        # visualization
        ixtsne = xtsne_matrix[i]
        ixpca = xpca2_matrix[i]
        ixisomap = isomap_matrix[i]

        jci = { "ID":cn, "typeLabel":cl, "T": ct}
        jci["TE"] = list(ixtsne)
        jci["PE"] = list(ixpca)
        jci["IE"] = list(ixisomap)
        vizCells.append(jci)

    genes=gene_names.tolist()

    # get eTFs
    # filter to keep TF expressions only
    cell_exps_exps = cell_exps[:, eTF_index]
    gene_names_exps = np.array(gene_names_u)[eTF_index]

    eTFs = getRTF(model, hid_var, cell_exps_exps, gene_names_exps, model_file, fcut)

    # output data.json
    output=open(model_file + "_data.json", "w")
    output.write("data=")
    output.write(json.dumps([nodes, edges, cells, genes, tfs, eTFs], sort_keys=True,indent=4))
    output.write('\n'+"data=JSON.stringify(data);")
    output.close()

    # output visualization json file
    output = open(model_file + "_CellViz.json", "w")
    output.write("drdata=")
    output.write(str(vizCells))
    output.write('\n' + "drdata=JSON.stringify(drdata);")
    output.close()

if __name__=="__main__":
    if len(sys.argv) < 4:
        print ("usage: python CSHMM_json.py data_file_name model_file_name tf_file "
               "number_of_genes_training fold_change_cut")
        sys.exit()

    data_file = sys.argv[1]
    model_file = sys.argv[2]
    tf_file = str(sys.argv[3])  # tf-DNA interaction file
    ng = int(sys.argv[4])
    fcut = float(sys.argv[5])
    use_time_as_label = str(sys.argv[6])
    print(data_file, model_file, tf_file, ng, fcut, use_time_as_label)

    if use_time_as_label == "True":
        use_time_as_label = True
    else:
        use_time_as_label = False

    # ng = 100000
    # if len(sys.argv)>3 and sys.argv[3]=='-ng':
    #     ng = int(sys.argv[4])

    cell_names,cell_day,cell_labels,cell_exps,gene_names=ML.load_data(data_file,ng)

    gene_names_u = [item.upper() for item in gene_names]

    # get TF full list
    tflistpath = 'input/tool/Human_mouse_TFList.txt'

    try:
        with open(tflistpath, 'r') as f:
            TFs = f.readlines()
            TFs = [item.strip().split()[0] for item in TFs]
    except:
        print("error! Please check your input TF List file")
        sys.exit(0)

    # get index in the expression for the TFs
    eTF_index = [gene_names_u.index(item) for item in gene_names_u if item in TFs]

    json_visualization(gene_names, cell_names, cell_exps, cell_day, model_file=model_file, tf_file=tf_file,
                       eTF_index=eTF_index, fcut=fcut, use_time_as_label=use_time_as_label)


