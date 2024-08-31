import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse
import networkx as nx
import scanpy as sc
from typing import Optional, Union
from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf
import matplotlib.ticker as ticker
import seaborn as sns

from .resources import TFs_human, TFs_mouse


def data_preparation(input_expData: Union[str, sc.AnnData, pd.DataFrame],
                     input_priorNet: Union[str, pd.DataFrame],
                     genes_DE: Optional[Union[str, pd.DataFrame, pd.Series]] = None,
                     additional_edges_pct: float = 0.01) -> dict[str: AnnData]:
    """
    Prepare the data object for DGRN.
    """

    print('[0] - Data loading and preprocessing...')

    ## [1] Single-cell RNA-seq data
    lineages = None
    if isinstance(input_expData, str):
        p = Path(input_expData)
        if p.suffix == '.csv':
            adata = sc.read_csv(input_expData, first_column_names=True)
        else:  # h5ad
            adata = sc.read_h5ad(input_expData)
    elif isinstance(input_expData, sc.AnnData):
        adata = input_expData
        lineages = adata.uns.get('lineages')
    elif isinstance(input_expData, pd.DataFrame):
        adata = sc.AnnData(X=input_expData)
    else:
        raise Exception("Invalid input! The input format must be '.csv' file or '.h5ad' "
                        "formatted file, or an 'AnnData' object!", input_expData)

    possible_species = 'mouse' if bool(re.search('[a-z]', adata.var_names[0])) else 'human'

    # Gene symbols are uniformly handled in uppercase
    adata.var_names = adata.var_names.str.upper()

    ## [2] Prior network data
    if isinstance(input_priorNet, str):
        netData = pd.read_csv(input_priorNet, index_col=None, header=0)
    elif isinstance(input_priorNet, pd.DataFrame):
        netData = input_priorNet.copy()
    else:
        raise Exception("Invalid input!", input_priorNet)
    # make sure the genes of prior network are in the input scRNA-seq data
    netData['from'] = netData['from'].str.upper()
    netData['to'] = netData['to'].str.upper()
    netData = netData.loc[netData['from'].isin(adata.var_names.values)
                          & netData['to'].isin(adata.var_names.values), :]
    netData = netData.drop_duplicates(subset=['from', 'to'], keep='first', inplace=False)

    # Transfer into networkx object
    priori_network = nx.from_pandas_edgelist(netData, source='from', target='to', create_using=nx.DiGraph)
    priori_network_nodes = np.array(priori_network.nodes())

    # in_degree, out_degree (centrality)
    in_degree = pd.DataFrame.from_dict(nx.in_degree_centrality(priori_network),
                                       orient='index', columns=['in_degree'])
    out_degree = pd.DataFrame.from_dict(nx.out_degree_centrality(priori_network),
                                        orient='index', columns=['out_degree'])
    centrality = pd.concat([in_degree, out_degree], axis=1)
    centrality = centrality.loc[priori_network_nodes, :]

    ## [3] A mapper for node index and gene name
    idx_GeneName_map = pd.DataFrame({'idx': range(len(priori_network_nodes)),
                                     'geneName': priori_network_nodes},
                                    index=priori_network_nodes)

    edgelist = pd.DataFrame({'from': idx_GeneName_map.loc[netData['from'].tolist(), 'idx'].tolist(),
                             'to': idx_GeneName_map.loc[netData['to'].tolist(), 'idx'].tolist()})

    ## [4] add TF information
    is_TF = np.ones(len(priori_network_nodes), dtype=int)
    if possible_species == 'human':
        TFs_df = TFs_human
    else:
        TFs_df = TFs_mouse
    TF_list = TFs_df.iloc[:, 0].str.upper()
    is_TF[~np.isin(priori_network_nodes, TF_list)] = 0

    # Only keep the genes that exist in both single cell data and the prior gene interaction network
    adata = adata[:, priori_network_nodes]
    if lineages is None:
        cells_in_lineage_dict = {'all': adata.obs_names}  # all cells are regarded as in on lineage if lineages is None
    else:
        cells_in_lineage_dict = {l: adata.obs_names[adata.obs[l].notna()] for l in lineages}

    print(f"Consider the input data with {len(cells_in_lineage_dict)} lineages:")
    adata_lineages = dict()
    for l, c in cells_in_lineage_dict.items():
        print(f"  Lineage - {l}:")

        adata_l = sc.AnnData(X=adata[c, :].to_df())
        adata_l.var['is_TF'] = is_TF
        adata_l.varm['centrality_prior_net'] = centrality
        adata_l.varm['idx_GeneName_map'] = idx_GeneName_map
        adata_l.uns['name'] = l

        ## [5] Additional edges with high spearman correlation
        if isinstance(adata_l.X, sparse.csr_matrix):
            gene_exp = pd.DataFrame(adata_l.X.A.T, index=priori_network_nodes)
        else:
            gene_exp = pd.DataFrame(adata_l.X.T, index=priori_network_nodes)

        ori_edgeNum = len(edgelist)
        edges_corr = np.absolute(np.array(gene_exp.T.corr('spearman')))
        np.fill_diagonal(edges_corr, 0.0)
        x, y = np.where(edges_corr > 0.6)
        addi_top_edges = pd.DataFrame({'from': x, 'to': y, 'weight': edges_corr[x, y]})
        addi_top_k = int(gene_exp.shape[0] * (gene_exp.shape[0] - 1) * additional_edges_pct)
        if len(addi_top_edges) > addi_top_k:
            addi_top_edges = addi_top_edges.sort_values(by=['weight'], ascending=False)
            addi_top_edges = addi_top_edges.iloc[0:addi_top_k, 0:2]
        edgelist = pd.concat([edgelist, addi_top_edges.iloc[:, 0:2]], ignore_index=True)
        edgelist = edgelist.drop_duplicates(subset=['from', 'to'], keep='first', inplace=False)
        print('    {} extra edges (Spearman correlation > 0.6) are added into the prior gene interaction network.\n'
              '    Total number of edges: {}.'.format((len(edgelist) - ori_edgeNum), len(edgelist)))

        adata_l.uns['edgelist'] = edgelist

        ## [6] Differential expression scores
        logFC = adata.var.get(l + '_logFC')
        if (genes_DE is None) and (logFC is None):
            pass
        else:
            if genes_DE is not None:
                if isinstance(genes_DE, str):
                    genes_DE = pd.read_csv(genes_DE, index_col=0, header=0)
            else:  # logFC is not None
                genes_DE = logFC

            genes_DE = pd.DataFrame(genes_DE).iloc[:, 0]
            genes_DE.index = genes_DE.index.str.upper()
            genes_DE = genes_DE[genes_DE.index.isin(priori_network_nodes)].abs().dropna()
            node_score_auxiliary = pd.Series(np.zeros(len(priori_network_nodes)), index=priori_network_nodes)
            node_score_auxiliary[genes_DE.index] = genes_DE.values
            node_score_auxiliary = np.array(node_score_auxiliary)
            adata_l.var['node_score_auxiliary'] = node_score_auxiliary
            genes_DE = None

        adata_lineages[l] = adata_l
        print(f"    n_genes × n_cells = {adata_l.n_vars} × {adata_l.n_obs}")

    return adata_lineages