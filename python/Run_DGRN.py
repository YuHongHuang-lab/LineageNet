import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import os
import DGRN
import DGCN_data_prepare
import joblib
import sys


adata_path = sys.argv[1]
prior_network_path = sys.argv[2]

adata = sc.read_h5ad(adata_path)
all_genes = np.arange(adata.shape[1])

# 随机选择5000个基因索引
random_genes = np.random.choice(all_genes, size=7000, replace=False)

# 创建一个新的 adata 对象，只包含随机选择的基因
adata = adata[:, random_genes].copy()

prior_network =  pd.read_csv(prior_network_path)

prior_network[['from', 'to']] = prior_network[['from', 'to']].apply(lambda col: col.str.replace(':','-'))
data = DGCN_data_prepare.data_preparation(adata, prior_network)

# 设置昇腾NPU设备ID
NPU = 0
lineagenet_results_dict = {}
for li, data_li in data.items():
    # 创建模型
    # 参数：
    #   epoch：  单词训练的Epoch数
    #   repeats：重复训练数
    #   cuda：   训练设备ID
    #   sead：   随机种子
    lineagenet_GRN_model = DGRN.NetModel(epochs=350, repeats=3, cuda=NPU, seed=-1)
    # 开启训练并选择最优模型
    lineagenet_GRN_model.run(data_li)

    # 获得模型结果
    lineagenet_results = lineagenet_GRN_model.get_lineagenet_results(edge_threshold_avgDegree=8)
    lineagenet_results_dict[li] = lineagenet_results