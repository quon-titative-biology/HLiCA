# (scpair) hongruhu@gpu-4-56:/group/gquongrp/workspaces/hongruhu/HLiCA/healthy_gamma$

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt



metadata = pd.read_csv('all_metadata.csv', index_col=0) # on_bad_lines='skip',
harmony_umap = pd.read_csv('all_harmony_umap.csv', index_col=0)
harmony_embedding = pd.read_csv('all_harmony_embedding.csv', index_col=0)
metadata.columns
metadata.STUDY.unique()
# array(['Mullen', 'DasGupta', 'Grun', 'Henderson', 'Scott', 'Toronto_1',
#        'Toronto_2', 'Toronto_3', 'Amit', 'Aronow', 'Cao', 'Chen', 'Fan',
#        'Fong', 'Heikenwalder', 'Li', 'Lleo', 'Qu', 'Schramm', 'Schwabe',
#        'Shi', 'Sun', 'Yan'], dtype=object)



harmony_umap = harmony_umap.loc[metadata.index]
harmony_embedding = harmony_embedding.loc[metadata.index]

obj = sc.AnnData(X=harmony_embedding.values, obs=metadata)
obj.obsm['X_umap'] = harmony_umap.values
obj.obsm['X_harmony'] = harmony_embedding.values

obj
# AnnData object with n_obs × n_vars = 844453 × 50
metadata.disease_status.value_counts()
# disease_status
# healthy    524699
# disease    319754

# drop ['old.ident'] from obs
obj.obs = obj.obs.drop(columns=['old.ident'])
obj.write("obj_harmony.h5ad")







#https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_06_celltyping.html
# 2.1 Label transfer


# def label_transfer(dist, labels, index):
#     lab = pd.get_dummies(labels)
#     class_prob = lab.to_numpy().T @ dist
#     norm = np.linalg.norm(class_prob, 2, axis=0)
#     class_prob = class_prob / norm
#     class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
#     # convert to df
#     cp_df = pd.DataFrame(
#         class_prob, columns=lab.columns
#     )
#     cp_df.index = index
#     # classify as max score
#     m = cp_df.idxmax(axis=1)
#     return m


# adata_ref = obj[obj.obs['disease_status'] == "healthy"]
# adata = obj[obj.obs['disease_status'] == "disease"]

# class_def = label_transfer(distances, adata_ref.obs['Gamma.Annotation'], adata.obs.index)
# obj.obs['predicted'] = class_def




# KNN
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

obj = sc.read("obj_harmony.h5ad")
# https://nbisweden.github.io/workshop-scRNAseq/labs/scanpy/scanpy_06_celltyping.html
# Ingest
adata_ref = obj[obj.obs['disease_status'] == "healthy"]
adata = obj[obj.obs['disease_status'] == "disease"]

sc.pp.neighbors(adata_ref, n_neighbors=15, use_rep='X_harmony')
sc.tl.umap(adata_ref)
sc.tl.ingest(adata, adata_ref, obs='Gamma.Annotation', embedding_method='umap')
sc.pl.umap(adata, color='Gamma.Annotation', wspace=0.5, save='_gamma_annotation_disease.png')
sc.pl.umap(adata_ref, color='Gamma.Annotation', wspace=0.5, save='_gamma_annotation_healthy.png')
adata.obs['pred_ingest'] = adata.obs['Gamma.Annotation'].copy()
adata.obs.to_csv('disease_metadata.csv')
sc.pl.umap(adata, color='pred_ingest', wspace=0.5, save='_pred_ingest_disease.png')

# SVM
import logging
import pickle
import numpy as np
from sklearn import svm
from sklearn.calibration import CalibratedClassifierCV

train_x = adata_ref.X
train_y = adata_ref.obs["Gamma.Annotation"].to_numpy()
classifier_dict = {
    "C": 1,
    "max_iter": 5000,
    "class_weight": "balanced",
}

clf = CalibratedClassifierCV(svm.LinearSVC(**classifier_dict))
clf.fit(train_x, train_y)

test_x = adata.X
adata.obs['svm_pred'] = clf.predict(test_x)

adata.obs['svm_pred' + "_probabilities"] = np.max(
                clf.predict_proba(test_x), axis=1
            )

sc.pl.umap(adata, color='svm_pred', wspace=0.5, save='_svm_pred_disease.png')
adata.obs[['pred_ingest', 'svm_pred']]
np.sum(adata.obs['pred_ingest'] == adata.obs['svm_pred'])
# 248754
np.sum(adata.obs['pred_ingest'] == adata.obs['svm_pred']) / len(adata.obs)
# 0.7779543023699469
adata.obs.to_csv('disease_metadata_ingest_svm.csv')
