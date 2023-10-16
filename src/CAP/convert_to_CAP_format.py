# Format based on this
# https://github.com/evanbiederstedt/cap_file_format_liver/blob/main/cap_anndata_schema.md#dataset-specific-metadata

"""

Required obs (Cell Metadata)
Dataset-specific
    donor_id
    organism_ontology_term_id
    organism
    disease_ontology_term_id
    disease
    assay_ontology_term_id
    assay
    tissue_type
    tissue_ontology_term_id
    tissue
Cell-specific
    cell_annotation_set
    cell_fullname
    cell_ontology_term_id
    cell_ontology_term
    rationale
    rationale_dois
    marker_gene_evidence
    synonyms

Required obsm
    "X + _ + [EMBEDDING_TYPE] + _ + [SUFFIX]"

Required uns
    cell_annotation_schema_version
    a lot of CAP
    cap_dataset_timestamp (after submission)


Required var
    index?

"""

import os, glob
import scanpy as sc

abspath = '/share/quonlab/HLiCA/czi_liver_atlas/integration/harmony/integrated_h5'

DATASETS = os.listdir(abspath)

# for dataset in DATASETS:
#     dataset_dir = os.path.join(abspath,dataset)
#     h5ad_files = glob.glob(os.path.join(dataset_dir,"*_RNA.h5ad"))
#     for h5ad_file in h5ad_files:
#         adata = sc.read_h5ad(h5ad_file)
#         break

adata = sc.read_h5ad(os.path.join(abspath,'all_harmonized_seurat.h5ad'))

sc.pp.neighbors(adata)
sc.tl.leiden(adata)


# Load metadata
adata.obs['organism'] = 'Homo sapiens'
adata.obs['disease']  = 'normal'
adata.obs['assay']    = "10x 3' v2"
adata.obs['tissue']   = 'liver'

adata_gene = sc.read_10x_h5('/share/quonlab/HLiCA/czi_liver_atlas/data/LiverNetworkData/alignment/ref_GRCh38p13_gencode_v42/submissions-czi004liv/dasgupta_2023/XHL338/XHL338/outs/raw_feature_bc_matrix.h5')
name2id = dict(zip(adata_gene.var.index, adata_gene.var.gene_ids))

adata.var.index = adata.var.index.map(name2id)
adata.var.index = adata.var.index.astype('str')
adata.write_h5ad(os.path.join(abspath,'all_harmonized_seurat_CAP.h5ad'))

sample = 'DasGupta_XHH012'
adata_subset = adata[adata.obs['orig.ident'] == sample]
adata_subset.write_h5ad(os.path.join(abspath,f'all_harmonized_seurat_CAP_subset-{sample}.h5ad'))
