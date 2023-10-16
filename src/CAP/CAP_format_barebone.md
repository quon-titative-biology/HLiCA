
## Cell Metadata (`obs`)

- The file must contain the following fields in `obs` with the column names exactly as listed below.
- These fields may not contain NA values.

| Field name  | Type  | Accepted values  |
|:----------|:----------|:----------|
| assay    | string    | <ul><li>This must be the specific assay, not a generic term such as scRNA-seq or 10x sequencing </li><li> Accepted value examples: `10x 3' v2`, `10x 3' v3`, `Smart-seq2`    |
| disease   | string    | <ul><li>This must be the specific disease term, e.g. glioblastoma or Alzheimer disease</li><li>For healthy samples use `normal` or `healthy`    |
| organism    | string    | <ul><li>Must be the Latin (Genus species) name, e.g. `Homo sapiens`, `Mus musculus`    |
| tissue    | string    | <ul><li>The most accurate anatomical term for where the sample was collected from </li><li> Accepted value examples: `retina`, `heart left ventricle`  |


### Optional `obs` field: clustering
- Cluster fields are not required but if they are in the dataset they must be denoted with the prefix `cluster_` and the field name may also include the `algorithm_type` and `suffix`.
- Examples of acceptable cliuster column names are: `cluster_leiden`,  `cluster_leiden_broad`,  `cluster_louvain_precise3`, `cluster_fine`, etc.

## Embedding (`obsm`)

- At least one embedding, tSNE, UMAP or PCA, is required, and more than one may be included.
- The embeddings(s) must be saved with the prefix `X_`, for example: `X_tsne`, `X_pca`, `X_umap`.

## Gene Metadata (`var`)
- CAP requires that gene names by ENSEMBL terms. These MUST be encoded in the [index](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.index.html) of the `var` fields following the [AnnData standard](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.html).
- `var` is a `pandas.DataFrame` object. The gene names MUST be used to index these rows, i.e. `pandas.DataFrame.index`.
- *Note:* the UI will convert the ENSEMBL terms to common gene names based on the organism specified. We currently support `Homo sapiens` and `Mus musculus`. If there are other species you wish to upload to CAP, please contact `support@celltpye.info` and we will work to accommodate your request.

## Count Matrix (`X`)

- The file must contain a raw count matrix saved as `.X` or `raw.X`.
- If the file contains a count matrix in `.X` the `.raw` layer must be empty.

## Dataset-wide Metadata (`uns`)
- There are no requirements for fields in uns for uploading to CAP.
