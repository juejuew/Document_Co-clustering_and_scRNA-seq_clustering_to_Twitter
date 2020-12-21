# Document_Co-clustering_and_scRNA-seq_Clustering_to_Twitter

Traditional topic modeling schemes including Latent Dirichlet allocation (LDA) and LatentSemantic Analysis (LSA) have been widely adopted in the social sciences, but are typically unable to cope well with short and noisy texts. Several document co-clustering algorithms havebeen proposed as potential resolutions, but their performances have seldom been comparedsystematically. Here, we assess the effectiveness ofthe existing document co-clustering algorithms using LDA, one of the most commonly usedtopic modeling techniques, as a baseline for comparison. 

In addition, we observe that microblog data such as Twitter data, once formatted as aword-document matrix, shares important properties with single-cell RNA sequencing data; in particular, both types of data are high dimensional, consist of small integer counts, and havesparse non-zero entries. Moreover, there are several advanced clustering methods designed toidentify cell types and gene expression patterns, so we apply single-cell RNA-seq data clustering methods to Twitter data. 

## Related Links:
### Topic Modeling algorithms:
- LDA: https://github.com/deekshakoul/LDA-Topic-modelling
- Spectral Co-clustering: https://scikit-learn.org/stable/auto_examples/bicluster/plot_bicluster_newsgroups.html#sphx-glr-auto-examples-bicluster-plot-bicluster-newsgroups-py
- Information Theoretical Co-clustering (ITCC): https://hal.archives-ouvertes.fr/hal-01804528/document
- Bregman Block Average Co-clustering (BBAC): https://github.com/felipeyanez/bbac
- Model-Based Co-clustering(LBM): https://cran.r-project.org/web/packages/blockcluster/blockcluster.pdf
- Directional Co-clustering with a Conscience (DCC): https://github.com/saghiles/dcc
- Ferg's Method: https://deepblue.lib.umich.edu/handle/2027.42/163193
### scRNA-Seq data clustering:
- Monocle3: https://cole-trapnell-lab.github.io/monocle3/docs/clustering/
- Seurat: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html

# References
Blei, David M., Andrew Y. Ng, and Michael I. Jordan. "Latent dirichlet allocation." Journal of machine Learning research 3.Jan (2003): 993-1022.

Dhillon, Inderjit S. "Co-clustering documents and words using bipartite spectral graph partitioning." Proceedings of the seventh ACM SIGKDD international conference on Knowledge discovery and data mining. 2001.

Dhillon, Inderjit S., Subramanyam Mallela, and Dharmendra S. Modha. "Information-theoretic co-clustering." Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining. 2003.

Iovleff, S., et al. "blockcluster: Coclustering Package for Binary, Categorical, Contingency and Continuous Data-Sets." R package version 4.2 (2015).

Role, François, Stanislas Morbieu, and Mohamed Nadif. "Coclust: a python package for co-clustering." (2018).

Banerjee, Arindam, et al. "A generalized maximum entropy approach to bregman co-clustering and matrix approximation." Journal of Machine Learning Research 8.Aug (2007): 1919-1986.

Salah, Aghiles, and Mohamed Nadif. "Model-based von Mises-Fisher Co-clustering with a Conscience." Proceedings of the 2017 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2017.

Levine, Jacob H., et al. "Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis." Cell 162.1 (2015): 184-197.

Ferg, Robyn. Modern Survey Estimation with Social Media and Auxiliary Data. Diss. 2020.
