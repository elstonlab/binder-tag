After running the cluster analysis pipeline, and the follow-up pooling and annotation, there are a few more post-processing steps for each quantity of interest.

Depending on the specific quantification, we:
- Filter out spurious clusters according to a few additional properties (too many simultaneous tracks; cluster present at start of imaging; detected incomplete segmentation).
- Incorporate additional data that wasn't directly tied to the cluster analysis output
- Apply additional statistical analysis (e.g. bootstrapping)

Explanation of folders:
- `codiff_in_clusters` estimates the percentage of co-diffusing tracks that are cluster associated, and computes a simulated randomized control (Fig. 6E)
- `lifetimes` estimates track lifetime for clustered and non-clustered tracks (Fig. S6R)
- `sizes` estimates the sizes of clusters (Fig. 6G, Fig. S6J)
- `sizes_slow` estimates the sizes of clusters after segregating by diffusional state (Fig. 6I, Fig. S6L-M) 
