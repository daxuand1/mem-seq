# MEM-SEQ
MEM-SEQ is a fast and scalable framework to select spatial eigenvectors from Moran's Eigenvector Maps (MEM) theory to adjust for spatial autocorrelation in spatial microbiome analysis.

This repository contains various MEM selection functions and simulation code for the MEM selection framework.

The key function is mem_seq(), which performs our proposed sequential selection approach. We also provide mem_fwd() and mem_glb() for the forward selection and global selection approaches.

We also provide simulate_spa_pat() function to simulate the spatial pattern described in our paper, and simulate_one() to run single simulation.

We provide a separate file toy_adh_seq.R to assess the SEQ method. To save time, it is based on a smaller sample size, different from the setting in our paper, which takes about 5 mins. The Type I error is preserved (about 0.03). You can replace the method argument to, for example, "fwd" to evaluate the FWD approach. It takes about 12 mins, and the Type I error is inflated (about 0.09).
