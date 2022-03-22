# Cluph

In drug discovery and drug development, docking-based virtual screening and pharmacophore-based virtual screening assist the hit discovery and lead optimization in a fast and cost-efficient way. These methods represent a key tool in the field, which enables to computationally screen large libraries of compounds. However, the selection of appropriate candidates as well as dealing with false-positive still represents a challenge in the post-screening evaluation. 
In this context, **Cluph** is a chemoinformatic tool built for the intended purpose of finding hits compounds that come from virtual screening, both pharmacophore-based and docking-based. 

To use Cluph, molecular descriptor for each chemical entity should be generated (PaDEL can be used for this purpose). Then, **Principal Component Analysis (PCA)** is performed to reduce the dimension of the data by eliminating the least significant one. Finally, compounds are clusterized using **Hierarchical Clustering** method
