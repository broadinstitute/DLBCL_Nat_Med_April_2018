This repository contains instructions on how to build and run a docker image containing the scripts used to generate clusters used to identify DLBCL subtypes in the publication:



Chapuy, Stewart, & Dunford et al. "Molecular subtypes of diffuse large B cell lymphoma are associated with distinct pathogenic mechanisms and outcomes." Nature Medicine, 30 April 2018. https://doi.org/10.1038/s41591-018-0016-8

Based on a modified consensus clustering algorithm first applied to expression data in the publication:

Brunet et al. "Metagenes and molecular pattern discovery using matrix factorization." PNAS March 23, 2004. 101 (12) 4164-4169; https://doi.org/10.1073/pnas.0308531101 

Please cite both publications if referring to this work.  


To run the module, first CD into the 'consensus_clustering' directory and build the docker container using the Dockerfile there:

cd consensus_clustering
docker build . -t dlbcl_consensus_clustering

Then run the provided wrapper script in the docker container using the following command:


docker run --mount type=bind,source=${PWD}/input_data,target=/Rscripts/input_data --mount type=bind,source=${PWD}/output_dir,target=/Rscripts/output_dir dlbcl_consensus_clustering Rscript run_consensus_clustering.R

