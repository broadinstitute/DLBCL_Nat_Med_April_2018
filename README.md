This repository contains instructions on how to pull and run a docker image containing the scripts used to generate clusters used to identify DLBCL subtypes in the publication:



"Genetic and Functional Drivers of Diffuse Large B Cell Lymphoma." Chapuy, B., Stewart C., Dunford A. et al.

To run the module,

docker run  --mount type=bind,source=${PWD}/input_data,target=/Rscripts/input_data  --mount type=bind,source=${PWD}/output,target=/Rscripts/output  consensus_clustering Rscript run_consensus_clustering.R
