This repository contains instructions on how to build and run a docker image containing the scripts used to generate clusters used to identify DLBCL subtypes in the publication:



"Genetic and Functional Drivers of Diffuse Large B Cell Lymphoma." Chapuy, B., Stewart C., Dunford A. et al.

To run the module, first CD into the 'consensus_clustering' directory and build the docker container using the Dockerfile there:

cd consensus_clustering
docker build . -t consensus_clustering

Then run the provided wrapper script in the docker container using teh following command:


docker run --mount type=bind,source=${PWD}/input_data,target=/Rscripts/input_data --mount type=bind,source=${PWD}/output_dir,target=/Rscripts/output_dir dlbcl_consensus_clustering Rscript run_consensus_clustering.R
