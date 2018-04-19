
docker build . -t dlbcl_consensus_clustering


docker run  --mount type=bind,source=${PWD}/input_data,target=/Rscripts/input_data  --mount type=bind,source=${PWD}/output,target=/Rscripts/output  dlbcl_consensus_clustering Rscript run_consensus_clustering.R
