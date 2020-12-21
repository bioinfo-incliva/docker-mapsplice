# docker-mapsplice
This repository contains a dockerfile to create a docker that will compile MapSplice_multi_threads_2.0.1.9 without errors with Ubuntu version 12.04 

Mapsplice is a software for mapping RNA-seq data to reference genome for splice junction discovery that depends only on reference genome, and not on any futher annotations. 

## Build docker
  
    sudo docker build -t  mapsplice .
 
## Run docker

    sudo docker run -ti -v /media/scratch/:/media/scratch mapsplice

## Copy the docker directory to the server
Since we want to run Mapsplice on the server, what we do is copy the compiled mapsplice folder to the server:

Get id_docker:
  
    sudo docker ps
    
Copy the docker directory to the server:

    sudo docker cp id_docker:/usr/local/src/MapSplice_multi_threads_2.0.1.9/ /usr/local/src/MapSplice_multi_threads_2.0.1.9
