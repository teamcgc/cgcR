---
title: "Kallisto-docker-container"
author: "Durga Addepalli"
date: "June 16, 2016"
output: html_document
---

```{r include=FALSE}
knitr::opts_chunk$set(eval = TRUE)
```

## INTRODUCTION
This guide will walk you through the steps to build a docker container for the tool "kallisto". kallisto is a program for quantifying abundances of transcripts from RNA-Seq data, or more generally of target sequences using high-throughput sequencing reads.You can read more about [Kallisto] here:(https://pachterlab.github.io/kallisto/)
(https://github.com/pachterlab/kallisto)

This tutorial assumes you have the docker installed and running.For documentation on docker installation refer to [Using Docker] (https://github.com/teamcgc/cgcR/blob/master/vignettes/UsingDocker.Part1.Rmd)

First lets copy the kallisto git repository from Github to any of your local folder!

```
git clone https://github.com/InSilicoDB/docker-kallisto.git kallisto-test

## Move into the dicretory created

cd kallisto-test

##  Build a new image from the source code in this directory witht he name 'kallisto' and version '0'
docker build -t kallisto:v0 .

## Once the build is done, you will need to get the Image id of the container and run the following command to be able to run kallisto
docker run -it kallisto:v0

## Push the image to docker hub
docker push durgaadd/kallisto:v0

```


## Mounting a volume (mounitng a local directory onto the docker container)

To mount a local directory onto the docker container, we need to edit the Docker file and add the following RUN command
```
vi Dockerfile 
RUN mkdir -p /tmp/files

## Specify the local directory with input data files

docker run -itv </Users/addepald/Documents/Data/>:/mnt/files kallisto:v0
```

Now inside the container (looks like root@<containerID>....), you can check to see the files/directory

```
cd /mnt/files

```

## Related Tutorials
[Install Docker Toolbox on Mac OS X](https://docs.docker.com/v1.10/mac/step_one/)

