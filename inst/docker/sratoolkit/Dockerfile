####################################################################################
#
# SRA Toolkit version 2.6.3 installed on ubuntu
#
# The SRA Toolkit tools are on the path, so you should be able to do:
#
#    docker run seandavi/sratoolkit fastq-dump SRR789392
#
# Resulting fastq file will end up in the / directory as SRR789392.fastq
#
####################################################################################
from ubuntu

MAINTAINER Sean Davis <seandavi@gmail.com>

RUN apt-get update  && apt-get install -y \
    wget

RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.6.3/sratoolkit.2.6.3-ubuntu64.tar.gz

RUN tar -xvzf sratoolkit.2.6.3-ubuntu64.tar.gz

ENV PATH $PATH:/sratoolkit.2.6.3-ubuntu64/bin
