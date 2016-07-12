---
title: "FireCloudWDL_RunTool_Kallisto"
author: "Durga Addepalli"
date: "July 12, 2016"
output: html_document
---
## Introduction
This tutorial provides the steps for running the tool 'kallisto' by creating a WiDdLe (WDL) file to run analysis on the FireCloud and follows the tutorial 'A basic introduction to the FireCloud' which can be found [here] {https://github.com/teamcgc/cgcR/blob/master/vignettes/Firecloud.intro.Rmd}

This tutorial assumes you have all the dependencies listed in the basic introduction to run the WDL tool and Cromwell.

To rerun the basic requirements run the following commands:

```{r}
brew install Caskroom
brew install cask

## Cromwell requires java 1.8+, which can be installed by running the command
brew cask install java

brew install cromwell
brew install sbt

git clone https://github.com/broadinstitute/cromwell.git cromewell
cd cromewell/
sbt assembly
cd ..
git clone https://github.com/broadinstitute/wdltool wdltool
cd wdltool/
sbt assembly

```

Docker needs to be running, so as to pull the docker image of the tool kallisto, and also mount the input files. 

If you have the Docker toolbox installed you could just run this tutotrial on the Docker quickstart terminal. If you have the docker machine installed run the following commands:

```{r}
docker-machine start dev
docker-machine regenerate-certs dev
docker-machine env dev
eval $(docker-machine env dev)
```

## Creating the WDL Script for running kallisto:

kallisto quantifies read files directly without the need for read alignment, but it does perform a procedure called pseudoalignment. Pseudoalignment requires processing a transcriptome file to create a “transcriptome index”. 

Below is a simple WDL script (kallisto.test.wdl) for running kallisto commands. The example here shows 
- Creating a task called kallisto
- the input file variable
- the command you want to run (kallisto index) with arguments and inputs
- the runtime block which has the information needed to pull the image from docker hub
- the output format
- Creating a workflow (called kalWF) which calls the task 'kallisto'

```{r}
task kallisto {

        File inkal

        command {
		             kallisto index -i kalout.idx ${inkal}
            		}
        runtime
                {
                docker : "durgaadd/kallisto:0.43.0"
                }
        output {
                File response = stdout()
                }

        }

workflow kalWF
        {
        File inkal

        call kallisto {
                        input:
                        inkal=inkal
                      }
        }
```

The input file for this kallisto command is a fasta file (inkal.fasta), which is like this:

```{r}
>GENSCAN00000000001 cdna:genscan chromosome:GRCh38:5:122151991:122153085:1 transcript_biotype:protein_coding
ATGGAAAGAGGAAAGAAGAAAAGAATTTCCAATAAGTTACAACAAACTTTTCACCATTCTAAAGAACCCACTTTCCTTATCAACCAAGCTGGGCTTCTCTCTAGTGACTCCTATTCTAGC
>GENSCAN00000000002 cdna:genscan chromosome:GRCh38:5:122675795:122676286:1 transcript_biotype:protein_coding
ATGATGAACAGAATGGCCCCAGAGAATTTCCAGCCAGACCCTTTCATCAACAGGAATGATTCCAACATGAAGTATGAAGAGCTAGAAGCTCTGTTTAGCCAGACTATGTTCCCAGATAGA
>GENSCAN00000000003 cdna:genscan chromosome:GRCh38:5:121876146:121916483:-1 transcript_biotype:protein_coding
ATGGATGACTCTAAGGGCAATGGAAAGAGGGCTAAGATTAGAGGTAAGGGTCCCAAGATATTCCTCAAGAGTCTCCTGGCCACACTGCCAAACACATCATATGTCTGTGCCTCAGAACCT
>GENSCAN00000000004 cdna:genscan chromosome:GRCh38:5:122174100:122189697:-1 transcript_biotype:protein_coding
ATGAAGGAATATCTGGATCATGGAGCACTCGAGTTTTTGCTCCAACAGAAACAGTGGAGCTGTTTTGACTCCACTGCGCAGTGGTGGGCAGAAGGTGGCAATGGAGACTGCAGAAGAAAC
>GENSCAN00000000005 cdna:genscan chromosome:GRCh38:5:122391135:122489431:1 transcript_biotype:protein_coding
ATGGAAGCCCCTGAATACCTTGATTTGGATGAAATTGACTTTAGTGATGACATATCTGACAATAGGAGTCAAGGGAACAGGCTACAAAAGCTTGGATTGGAGGACACAGACAGGGAAGATGCAATGGGCTTTGGTTCCCATAGGGCCAAACTGACAGTAGTTGCTGCCCTGGGAGCTTGC
>GENSCAN00000000006 cdna:genscan chromosome:GRCh38:5:121851964:121852692:1 transcript_biotype:protein_coding
ATGCTGTCCTGCTTCAGGCTCCTCTCCAGGCACATCAGCCCTTCGCTGGCGTCTCTGCGCCCGGTGCGCTGCTGCTTCGCGCTCCCGCTGCGTTGGGCCCCGGGGCGCCCCTTGGACCCC
>GENSCAN00000000007 cdna:genscan chromosome:GRCh38:5:122307840:122334010:1 transcript_biotype:protein_coding
ATGAAGAACAGATCAGGTGAGGAATGTCAAATTGTTAGTCACATCCAGCTCGAGCACACAGGGACACACATGGCAAGAAAGATGTGCCCTGGCCTAAAGAAAAAGATAAGGGAGTTTATGCTGGAACCCTGGAATGGTCTGGGGACCAATGAGATGGCAGCGGTTAACGCTTGGATCACG
>GENSCAN00000000008 cdna:genscan chromosome:GRCh38:5:121975382:122027126:1 transcript_biotype:protein_coding
ATGAAGGAATTGAAACCTGACATAGTAACTAAATCTGCTCTTGGTGATGATATCAACTTTGAAAAAATCTTCAAAAAGCCAGATTCTACTGCAACTGAAAGAGCAATTGCCAGACTAGCAGTACATCCTCTTCTGAAGAAAAAGATAGATGTGCTAAAAGCTGCTGTACAAGCCTTTAAAGAAGCAAGACAAAATGTTGCTGAAGTTGAGTCATCAAAGAATGCTTCAGAGGACAATCAT
>GENSCAN00000000009 cdna:genscan chromosome:GRCh38:5:121746572:121763132:-1 transcript_biotype:protein_coding
ATGGCCAGGGGAAGGATCGGCAAATATTTTGGGAAAGATTCGAGACAGGAGCCAGTGGCACCCATAGAAACAGAATGTAGAGTTCAGAAAGAACAATAA
>GENSCAN00000000010 cdna:genscan chromosome:GRCh38:5:122498332:122638569:-1 transcript_biotype:protein_coding
ATGGTCTCCAGAGAGCCTGAGCGTCTGCTGGGGAACTTCCAGAAGCTTCTGTCTGAGATTGGTTCCATGCTTGTAACTGGGACCAGCAAATGGCAAAATAAAGGACAGATACGGAGTACT
>GENSCAN00000000011 cdna:genscan chromosome:GRCh38:5:122065891:122077985:-1 transcript_biotype:protein_coding
ATGCGCTTCGCCTGGACCGTGCTCCTGCTCGGGCCTTTGCAGCTCTGCGCGCTAGTGCACTGCGCCCCTCCCGCCGCCGGCCAACAGCAGCCCCCGCGCGAGCCGCCGGCGGCTCCGGGCGCCTGGCGCCAGCAGATCCAATGGGAGAACAACGGGCAGGTGTTCAGCTTGCTGAGCCTGGGCTCACAGTACCAGCCTCAGCGCCGCCGGGACCCGGGCGCCGCCGTCCCTGGTGCAGCC
>GENSCAN00000000012 cdna:genscan scaffold:GRCh38:KI270751.1:151:32727:1 transcript_biotype:protein_coding
ATGGCACAAGTTGCAGTTTCCACCCTGCCCATTGAAGATGAGGAGTCTGTTGAAGATGAGGAGTCCTTGGAGAGCAGGATGGTGGTGACATTCCTGTCAGCTCTCGACTCCATGAAAGTGTCAGAGCCTTTACGTGGACCTTCTCATGAAAACGGAAACAGAATAGTCAATGGAAAAGGAGAAGAAACAAATGCTGTCCTTGAAAAGTATATAAAACTCAATGAGGAATTGATAACAATA
```

Once you have the docker running, the kallisto.test.wdl file and the input file, run the following command

```{r}
java -jar wdltool/target/scala-2.11/wdltool-0.5.jar inputs kallisto.test.wdl > kallisto.test.json
```

The wdltool here creates a json file which defines the inputs, it includes theworkflow name and the inputs. It looks like below:

```{r}
{
  "kalWF.inkal": "File"
}

```

THE Json NEEDS TO BE MODIFIED and the input file specified (I used the vi editor):
```{r}
{
  "kalWF.inkal": "inkal.fasta"
}
```

Run the WDL
```{r}
java -jar cromwell/target/scala-2.11/cromwell-0.20.jar run kallisto.test.wdl kallisto.test.json
```

A successful run would have the following lines at the end of the run:
```{r}
[2016-07-12 11:01:55,21] [info] WorkflowManagerActor Workflow 989b788d-b7d2-4ef6-bce3-6d4c2372d17c succeeded!
[2016-07-12 11:01:55,21] [info] WorkflowManagerActor WorkflowActor-989b788d-b7d2-4ef6-bce3-6d4c2372d17c is in a terminal state: WorkflowSucceededState
{
  "outputs": {
    "kalWF.kallisto.response": "/Users/addepald/Documents/Docker/cromwell-executions/kalWF/989b788d-b7d2-4ef6-bce3-6d4c2372d17c/call-kallisto/stdout"
  },
  "id": "989b788d-b7d2-4ef6-bce3-6d4c2372d17c"
}
[2016-07-12 11:01:59,62] [info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
```

Kallisto index - produces a index file (with file name 'kalout.idx' specified in the command) in the output directory under /cromwell-executions mentioned inthe last few lines of the run, which can then be used to quantify abundances of the transcripts.
See [kallisto] {https://pachterlab.github.io/kallisto/} for furhter details.

---- More coming soon
