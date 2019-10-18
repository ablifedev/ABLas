## Introduce

ABLas is a workflow for alternative splicing analysis of RNA-seq sequencing datasets.

Website： https://github.com/ablifedev/ABLas/

## Download

Download current release: [ABLas-v0.2](https://github.com/ablifedev/ABLas/archive/v0.2.tar.gz)

Download example data:    [ABLas-example-data](https://github.com/ablifedev/ABLas-example-data)

## Installation

OS Requirements
This package is supported for *Linux*. The package has been tested on the following systems:
+ Linux: Centos 6.5

Required software:

* Perl (>=5.10, https://www.perl.org/get.html)

* Python2 (>=2.7, https://www.python.org/downloads/)

* R (>=3.2.0, https://cloud.r-project.org/)

* Samtools (>=1.2, http://www.htslib.org/download/)

* TopHat2 (>=2.0.13, http://ccb.jhu.edu/software/tophat/index.shtml)

Required R packages:

* ggplot2  

* reshape2  

* RColorBrewer  

Required Perl packages:

* Bio::DB::Sam  

* GD::Simple  

* GD  

* YAML::XS 

After all of the requirments being installed, you can install the program as followed:   

1. Unzip the file `ABLas.tar.gz`:   
    `tar -zxf ABLas.tar.gz`   

## Usage

First time to run the program,you can by the fowlling steps:  

1. Step1:  
Copy the configuration file `as_config` to you workdir:  
    `cp as_config <workdir>` 
    
2. Step2:  
Fill in the configuration file correctly fill in the parameters according to the information in the configuration file, specific requirements will be detailed in the annotation. 

3. Step3:  
Run the ABLas : ou can run the ABLas like this:  
    `perl runABLas.pl -config as_config`

## Test data

In the `example/test_reads/`, there are junctions reads files of two RNA-seq data, ‘Test’ is the PTB-knockdown RNA-seq data in HeLa cells, and the other "Ctrl" is the control, you could test the ABLas with them or download raw data from GSE19323.
In the `example/report/`, there is a complete output of the ABLas program, it can be used as a reference after the program is run.  

## Config file template

You can copy and modify this config template.

```
# This is the config file to do AS for the RNA-Seq data

# the species
SPECIES human

# the gff file
GFF ~/ABLas-example-data/human.gff

# the fasta file
FASTA   human.fa

#whether to do the RAS step,if 0,exit;if 1 do the RAS step; default 1
RAS 1

#the in directory of the file
INDIR   ~/ABLas-example-data/

#the sample junctions file and sample sam file and sample name,the file and the name are separated by colon
SAMPLE1 Test/junctions.bed:Test/accepted_hits.uniq.sam:Test
SAMPLE2 Ctrl/junctions.bed:Ctrl/accepted_hits.uniq.sam:Ctrl

#the comparing group to do RAS, the group name is the sample file name,the group names are separated by colon
GROUP1  Test:Ctrl

#the out directory of the result,please input the absolute path(start with "/")

OUTDIR  ~/ABLas-example-data/result/

```

## Running time

<1 hour in 16 cores computer
