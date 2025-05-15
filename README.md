# HiC_pipeline
Automated pipline for chromatin conformation data analysis

# Installation and configuration
All the packages required to run the pipeline have already been installed within a Docker image available on our repository as hedgelab/bulk_image:image15.
The use of Docker is recommended to ensure reproducibility of the results. In case you don't have Docker installed on your computer, you can install it following the instructions: https://docs.docker.com/get-started/get-docker/
Otherwise, it is possible to recreate the working environment following the installation commands in the provided Dockerfile (optimized for a bash shell).

Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run only the analysis of already aligned matrices, 30 GB of RAM will be enough. Else, if you need to perform either indexing or alignment, depending on the size of your data you could need between 50 and 100 GB. To set the proper RAM limits for your analysis, some versions of Docker Desktop allow you to directly edit these parameters in the settings, while others require the use of an external settings file. In the latter case, create a .wslconfig file in the User/user_name folder of your computer. You can follow the example template provided. After making these changes, restart the Docker engine to apply them effectively. 

Once Docker is installed and properly configured, activate the Docker engine directly opening the Docker Desktop application or by running on your terminal:
start docker

# Preparation of the input files
The first step is to choose a working directory on your computer. This directory will be the one shared with the Docker container, so anything needed for the analysis needs to be located here, either directly or by organizing your input files in subdirectories. Also the final output will be written in the same directory. The container will be connected to this directory and will recognise it from a default internal path assigned when the container is created, which is /home/shared_folder. This means that to make your files recognisable within the container, all the subsequently required paths need to start with the same prefix /home/shared_folder/ folowed by the internal paths within your working directory. For instance, if a fastq file is saved in D:/Data/Working_dir/Fastq/sample1.fastq.gz, and I decide to use Working_dir as my working directory, the path of this file will be /home/shared_folder/Fastq/sample1.fastq.gz.

In your working directory you need to have:
- your input files, either zipped fastq raw files, aligned matrices and other associated files if required (.bed or .biases)
- genome information:
    -  the FASTA (.fa or .fasta) file if you need to perform the indexing. This files can be easily downloaded from the UCSC site: for the human genome, the latest version released is hg38, and you can directly download the FASTA at this link https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz. if you wish to use hg19 you can download the FASTA at this link https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz.
    -  the genome folder if you already have an indexed genome and you need to perform the alignment (the use of an already indexed genome is possible only in the case it has been created with Bowtie-2 for HiC-pro alignment and BWA for Arima alignment)
    -  the chrom.sizes file, which is needed for all the steps of the analysis. It can be downloaded from the UCSC site and must correspond to the same version of the genome used for the alignment index. hg19 and hg38 files for all the human chromosomes are provided, you can remove the lines related to some chromosomes which you do not want to ionclude in your analysis (ex: chrY for female samples)
- the table.xlsx file with all the information of your samples. You can find two examples of this file in the provided folders. The name of this file must always be table.xlsx and must be stored directly in the working directory (not in subfolders). This is important since the pipeline directly llok for this file as /home/shared_folder/table.xlsx. This file is devided in different columns, each with a specific name needed to recognise the corresponding information:
    - sample_ID: contains a unique identifier of each sample
    - condition: contains the information related to the biological condition of each sample (ex: control, treated). Two samples in the same condition must have the same name in this column
    - batch: biological replicate identifier, needed to correct, if specified, for the batch effect.
    - fastq1: path of fastq files read1 within the shared folder, if you want to start from the alignment
    - fastq2: path of fastq files read2 within the shared folder, if you want to start from the alignment
    - matrix: path of the matrix file within the shared folder if you have your data already aligned. It can be provided either in ICE normalised or raw format and it is necessary to specify this feature in the settings.xlsx file (normalized = TRUE/FALSE). The matrix is accepted in HiC-Pro .matrix format, .hic, .homer, .h5, .cool and .mcool formats. it is not strictly necessary that all the samples are provided in the same format unless HiC-Pro format is present in at least one sample.
    - bed: if the matrix is provided in the .matrix HiC-pro format, it is necessart to associate the corresponding bed file
    - bias: if the matrix is provided in ICE normalized format it is possible to add the file with the normalisation biases. This file is directly provided with the HiC_pro alignment otgether with the normalised matrix, and it is made up of one column with the bin ID and one column with the bias value for that bin. Not providing this file for a normalized input will result in skipping the differential loops analysis which specifically requires this file
- the settings.xlsx file, which specifies the steps of the analysis that you want to perform and the parameters needed to run each step. Two examples of this file are provided. Similarly to the table.xlsx also for this file the path must be fixed as /home/shared_folder/settings.xlsx, meaning that you cannot change the name of this file or move it in subdirectories. This file is made up of at least three columns:
    - the first contains a unique identifier of the parameters which must never be changed
    - the second contains a brief explenation about how to compile this parameter
    - the other ones contain the information needed for this specific parameter and can be edited. All the information must be added in order without leaving intermediate empty columns: for example, if only one info is needed (ex: alignment = TRUE/FALSE) the first column must be used, if more you should start top compile from the first column, then the second, and so on. The parameters which could need more columns are binning, res_untargeted, res_v4C, res_transC and seeds, which require to insert a list of bin sizes or chromosomal loci (one for each column). In this file it is also required to add the paths of the genome information files. All the possible parameters are provided, but once you have chosen with the TRUE/FALSE parameters to selectively run only some specific steps of your analysis, all the parameters related to deselected steps can be left empty (they will just be ignored by the pipeline). For example if you do not need to run the alignment step and you set it to FALSE, all the parameters related to this step, like aligner, trimming, binning and RE can be left empty.
- Rscript.R which you can directly download in its latest version from this GitHub (you can find it in the docker_creation folder)

# Running the pipeline
Once all the files have been put in the working directory and all the paths have been indicated in the settings and table files, open your terminal to run the analysis.
First you have to create the connection between your folder and the container, and to do this you have to run the following command:

docker run -d -v [Path/To/Your/Folder]:/home/shared_folder --name [container_name] hedgelab/hic_pipeline:image15

what you need to edit is everything between squared brackets [], adding the path of your folder and assigning a unique name to your container. If you have not downloaded the image yet on your local Docker repository it will be directly downloaded from the Docker Hub. Once this download has been done the creation of futher containers for other analyses will be faster.
The next step consists in launching an Rscript which is saved in your shared folder, and it is also provided in the GitHub in case you want to edit and replace it. This script can performe all the steps of the analysis indicated in the settings file using the files provided in the table. It can then save all the outputs within a folder named Output_yymmdd in your working directory. To run this script you just have to run the following command:

docker exec -it [container_name] Rscript /home/shared_folder/Rscript.R

The only part that you need to edit is the container_name which should match the one assigned in the previous step.

# Output interpretation
In your working directory it will be created a new folder system:
Output_yymmdd


