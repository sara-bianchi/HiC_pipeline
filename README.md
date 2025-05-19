# HiC_pipeline
Automated pipeline for chromatin conformation data analysis

# Installation and configuration
All the packages required to run the pipeline have already been installed within a Docker image available on our repository as hedgelab/bulk_image:image15.
The use of Docker is recommended to ensure reproducibility of the results. In case you don't have Docker installed on your computer, you can install it following the instructions: https://docs.docker.com/get-started/get-docker/
Otherwise, you can recreate the working environment by following the installation commands in the provided Dockerfile (optimized for a bash shell).

Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run only the analysis of already aligned matrices, 30 GB of RAM will be enough. Else, if you need to perform either indexing or alignment, depending on the size of your data, you could need between 50 and 100 GB. To set the proper RAM limits for your analysis, some versions of Docker Desktop allow you to directly edit these parameters in the settings. In contrast, others require the use of an external settings file. In the latter case, create a .wslconfig file in the User/user_name folder of your computer. You can follow the example template provided. After making these changes, restart the Docker engine to apply them effectively. 

Once Docker is installed and properly configured, activate the Docker engine by directly opening the Docker Desktop application or by running in your terminal:
start docker

# Preparation of the input files
The first step is to choose a working directory on your computer. This directory will be the one shared with the Docker container, so anything needed for the analysis should be located here, either directly or by organizing your input files in subdirectories. The final output will be written in the same directory. The container will be connected to this directory and will recognise it from a default internal path assigned when the container is created, which is /home/shared_folder. This means that to make your files recognisable within the container, all the subsequently required paths need to start with the same prefix /home/shared_folder/followed by the internal paths within your working directory. For instance, if a fastq file is saved in D:/Data/Working_dir/Fastq/sample1.fastq.gz, and I decide to use Working_dir as my working directory, the path of this file will be /home/shared_folder/Fastq/sample1.fastq.gz.

In your working directory, you need to have:
- your input files, either zipped fastq raw files, aligned matrices, and other associated files if required (.bed or .biases)
- genome information:
    -  the FASTA (.fa or .fasta) file if you need to perform the indexing. These files can be easily downloaded from the UCSC site: for the human genome, the latest version released is hg38, and you can directly download the FASTA at this link https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz. If you wish to use hg19, you can download the FASTA at this link https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz.
    -  the genome folder if you already have an indexed genome and you need to perform the alignment (the use of an already indexed genome is possible only in the case it has been created with Bowtie-2 for HiC-pro alignment and BWA for Arima alignment)
    -  the chrom.sizes file, which is needed for all the steps of the analysis. It can be downloaded from the UCSC site and must correspond to the same version of the genome used for the alignment index. hg19 and hg38 files for all the human chromosomes are provided; you can remove the lines related to some chromosomes that you do not want to include in your analysis (ex: chrY for female samples)
- the table.xlsx file with all the information about your samples. You can find two examples of this file in the provided folders. The name of this file must always be table.xlsx and must be stored directly in the working directory (not in subfolders). This is important since the pipeline directly looks for this file as /home/shared_folder/table.xlsx. This file is divided into different columns, each with a specific name needed to recognise the corresponding information:
    - sample_ID: contains a unique identifier of each sample
    - condition: contains the information related to the biological condition of each sample (ex: control, treated). Two samples in the same condition must have the same name in this column
    - batch: biological replicate identifier, needed to correct, if specified, for the batch effect.
    - fastq1: path of fastq files read1 within the shared folder, if you want to start from the alignment
    - fastq2: path of fastq files read2 within the shared folder, if you want to start from the alignment
    - matrix: path of the matrix file within the shared folder, if you have your data already aligned. It can be provided either in ICE normalised or raw format, and it is necessary to specify this feature in the settings.xlsx file (normalized = TRUE/FALSE). The matrix is accepted in HiC-Pro .matrix format, .hic, .homer, .h5, .cool, and .mcool formats. It is not strictly necessary that all the samples are provided in the same format unless the HiC-Pro format is present in at least one sample.
    - bed: if the matrix is provided in the .matrix HiC-pro format, it is necessary to associate the corresponding bed file
    - bias: if the matrix is provided in ICE normalized format, it is possible to add the file with the normalization biases. This file is directly provided with the HiC_pro alignment together with the normalised matrix, and it is made up of one column with the bin ID and one column with the bias value for that bin. Not providing this file for a normalized input will result in skipping the differential loops analysis, which specifically requires this file
- the settings.xlsx file, which specifies the steps of the analysis that you want to perform and the parameters needed to run each step. Two examples of this file are provided. Similarly to the table.xlsx, also for this file the path must be fixed as /home/shared_folder/settings.xlsx, meaning that you cannot change the name of this file or move it to subdirectories. This file is made up of at least three columns:
    - the first contains a unique identifier of the parameters, which must never be changed
    - the second contains a brief explanation about how to compile this parameter
    - the other ones contain the information needed for this specific parameter and can be edited. All the information must be added in order without leaving intermediate empty columns: for example, if only one info is needed (ex: alignment = TRUE/FALSE), the first column must be used, if more you should start to compile from the first column, then the second, and so on. The parameters that could need more columns are binning, reference, res_untargeted, res_v4C, res_transC, viewpoint, and seeds, which require inserting a list of bin sizes or chromosomal loci (one for each column). In this file, it is also required to add the paths of the genome information files. All the possible parameters are provided, but once you have chosen the TRUE/FALSE parameters to selectively run only some specific steps of your analysis, all the parameters related to deselected steps can be left empty (they will just be ignored by the pipeline). For example, if you do not need to run the alignment step and you set it to FALSE, all the parameters related to this step, like aligner, trimming, binning, and RE, can be left empty.
- Rscript.R, which you can directly download in its latest version from this GitHub (you can find it in the docker_creation folder)

# Running the pipeline
Once all the files have been put in the working directory and all the paths have been indicated in the settings and table files, open your terminal to run the analysis.
First, you have to create the connection between your folder and the container, and to do this, you have to run the following command:

docker run -d -v [Path/To/Your/Folder]:/home/shared_folder --name [container_name] hedgelab/hic_pipeline:image15

What you need to edit is everything between squared brackets [], adding the path of your folder and assigning a unique name to your container. If you have not downloaded the image yet on your local Docker repository, it will be directly downloaded from the Docker Hub. Once this download has been done, the creation of further containers for other analyses will be faster.
The next step consists of launching an Rscript that is saved in your shared folder, and it is also provided on this GitHub in case you want to edit and replace it. This script can perform all the steps of the analysis indicated in the settings file using the files provided in the table. It can then save all the outputs within a folder named Output_yymmdd in your working directory. To run this script, you just have to run the following command:

docker exec -it [container_name] Rscript /home/shared_folder/Rscript.R

The only part that you need to edit is the container_name which should match the one assigned in the previous step.

# Output interpretation
In your working directory it will be created a new folder system:
output_yymmdd
- bwa_index / bowtie_index: indexed genome if indexing is performed within the pipeline (BWA for Arima and Bowtie-2 for HiC-Pro)
- trimmed_fastq: trimmed FASTQ files if trimming has been performed
- HiC-pro: HiC pro intermediate output, statistics, and final matrices are allignment is performed with HiC-Pro
- Arima: Arima intermediate output, statistics, and final matrices are aligned with Arima
- cool_raw / cool_iced / h5_iced / mcool_raw / mcool_iced: folder with raw or normalized matrices in different formats based on the starting input and types of analyses required, each at the different binning sizes required in the settings
- ICE_normalisation: HiC-pro output obtained after ICE normalization of the provided raw matrices
- targeted_analysis:
    - virtual_4C: bedGraph tracks of the virtual 4C interaction score for each viewpoint (calculated with the Bioconductor HiContacts package) and .pdf plot of the 10 upstream and downstream bins from the viewpoint at different zoom levels (100 and 1000 maximal interaction values)
    - transC: output of trans interacting network analysis (https://github.com/Noble-Lab/trans-C)
        - transC normalized matrices (transC_matrix)
        - raw bin score values (bin_acores_all)
        - .bedGraph files of bin scores of trans chromosomes only (not the same chromosomes as the loci in the initial clique)
        - .pdf plot of trans bin scores interaction (bin_scores_trans)
        - .pdf plot of the bin score density distribution with the 1%, 5%, and 10% thresholds in red (density_plot)
        - .pdf circular plot with the top-ranked interactions from each of the starting seeds to the corresponding bin, colored by percentile classification (1% is red, 1-5% is purple, and 5-10% is light blue)
- untargeted_analysis:
    - autocorrelation: ct_perc plots with cis/trans interaction distribution across all the chromosomes and interaction map + autocorrelation matrix for each chromosome in each sample (HiContacts Bioconductor package)
    - correlation: Pearson correlation and PCA of the eigenvector values calculated for each sample
    - dcHiC: dcHiC output for each of the possible comparisons with the reference condition indicated in the settings.xlsx file (PCA calculation, best PCA selection, differential compartments and subcompartments analysis, B to A genes GSEA, fitHiC loops identification, differential loops analysis, and IGV summary)
    - diff_loops: differential loops analysis performed with limma on fitHiC output after filtering for significant loops over the background in all the samples of the same condition. The diff_loops files contain the limma output with logFC and adjusted p-value for each loop in each comparison.
    - genom_tracks: pyGenometracks for each chromosome for the individual samples (compartments, HiC matrix and TADs pattern), and for each comparison over the reference (compartments of the 2 conditions, mahalanobis distance and log10 adjusted-p-value of different compartments, subcompartments for the two conditions and log10 p-adjusted p-value, significantly identified loops in each condition and significantly identified loops in either conditions with width proportional to the -log10 adjusted p-value)
    - gsea: gsea tables and enrichment plots for each ontology for each comparison of the genes in significantly different B to A compartments
    - HiCrep: SCC score (stratified correlation coefficient) among the total HiC matrices and each chromosome. The scc_score_all plot shows a summary of the distribution of the SCC score in all the chromosomes for all the possible pairs of samples
    - IGV:
          - IGV .bedGraph tracks of the compartments and subcompartments for each sample, and log10 p-adjusted of the compartments and subcompartments for each pair of conditions compared
          - IGV .bedpe tracks of significantly different loops and log10 p-adjusted for each comparison
    - ratio: for each chromosome side to side and ratiometric matrices of each pair of samples of the same batch and different conditions (based on the selected references)
    - TADs: TADs boundaries, score, and domain classification for each sample


