FROM ubuntu:24.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -y update && apt-get -y upgrade
RUN apt-get update \
  && apt-get install -y build-essential \
  wget \
  unzip \
  bzip2 \
  gcc \
  g++ && apt-get clean
RUN apt-get install xxd
RUN apt-get -y install make
RUN apt-get install zlib1g-dev

#Bioconductor
RUN apt-get -y update && apt-get -y upgrade && apt-get -y install r-base
RUN apt-get -y install gdebi-core
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2023.12.1-402-amd64.deb
RUN gdebi -n rstudio-server-2023.12.1-402-amd64.deb
RUN useradd rstudio -p "\$y\$j9T\$/.6YKeUOB4ifaPjuG/xaC1\$0162SW98NtTo5c6I7uXbwlNlKGuu9LTcUanCzz6DF/C" -d /home/rstudio -m
RUN apt-get update
RUN apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
RUN apt-get update && apt install -y libudunits2-dev libgdal-dev
RUN apt-get update
RUN apt-get -y install gfortran
RUN apt-get -y install build-essential
RUN apt-get -y install fort77
RUN apt-get -y install xorg-dev
RUN apt-get -y install liblzma-dev  libblas-dev gfortran
RUN apt-get -y install gobjc++
RUN apt-get -y install aptitude
RUN apt-get -y install libbz2-dev
RUN apt-get -y install libpcre3-dev
RUN aptitude -y install libreadline-dev
RUN apt-get -y install libcurl4-openssl-dev
RUN apt install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN apt-get install -y libgit2-dev
RUN apt-get install -y libharfbuzz-dev
RUN apt-get install -y libfribidi-dev
RUN apt-get install -y cmake
RUN apt-get install -y libcairo2-dev
RUN Rscript -e 'install.packages(c("ggplot2", "dplyr", "BiocManager", "tidyverse", "ggplotify", "stringr", "viridis", "ggthemes", "ggrepel", "igraph", "readxl", "pheatmap", "reshape2", "pals", "devtools"), dependencies = TRUE)'
RUN  Rscript -e 'BiocManager::install(c("limma", "Glimma", "edgeR", "scran", "fgsea", "DESeq2", "ensembldb", "topGO", "org.Hs.eg.db", "clusterProfiler", "enrichplot",  "AnnotationHub", "ADImpute", "pathview", "TCseq"), update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("HiCExperiment", update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("HiCool", update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("HiContacts", update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("HiContactsData", update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install("fourDNData", update = TRUE, ask = FALSE)'
RUN Rscript -e 'BiocManager::install(c("impute", "preprocessCore"), update = TRUE, ask = FALSE)'
RUN Rscript -e 'install.packages(c("WGCNA"), dependencies = TRUE)'
RUN Rscript -e 'install.packages(c("factoextra"), dependencies = TRUE)'

#HiC-pro
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -p /usr/local/anaconda
RUN rm ~/anaconda.sh
ENV PATH /usr/local/anaconda/bin:$PATH
COPY environment_noR.yml /
RUN conda env create -f /environment_noR.yml && conda clean -a
ENV PATH /usr/local/anaconda/envs/HiC-Pro_v3.1.0/bin:$PATH
RUN cd /tmp && \
    echo "master.zip" | wget --base=http://github.com/nservant/HiC-Pro/archive/ -i - -O hicpro_latest.zip && \
    unzip hicpro_latest.zip && \
    cd HiC-Pro-master  && \ 
    make configure prefix=/ && \
    make install && \
    cd .. && \
    rm -fr HiC-Pro*

#HiCexplorer
RUN wget https://repo.continuum.io/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh -O ~/anaconda.sh
RUN bash ~/anaconda.sh -b -u -p /usr/local/anaconda
RUN rm ~/anaconda.sh
ENV PATH /usr/local/anaconda/bin:$PATH
RUN conda install -y python
RUN conda install -y numpy
RUN conda install -y scipy
RUN conda install -y matplotlib-base
RUN conda install -y ipykernel
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda install -y pysam
RUN conda install -y intervaltree
RUN conda install -y biopython
RUN conda install -y pandas
RUN conda install -y pybigwig
RUN conda install -y jinja2
RUN conda install -y unidecode
RUN conda install -y hicmatrix
RUN conda install -y hic2cool
RUN conda install -y psutil
RUN conda install -y pygenometracks
RUN conda install -y fit_nbinom
RUN conda install -y cooler
RUN conda install -y krbalancing
RUN conda install -y pybedtools
RUN conda install -y future
RUN conda install -y tqdm
RUN conda install -y python-graphviz
RUN pip install hyperopt
RUN conda install hicexplorer

#dcHiC
RUN apt-get update
RUN apt-get install -y git
RUN git clone https://github.com/ay-lab/dcHiC
RUN conda env create -f /dcHiC/packages/dchic.yml
RUN Rscript -e 'install.packages(c("reticulate"), dependencies = TRUE)'
RUN conda init bash
RUN rm /dcHiC/dchicf.r
COPY dchicf.r /dcHiC/dchicf.r
COPY script_dcHiC.sh /script_dcHiC.sh
RUN Rscript -e 'library(reticulate); conda_run2(envname = "dchic", cmd_line = "/script_dcHiC.sh", echo = TRUE)'

#fitHiC
RUN rm /fithic/fithic/fithic.py
COPY fithic.py /fithic/fithic/fithic.py
COPY preprocess.py /preprocess.py

#transC
RUN apt-get update
RUN git clone https://github.com/Noble-Lab/trans-C.git
RUN conda install -y scikit-learn-intelex
RUN Rscript -e 'install.packages(c("circlize"), dependencies = TRUE)'
RUN Rscript -e 'install.packages(c("Rlab"), dependencies = TRUE)'

#genomeTracks
COPY individual_samples_track /individual_samples_track
COPY comparison_track /comparison_track
COPY 4C_track /4C_track

#HiC_rep
RUN Rscript -e 'library(devtools); devtools::install_github("TaoYang-dev/hicrep")'

#Arima
RUN git clone https://github.com/ArimaGenomics/mapping_pipeline.git
RUN chmod 777 /mapping_pipeline/filter_five_end.pl
RUN chmod 777 /mapping_pipeline/arima_mapping_pipeline.sh
RUN chmod 777 /mapping_pipeline/get_stats.pl
RUN chmod 777 /mapping_pipeline/two_read_bam_combiner.pl

#bwa installation
RUN apt-get update && apt-get install -y git && git clone https://github.com/lh3/bwa.git && cd bwa && make

#java
COPY jdk-17.0.12_linux-x64_bin.tar.gz /home/jdk-17.0.12_linux-x64_bin.tar.gz
RUN mkdir -p /usr/lib/jvm
RUN tar zxvf /home/jdk-17.0.12_linux-x64_bin.tar.gz -C /usr/lib/jvm
RUN update-alternatives --install "/usr/bin/java" "java" "/usr/lib/jvm/jdk-17.0.12/bin/java" 1
RUN update-alternatives --set java /usr/lib/jvm/jdk-17.0.12/bin/java

#picard-tools
RUN git clone https://github.com/broadinstitute/picard.git
WORKDIR /picard
RUN ./gradlew shadowJar
WORKDIR /
RUN chmod 777 /picard/build/libs/picard.jar

#joicer tools
COPY juicer_tools_1.13.02.jar /home/juicer_tools_1.13.02.jar

#pipeline
COPY base.py /trans-C/code/base.py
COPY trans_C.py /trans-C/code/trans_C.py
COPY Rscript.R /Rscript.R

EXPOSE 8787
ENTRYPOINT ["tail"]
CMD ["-f","/dev/null"]