bootstrap: docker
From: fedora:30

%help
    This container provides portable & reproducible components for SRAssembler:
    Selective Recursive local Assembler from the Brendel Group.
    Please see https://github.com/BrendelGroup/SRAssembler for complete documentation.
    To run single-threaded, use `singularity run SRAssembler.simg` with appropriate arguments.
    To run multi-threaded, use `singularity exec SRAssembler.simg mpirun -np $NUMBER_OF_PROCESSORS SRAssembler_MPI` with appropriate arguments.
    
%post
    dnf -y update
    dnf -y install @development-tools
    dnf -y install gcc-c++
    dnf -y install bc git tcsh tzdata unzip zip wget which bzip2
    dnf -y install nano
    dnf -y install python2-setuptools
    dnf -y install python3-setuptools

    echo 'Installing Abyss assembler version 2.1.5'
    #### Prerequisites
    dnf -y install boost-devel
    dnf -y install autoconf automake
    dnf -y install sparsehash-devel
    dnf -y install openmpi openmpi-devel
    export PATH=$PATH:/usr/lib64/openmpi/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
    #### Install
    cd /opt
    git clone --branch 2.1.5 --depth 1 https://github.com/bcgsc/abyss.git
    cd abyss
    ./autogen.sh
    ./configure
    make
    
    echo 'Installing SOAPdenovo2 assembler version r241'
    #### Prerequisites
    #### Install
    cd /opt
    git clone --branch r241 --depth 1 https://github.com/aquaskyline/SOAPdenovo2.git
    cd SOAPdenovo2
    make
    
    echo 'Installing VMATCH aligner version 2.3.0 from http://vmatch.de '
    #### Prerequisites
    #### Install
    cd /opt
    wget http://vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
    tar -xzf vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
    rm vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
    ln -s vmatch-2.3.0-Linux_x86_64-64bit VMATCH
    
    echo 'Installing BLAST+ version 2.9.0 from NCBI '
    #### Prerequisites
    #### Install
    cd /opt
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
    tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz
    rm ncbi-blast-2.9.0+-x64-linux.tar.gz
    ln -s ncbi-blast-2.9.0+ BLAST

    echo 'Installing GeneSeqer spliced aligner '
    #### Prerequisites
    #### Install
    cd /opt
    git clone https://github.com/BrendelGroup/GeneSeqer.git
    cd GeneSeqer/src
    make linux
    
    echo 'Installing GenomeThreader version 1.7.0 spliced aligner '
    #### Prerequisites
    #### Install
    cd /opt
    wget http://genomethreader.org/distributions/gth-1.7.0-Linux_x86_64-64bit.tar.gz
    tar -xzf gth-1.7.0-Linux_x86_64-64bit.tar.gz
    rm gth-1.7.0-Linux_x86_64-64bit.tar.gz
    ln -s gth-1.7.0-Linux_x86_64-64bit GENOMETHREADER

    echo 'Installing SNAP gene finder '
    #### Prerequisites
    #### Install
    cd /opt
    git clone https://github.com/KorfLab/SNAP.git
    cd SNAP
    make
        
    echo 'Installing SRAssembler from https://github.com/BrendelGroup/SRAssembler.git '
    #### Prerequisites
    #### Install
    cd /opt
    git clone https://github.com/BrendelGroup/SRAssembler.git
    cd SRAssembler/src
    make
    make clean
    make mpi
    make clean

    echo 'Installing bioawk '
    dnf -y install byacc zlib-devel
    #### Install
    cd /opt
    wget https://github.com/lh3/bioawk/archive/master.zip
    unzip master.zip 
    cd bioawk-master
    make
    cd ..
    rm master.zip

    echo 'Installing SAMTOOLS version 1.9 '
    #### Prerequisites
    dnf -y install libncurses.so.5 ncurses ncurses-devel bzip2-devel xz-devel
    #### Install
    cd /opt
    wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
    tar -xf samtools-1.9.tar.bz2
    rm samtools-1.9.tar.bz2
    cd samtools-1.9
    ./configure
    make && make install

    echo 'Installing BOWTIE2 version 2.3.4.3 '
    #####
    cd /opt
    wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip
    unzip bowtie2-2.3.4.3-linux-x86_64.zip
    rm bowtie2-2.3.4.3-linux-x86_64.zip

    echo 'Installing Integrated Genome Viewer version 2.4.14'
    #### Dependencies
    dnf -y install java-1.8.0-openjdk xorg-x11-server-Xvfb
    ####
    cd /opt
    wget http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip
    unzip IGV_2.4.14.zip
    rm IGV_2.4.14.zip
    # Patch the igv script to run headless. No need for X display.
    printf '#!/bin/sh\n(Xvfb :10 &) && DISPLAY=:10 java -Xmx4000m -Dapple.laf.useScreenMenuBar=true -Djava.net.preferIPv4Stack=true -jar /opt/IGV_2.4.14/lib/igv.jar "$@" && pkill Xvfb\n' > /usr/local/bin/igv
    chmod 777 /usr/local/bin/igv

    echo 'Installing MUSCLE aligner version 3.8.31 '
    #### Dependencies
    ####
    cd /opt
    wget https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
    tar -xf muscle3.8.31_i86linux64.tar.gz
    rm muscle3.8.31_i86linux64.tar.gz
    ln -s muscle3.8.31_i86linux64 muscle

    echo 'Installing InterMine Python package '
    easy_install intermine
        
    echo 'Installing R '
    #### Dependencies
    dnf -y install udunits2-devel libcurl libcurl-devel geos-devel
    ####
    cd /opt
    dnf -y install R-base R-devel
    R CMD javareconf

    echo 'Installing R packages'
    echo 'repo <- "http://ftp.ussg.iu.edu/CRAN"' > R2install
    #### Dependencies
    dnf -y install cairo-devel libXt-devel
    ####
    echo 'install.packages("R.devices", repos = repo, quick = TRUE)' >> R2install
    #### Dependencies
    dnf -y install openssl-devel mysql-devel postgresql-devel 
    ####
    echo 'install.packages("dplyr", repos = repo, quick = TRUE)' >> R2install
    echo 'install.packages("ggplot2", repos = repo, quick = TRUE)' >> R2install
    #### Dependencies
    dnf -y install libjpeg-turbo-devel libxml2-devel ImageMagick-c++-devel
    ####
    echo 'install.packages("knitr", repos = repo, quick = TRUE)' >> R2install
    echo 'install.packages("readr", repos = repo, quick = TRUE)' >> R2install
    echo 'install.packages("BiocManager", version = "1.30.4", repos = repo, quick=TRUE)' >> R2install
    echo 'BiocManager::install(c("GenomicRanges"), version = c("3.8"))' >> R2install
    
    Rscript R2install



    
%environment
    export LC_ALL=C
    export PATH=$PATH:/opt/VMATCH
    export PATH=$PATH:/opt/SOAPdenovo2
    export PATH=$PATH:/opt/GeneSeqer/bin
    export PATH=$PATH:/opt/BLAST/bin
    export PATH=$PATH:/opt/GENOMETHREADER/bin
    export PATH=$PATH:/opt/SNAP
    export PATH=$PATH:/opt/abyss/bin
    export PATH=$PATH:/opt/abyss/ABYSS
    export PATH=$PATH:/opt/bioawk-master
    export PATH=$PATH:/opt/samtools
    export PATH=$PATH:/opt/bowtie2-2.3.4.3-linux-x86_64
    export PATH=$PATH:/usr/lib64/openmpi/bin
    export PATH=$PATH:/opt
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
    export BSSMDIR="/opt/GENOMETHREADER/bin/bssm"
    export GTHDATADIR="/opt/GENOMETHREADER/bin/gthdata"
    export ZOE="/opt/SNAP"
    export PATH=$PATH:/opt/SRAssembler/bin

%runscript
    exec SRAssembler "$@"

%test
    export LC_ALL=C
    export PATH=$PATH:/opt/VMATCH
    export PATH=$PATH:/opt/SOAPdenovo2
    export PATH=$PATH:/opt/GeneSeqer/bin
    export PATH=$PATH:/opt/BLAST/bin
    export PATH=$PATH:/opt/GENOMETHREADER/bin
    export PATH=$PATH:/opt/SNAP
    export PATH=$PATH:/opt/abyss/bin
    export PATH=$PATH:/opt/abyss/ABYSS
    export PATH=$PATH:/usr/lib64/openmpi/bin
    export PATH=$PATH:/opt
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
    export BSSMDIR="/opt/GENOMETHREADER/bin/bssm"
    export GTHDATADIR="/opt/GENOMETHREADER/bin/gthdata"
    export ZOE="/opt/SNAP"
    export PATH=$PATH:/opt/SRAssembler/bin
    cd /opt/SRAssembler/demo
    ./xtest defaultpath


%labels
    Maintainer vpbrendel
    Version v1.1.0
