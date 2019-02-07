# SRAssembler Installation and Setup

## Obtaining SRAssembler

Presumably you are reading this file on our github site and thus you are likely to know that the following commands on your local machine should get you going:

```bash
git clone https://github.com/BrendelGroup/SRAssembler
cd SRAssembler/
```

That said, an implicit assumption is that your local machine runs some version of Linux.

If you would like to avoid the hassle of installing and configuring a bunch of dependencies, we highly recommend using the containerization software [Singularity](https://www.sylabs.io/singularity/) and running `singularity pull shub://BrendelGroup/SRAssembler` to download a convenient SRAssembler container.

## Prerequisite programs
__SRAssembler__ is a workflow that invokes easily available third-party software, which must be pre-installed on your system.

- string matching and read mapping:
  [Vmatch](http://www.vmatch.de)

- at least one short read assembler:
  [SOAPdenovo2](http://soap.genomics.org.cn/soapdenovo.html)
  [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss)

- at least one program for spliced-alignment:
   [GeneSeqer](https://github.com/BrendelGroup/GeneSeqer)
   [Genomethreader](http://www.genomethreader.org)

- (optional) a gene finder:
  [SNAP](http://korflab.ucdavis.edu/software.html)

All installed programs must be accessible in your binary search path ($PATH, typically set in _~/.bashrc_ or _~/.profile_).
To run the parallel version of __SRAssembler__, you will need to have an MPI version such as [Open MPI](http://www.open-mpi.org/) installed.
You may find you need to set the LD_LIBRARY_PATH for proper Open MPI function:
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib
```

The following is an executable record of installation on an Ubuntu system that should serve as a template for installation on any flavor of Linux.

```bash
# Change this if you want to install dependencies in a different location.
source_location="/usr/local/src"
# Change this if you want the executable binaries in a different location.
binary_location="/usr/local/bin"
pushd ${source_location}
	
#SYSTEM
	apt-get install openmpi-bin openmpi-doc libopenmpi-dev
	apt-get install libboost-all-dev

#ABYSS
	apt-get install abyss
	
#DUSTMASKER
	wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz
	tar -xzf ncbi-blast-2.8.1+-x64-linux.tar.gz
	rm ncbi-blast-2.8.1+-x64-linux.tar.gz
	ln -i -t "${binary_location}" ncbi-blast-2.8.1+/bin/*

#GENESEQER
	git clone https://github.com/brendelgroup/GeneSeqer
	pushd GeneSeqer/src
		make linux
	popd
	ln -i GeneSeqer/bin/GeneSeqer "${binary_location}"

#GENOMETHREADER
	wget http://genomethreader.org/distributions/gth-1.7.1-Linux_x86_64-64bit.tar.gz
	tar -xzf gth-1.7.1-Linux_x86_64-64bit.tar.gz
	rm gth-1.7.1-Linux_x86_64-64bit.tar.gz
	ln -i gth-1.7.1-Linux_x86_64-64bit/bin/gth "${binary_location}"
	# GenomeThreader requires its data to be in the same directory as the binary.
	cp -ir gth-1.7.1-Linux_x86_64-64bit/bin/gthdata "${binary_location}"
	cp -ir gth-1.7.1-Linux_x86_64-64bit/bin/bssm "${binary_location}"

#SNAP
	git clone https://github.com/KorfLab/SNAP.git
	pushd SNAP
		make
	popd
	ln -i SNAP/snap "${binary_location}"
	# SNAP requires the ZOE environment variable (see its manual).
	# If you intend to use SNAP, adjust the following line
	# and add it to your .bashrc, .bash_profile, or .profile file.
	export ZOE="${source_location}/SNAP"

#SOAPdenovo2
	git clone https://github.com/aquaskyline/SOAPdenovo2
	pushd SOAPdenovo2/
		make
	popd
	ln -i -t "${binary_location}" SOAPdenovo2/SOAPdenovo*

#VMATCH
	wget http://www.vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
	tar -xzf vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
	rm vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
	pushd vmatch-2.3.0-Linux_x86_64-64bit
		ln -i -t "${binary_location}" chain2dim cleanpp.sh matchcluster mkdna6idx mkvtree repfind.pl vendian vmatch vmatchselect vmigrate.sh vseqinfo vseqselect vstree2tex vsubseqselect
	popd
	
popd
```

## Serial version installation
```bash
cd src
make
# ... will put binary SRAssembler into the ../bin folder
make install
# ... will put the contents of ../bin into /usr/local/bin [default]
```

## Parallel version installation
```bash
cd src
make mpi
# ... will put binary SRAssembler_MPI into the ../bin folder
make install
# ... will put the contents of ../bin into /usr/local/bin [default]
```

## Installation of both versions
```bash
cd src
make
make clean
make mpi
make install
```

## C++ Boost libraries option
In most Linux systems, the C++ Boost libraries are already installed.
If this is not the case on your system, then please see about installation at [boost.org](http://www.boost.org/).
If the boost header files are not installed in the usual _/usr/include/boost_ directory, then use the _with-boost_ option to make to specify the location:

```bash
make with-boost=boost_path
#or
make mpi with-boost=boost_path
```

## Finally
Go to the _demo_ directory and examine and execute the _xtest_ script.
