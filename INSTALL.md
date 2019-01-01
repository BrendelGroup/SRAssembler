# SRAssembler Installation and Setup

## Obtaining SRAssembler

Presumably you are reading this file on our github site and thus you are
likely to know that the following commands on your local machine should get
you going:

```bash
git clone https://github.com/BrendelGroup/SRAssembler
cd SRAssembler/
```

That said, an implicit assumption is that your local machine runs some version
of Linux.

## Prerequisite programs
__SRAssembler__ is a workflow that invokes easily available third-party software,
which must be pre-installed on your system.

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

To run the parallel version of __SRAssembler__, you will to have an MPI version such as [Open MPI](http://www.open-mpi.org/) installed.
All installed programs should of course be accessible in your binary search path (typically set in _~/.bashrc_ or _~/.profile_).
For Open MPI you may need to set the _LD_PRELOAD_ variable.
For example, in bash:
```bash
     export LD_PRELOAD=/usr/lib64/openmpi/lib/libmpi_cxx.so
```

The following is an executable record of installation on an Ubuntu system that should serve as a template for installation on any flavor of Linux.

```bash
#SYSTEM
	cd /usr/local/src
#... change the above if you want to install in a different location
#
	apt-get install openmpi-bin openmpi-doc libopenmpi-dev
	apt-get install libboost-all-dev

#ABYSS
	apt-get install abyss

#GENESEQER
	mkdir GENESEQER
	cd GENESEQER
	git clone https://github.com/brendelgroup/GeneSeqer
	cd GeneSeqer/
	cd src
	make linux
	make clean
	cp makefile.lnxMPI makefile.lnxMPIorig
	sed -e "s/^#MPICC/MPICC/" makefile.lnxMPI | sed -e "0,/^MPICC/s/^MPICC/#MPICC/" > makefile.lnxMPIu
	make -f makefile.lnxMPIu
	make clean
	make install
	cd ../../..

#GTH
	mkdir GTH
	cd GTH
	wget http://genomethreader.org/distributions/gth-1.6.6-Linux_x86_64-64bit.tar.gz
	tar -xzf gth-1.6.6-Linux_x86_64-64bit.tar.gz
	cd gth-1.6.6-Linux_x86_64-64bit
	cd bin
##note: To use gth, the user will have to get a licence file from http://genomethreader.org/cgi-bin/download.cgi
##      and either put the file into their home directory or replace the following symbolic link:
	ln -s ~/gth.lic ./
	cd ../../..

#SNAP
	mkdir SNAP
	cd SNAP
	wget http://korflab.ucdavis.edu/Software/snap-2013-11-29.tar.gz
	tar -xzf snap-2013-11-29.tar.gz
	cd ..

#SOAPdenovo2
	git clone https://github.com/aquaskyline/SOAPdenovo2
	cd SOAPdenovo2/
	make
	cp SOAP* /usr/local/bin
	make clean
	cd ..

#VMATCH
	mkdir VMATCH
	cd VMATCH
	wget http://www.vmatch.de/distributions/vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
	tar -xzf vmatch-2.3.0-Linux_x86_64-64bit.tar.gz
	ln -s vmatch-2.3.0-Linux_x86_64-64bit vmatch.distribution

	chown -R root:users vmatch.distribution
	chmod a-s vmatch.distribution
	chmod a+x vmatch.distribution
	chmod a+r vmatch.distribution

	cd vmatch.distribution
	chmod a+r chain2dim cleanpp.sh matchcluster mkdna6idx mkvtree repfind.pl upgradeprj.pl vendian vmatch vmatchselect Vmatchtrans.pl vmigrate.sh vseqinfo vseqselect vstree2tex vsubseqselect
	chmod a+x chain2dim cleanpp.sh matchcluster mkdna6idx mkvtree repfind.pl upgradeprj.pl vendian vmatch vmatchselect Vmatchtrans.pl vmigrate.sh vseqinfo vseqselect vstree2tex vsubseqselect
	chmod a+r *pdf README.distrib SELECT TRANS
	chmod a+x SELECT TRANS
	cd SELECT
	chmod a+r *
	cd ../TRANS
	chmod a+r *
	cd ..

	\cp chain2dim cleanpp.sh matchcluster mkdna6idx mkvtree repfind.pl vendian vmatch vmatchselect vmigrate.sh vseqinfo vseqselect vstree2tex vsubseqselect  /usr/local/bin/
	cd ../..
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
# ... will put binary SRAssemblerMPI into the ../bin folder
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
go to the _demo_ directory and study and execute the _xtest_ script.
