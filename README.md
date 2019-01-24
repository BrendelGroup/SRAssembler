# SRAssembler Manual

The **SRAssembler** repository is the distribution point of our code for selective
local assembly of genomic regions matching a homologous query protein or DNA sequence.

SRAssembler (Selective Recursive local Assembler) is a modular pipeline program that can assemble genomic DNA reads into contigs that are homologous to a query DNA or protein sequence. SRAssembler first collects the reads that can be locally aligned to the query sequence and assembles them into contigs. Additional reads are then found by aligning against the contigs assembled in the previous round. This _in silico_ chromosome walking strategy recursively gathers reads that are associated with regions of interest. The gene structures of the final contigs are predicted by spliced alignment and _ab initio_ gene finding programs.

Hopefully this code will be useful to you if you 
- have a set of genomic reads for some organism and want to find out whether that set of reads represents deep enough sequencing to assemble genomic regions containing or encoding a homolog to a known gene or protein you are providing as query, without having to completely assemble the genome first.
- have sequencing reads for several organisms and want to quickly assemble contigs containing the homologs of a query gene in each organism, without having to completely assemble each genome first.
- have a set of genomic reads for some organism and want to identify several potentially homologous genes of one particular query, without having to completely assemble the genome first.

## Installation

Detailed instructions for installation from source are available in the [INSTALL](./INSTALL.md) document. You can also download and use our [Singularity image](https://singularity-hub.org/collections/1653) directly.

## Running SRAssembler

SRAssembler is run at the command line with `path/to/SRAssembler [options]`. It can also take advantage of multiple processors on machines with OpenMPI by using `mpirun -np $number_of_processors path/to/SRAssembler_MPI [options]`.

If you are using the Singularity image instead of compiling the SRAssembler binaries, you can run the SRAssembler program with `singularity run --cleanenv --bind $(pwd) path/to/SRAssembler.simg [options]`. To use the Singularity image with multiple processors, run `singularity exec --cleanenv --bind $(pwd) path/to/SRAssembler.simg mpirun -np $number_of_processors SRAssembler_MPI [options]`.

### Command-line options

SRAssembler is a workflow program that controls and directs the activity of other programs (e.g., a read aligner, an assembler, a spliced aligner, an _ab initio_ gene finder). SRAssembler requires both a parameter configuration file and command-line arguments. The parameter file contains settings for the programs under SRAssembler's control, while the command-line options have broader effects on SRAssembler's behavior.

**-q** *query_file* : FASTA file containing the DNA or protein query sequence.
Required unless only pre-processing reads (-P).

**-t** *query_type* : Query sequence type.
Accepts "protein" or "dna" (default: "protein").

**-p** *parameter_file* : Configuration file for the programs used by SRAssembler.
SRAssembler allows user to specify the parameters of the programs it controls. The parameter configuration file has sections for each program, marked by a header line in square brackets. In some cases, such as the aligner Vmatch, there are multiple headers for a program, distinguishing different parameters for different usages. The example [parameter file](demo/SRAssembler.conf) contains all of the program sections that SRAssembler recognizes, even though not all of them will be used in any one run of SRAssembler (e.g., GeneSeqer and GenomeThreader are mutually exclusive). The example file also contains comment lines marked with '#' that explain what each section is for. Copying and modifying this example file for your own runs is encouraged. A parameter file is required unless only pre-processing reads (-P flag).

**-o** *output\_directory* : Output directory (default: "./SRAssembler_output").

**-T** *temp_directory* : Directory for temporary files.
SRAssembler writes and reads many small files that are removed at the end of a run. Default is "/tmp", but you may get speed enhancements by using "/dev/shm" to keep these files entirely in memory.

#### Read libraries

**-l** *library_file* : Library definition file for sequencing reads.
Required if the **-1** option is not used. Paired-end reads must be in two sorted files rather than one interleaved file. Each library is specified in its own section, headed by a "\[LIBRARY\]" line. See the [demo](./demo) directory for [example](./demo/libraries_200bp.conf) [files](./demo/libraries_200bp_1kb.conf). Each library includes the following items, although default values will be used for some items if they are absent: 
- library_name: The name used for the directory containing the processed reads of this library. If unspecified, will default to the basename of the left reads file.
- insert_size: The insert size of this sequencing experiment. Only used for paired-end reads. Defaults to "300".
- direction: The sequencing direction of paired-end reads. Use "0" for forward-reverse, "1" for reverse-forward. Defaults to "0".
- r1: The left-end reads file or the single-end reads file. Specify the file path relative to the directory in which you will run SRAssembler, or use an absolute path.
- r2: The right-end reads file. Do not specify if this library is single-end.
- format: "fastq" or "fasta". Defaults to "fastq".

**-1** *left\_end\_file* : The single-end reads file, or the left-end reads file for paired-end reads.
Required if **-l** option is not used.

**-2** *right\_end\_file* : The right-end reads file for paired-end reads.

**-Z** *insert_size* : The insert size of the paired-end reads specified by **-1** and **-2**, required if **-2** option is used.

**-r** *preprocessed_directory* : Directory in which to store or from which to retrieve processed reads.
Before the recursive search and assembly of SRAssembler, the input read files specified for each Library are split into multiple files and indexed by the alignment program, Vmatch. This pre-processing step only needs to happen once. When SRAssembler is run multiple times, pre-processed reads in this directory that correspond to the input libraries will be used again without needing to be processed.
If not specified, defaults to "*output\_directory*/processed\_reads".

**-R** *num\_reads\_per\_split* : Number of reads in each split file.
While using larger split file parts tends to be faster, there should be enough parts to take advantage of all the processors available if you are using MPI parallelization. Note that if pre-processed reads already exist for a library, they will not be re-processed into new sizes due to this option. Defaults to "500000".

**-P** : Run the read processing step only.

#### Sub-process programs

**-A** *assembler_program* : The read assembler program used by SRAssembler.
Use "0" for SOAPdenovo2 or "1" for ABySS, defaults to "0".

**-k** *kmer_sizes* : The k-mer sizes used by the assembler.
SRAssembler will test multiple k-mer sizes and select the best for each round based on spliced-alignment results. The format for specifying which k-mers to test is *start\_k*:*interval*:*end\_k*. The *start\_k* and *end\_k* values must be odd, and the *interval* value must be even. For example, the default of "15:10:45" means k-mer values 15, 25, 35, and 45 will be tested.

**-S** *spliced\_alignment\_program* : The spliced alignment program used by SRAssembler.
Use "0" for GeneSeqer or "1" for GenomeThreader, defaults to "1".

**-s** *species* : Species model for spliced alignment.
Options are "human", "arabidopsis", "mouse", "rat", "chicken", "drosophila", "nematode", "fission_yeast", "aspergillus", "maize", "rice", and "medicago". Defaults to "arabidopsis".

**-G** *gene\_finding\_program* : The *ab initio* gene finding program used by SRAssembler.
Use "0" to not use a gene finder or "1" for SNAP, defaults to "0".

#### Completion thresholds

**-i** *query\_contig\_minimum* : Minimum length of contigs used for finding more reads during chromosome walking.
Increasing this may make SRAssembler searches slower, but more stringent. Defaults to "200".

**-m** *minimum\_contig\_length* : Minimum length of contig to accept as a potential hit for the query.
Defaults to "200".

**-M** *maximum\_contig\_length* : Maximum length of contig to assemble before discontinuing extension of that contig.
If an assembled contig is larger than the *maximum\_contig\_length*, it is trimmed at both ends to equal the *maximum\_contig\_length* and copied to the candidate-long-contig file. The candidate-long-contigs are aligned to the assembled contigs of the following round, and if they match (confirming a correct assembly) they are moved to the permanent long-contigs file and removed from the query file for the next round. The current set of matched reads are aligned to the long-contigs, and matching reads are removed (not used for assembly in future rounds). Defaults to "10000".

**-e** *minimum\_score* : Minimum spliced alignment score to accept as a potential hit for the query.
Defaults to "0.5".

**-c** *minimum\_coverage* : Minimum fraction of coverage of the query by the spliced alignment.
Defaults to "0.8".

#### SRAssembler behavior

**-n** *num_rounds* : Maximum number of rounds of chromosome walking.
Defaults to "10".

**-E** *extra_rounds* : Number of additional rounds of recursion to perform once a hit contig is found.
This option may be useful when automating many SRAssembler runs in which extending the contig beyond the bounds of the query (e.g., into the UTRs surrounding a CDS) is desirable, but the number of rounds necessary is unpredictable. This does not always lead to longer hit contigs, and may in fact result in a worse assembly. This option will not cause SRAssembler to exceed the number of rounds set by **-n**. Defaults to "0".

**-a** *assemble_round* : The round in which to start read assembly, defaults to "1".
If read coverage depth is good, starting in round 1 is fine. In some cases, however, using the reads found in round 1 as queries to find additional reads before doing contig assembly can lead to better results.

**-b** *clean_frequency* : The frequency with which to remove unrelated contigs and reads from the recursive search.
SRAssembler will periodically remove assembled contigs that do not align to the query sequence, and remove reads that are not being assembled into contigs. More frequent cleaning can make SRAssembler run faster, but may aggressively remove reads that would be relevant to correct assembly of a homolog of the query. The default of "4" means that SRAssembler will purge contigs and reads after three rounds of not doing so. This is not the same as doing so on every round evenly divisible by four, because the count from the last cleaning round can be reset by a forced cleaning when too many contigs are assembled, as explained under **-d**, or when continuing a run from pre-existing results. 

**-d** *contig_limit* : The minimum number of assembled contigs to automatically trigger removal of unrelated contigs and reads.
This can prevent slowdown or crashing in cases where an excessive number of reads are found due to things like common domains or transposable elements. If set to "0", do not remove unrelated contigs and reads except as scheduled by **-b**. Defaults to "500".

**-j** *taboo_file* : FASTA file containing DNA or protein sequences used to taboo reads and prevent them from being used for assembly.
SRAssembler may struggle with reads from genomes with lots of transposons or other repetitive elements, or with assembling homologs to a query that contains a very common domain. With this option, SRAssembler will start by searching for reads that match the *taboo_file*, and exclude them from all assemblies.

**-J** Taboo file type.
Accepts "protein" or "dna" (default: "protein").

**-z** *masking_round* : The round in which to start masking low-complexity regions of intermediate contigs before searching for reads.
Reads with low-complexity regions that are assembled into contigs can wreak havoc on SRAssembler function. The number of irrelevant reads found in the next round due to shared low-complexity sequences can dwarf the number of relevant reads, and can slow down or crash the assembler. By default this option is set to "1", and SRAssembler uses NCBI's [DustMasker](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/) to mask low-complexity regions in the intermediate contigs in every round. If you want to disable this feature, set this option to "0".

**-x** *tip\_search\_round* : The round in which to start masking the center of query contigs.
If you find SRAssembler runs slowly due to finding excessive numbers of new reads, it can be helpful to mask the center of the intermediate contigs, leaving only the contig 'tips' to match against new reads. This causes SRAssembler to prioritize finding reads that can extend the contigs, at the cost of possibly missing some reads that might contribute to a hit contig. When set to "0", no masking of the center of contigs will occur. Defaults to "2", a value that ensures that contigs derived from the inital round of searching using the query sequence to find reads (rather than intermediate contigs) are used fully before masking. If using center masking, actively experiment with this option, as it can have a significant effect on whether optimal hit contigs are found.

**-X** *tip\_search\_length* : Length of intermediate contigs to leave unmasked at both ends as queries to find new reads.
When using this option, it should at minimum be as long as the match length requirements set in the config file input with **-p**. When set to the default of "0", no masking of the center of contigs will occur.

**-f** : Do not check for contigs that meet the completion thresholds.
Force SRAssembler to continue chromosome walking for the number of rounds specified with **-n**. This is useful if you want to build contigs longer than just the coding region in order to capture promoter regions or nearby elements.

**-y** : Do not resume chromosome walking from existing SRAssembler results in the output directory specified by **-o**.
By default, SRAssembler will continue the previous run if the **-n** option is larger than the final round in the output directory. This behavior is useful if you want to extend the hit contigs for only a few rounds or change your completion thresholds to be more selective. Using the **-y** option will cause SRAssembler to overwrite the results in the output directory and start from round 1.

#### Help

**-v** : Verbose mode.
Not only is output to stdout more verbose, but several files in the aux/ directory that are normally removed to keep things tidy are left behind. This mode is typically only used for debugging.

**-h** : Print an SRAssembler usage synopsis and exit.

## Output

If no output directory was specified with **-o**, SRAssembler will create a directory called "SRAssembler\_output" within the working directory. 
If no processed reads directory was specified with **-r**, the processed reads are created in a directory called "processed\_reads" within the output directory.
If no temp directory was specified with **-t**, SRAssembler will use "/tmp".
Assuming these defaults, SRAssembler will create the following directories and files:

- **SRAssembler_output/summary.html** : A summary of the state of the run at each round, and the final results.
- **SRAssembler\_output/all\_contigs.fasta** : All of the contigs longer than the *query\_contig\_minimum* (set with **-i**) that were assembled in the final round.
- **SRAssembler\_output/hit\_contigs.fasta** : The contigs assembled in the final round that meet the completion thresholds (set with **-m**, **-M**, **-e**, and **-c**).
- **SRAssembler\_output/output.aln** : The report on spliced alignment of the query to the contigs in *SRAssembler\_output/all\_contigs.fasta*.
- **SRAssembler\_output/output.ano** : The *ab initio* gene prediction report, if a gene-finder was used.
- **SRAssembler\_output/msg.log** : The detailed log file of the SRAssembler run. 
Contains, among other things, all of the command-line options and the contents of the parameter file used for this run.
- **SRAssembler\_output/intermediates/** : This directory contains FASTA files of the assembled contigs from each round.
It also contains files that end ".masked", showing the effects, if any, of using the contig masking options, as well as files that end ".beforeclean" that correspond to the original assembled contigs in "cleaning" rounds (see option **-b**) in which unrelated contigs were purged.
- **SRAssembler\_output/processed\_reads/** : This directory contains subdirectories, one for each library specified in the run.
Within these subdirectories are FASTA files containing subsets of the reads from the library (number of reads per file set by **-R**) and the aligner indices of those read files. This directory will also contain symbolic links named "lib1", "lib2", etcetera, linking to each of the libraries specified for this run. Because of this, different read libraries can be used in different SRAssembler runs without the need to preprocess the reads multiple times.
- **SRAssembler\_output/aux/** : This directory stores temporary files produced during the SRAssembler run.
These can serve as checkpoint files so that SRAssembler can continue previous assembling. This directory will have more files in it if SRAssembler was run in verbose mode (with **-v**), which may be helpful for debugging.
- **/tmp/SRAssemblertemp*0123456789*/** : Rather than having a specific name, this directory will take the numeric string of its name from the unique process ID of the SRAssembler run.
This directory is entirely for the temporary storage of many very small files that are produced during an SRAssembler run. At the beginning of a run, a symbolic link to this directory is also created in SRAssembler\_output/. Upon completion of an SRAssembler run, both this directory and the symbolic link will be deleted, except during a verbose (**-v**) run, in which case the symbolic link will be deleted and this directory will be moved into SRAssembler\_output/. If SRAssembler is interrupted, this directory and the symbolic link will be left behind. The symbolic link is to remind the user to manually remove this directory in these cases.

## Testing installation

The [demo/](demo/) directory contains the necessary files to test the functionality of most of SRAssembler's options. The bash script [xtest](demo/xtest) runs SRAssembler several times using various command-line options and compares the results to the summary files in [demo/standards/](demo/standards). This confirms the installation and functionality of all the SRAssembler dependencies. The bash script xquicktest runs just tests 1 and 4, confirming the basic functionality of the SRAssembler and SRAssembler_MPI executables. The demo/ directory also contains two library configuration files and a sample parameter configuration file. Examining these files is a good way to familiarize yourself with setting up SRAssembler runs.

Within [demo/input/](demo/input/) are two FASTA query files and four FASTQ reads files. The query files are the transcript sequence of At1g01950.1, an *Arabidopsis thaliana* gene located at Chr1:325316-330619, and the amino acid sequence of the homologous rice protein LOC_Os06g04560.1. We used SAMTools' wgsim program to simulate two libraries from region 300,000..400,000 of chromosome 1 of Arabidopsis, with insert sizes 200bp and 1kb. Each library has 50,000 paired reads. The base error rate is 0.002 and the fraction of indels is 0.001.

Running the tests is as simple as executing `bash xtest`. Each test will output progress to standard out, and after each test the script will report whether the test's summary file matches the one in the standards/ directory. At the end of the script, the comparison of all the tests' results to the standards will again be reported for easy viewing. To test the functionality of the Singularity image, run `singularity exec -e -B $(pwd) path/to/SRAssembler.simg ./xtest defaultpath` from within the demo/ directory. This will bind the demo directory to the Singularity instance, allowing it to read and write, but cause it to use only its internal SRAssembler software to execute the tests.

## License

The SRAssembler package conforms to our [RAMOSE](https://github.com/BrendelGroup/) philosophy and promise: we made every attempt to make the software such that you can produce **reproducible**, **accurate**, and **meaningful** results; our software is **open** (source), designed to be **scalable**, and **easy** to use (so that a typical genome laboratory should be able to run it).  And, of course, our software should be _ramose_ ("having many branches") - it's yours to use, modify, improve, and enhance under the GNU General Public License v3.0.

## Reference

Thomas W. McCarthy, Hsien-Chao Chou, and Volker P. Brendel (2018) **SRAssembler:**
_Selective Recursive local Assembly of homologous genomic regions._ To be submitted.

## Contact

Please direct all comments and suggestions to
[Volker Brendel](<mailto:vbrendel@indiana.edu>)
at [Indiana University](http://brendelgroup.org/).

