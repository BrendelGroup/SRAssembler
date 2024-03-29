/*
 * SRAssembler.cpp
 *
 *  Created on: Oct 12, 2011
 *      Authors: Hsien-chao Chou (first version); Thomas McCarthy and Volker Brendel (modifications)
 */

#include <ctime>
#include "SRAssembler.h"
#include "SRAssemblerMaster.h"
#include "SRAssemblerSlave.h"
#include "SRAssemblerSlave.h"
#include <vector>
#include <regex>


SRAssembler* SRAssembler::_srassembler = NULL;

SRAssembler::SRAssembler() {
	// Auto-generated constructor stub
}

SRAssembler::~SRAssembler() {
	// Auto-generated destructor stub
}

int SRAssembler::init(int argc, char * argv[], int rank, int mpiSize) {
	// Set the default values.
	// TODO output these to the logger.
	init_match_length = INIT_MATCH_LENGTH_PROTEIN;
	recur_match_length = RECUR_MATCH_LENGTH;
	mismatch_allowed = MISMATCH_ALLOWED_PROTEIN;
	num_rounds = NUM_ROUNDS;
	verbose = VERBOSE;
	assembly_round = ASSEMBLY_ROUND;
	ignore_contig_explosion = false;
	reads_per_file = READS_PER_FILE;
	fastq_format = FASTQ_INTERLEAVED;
	out_dir = OUT_DIR;
	data_dir = "";
	tmp_loc = TMP_LOC;
	probe_type = QUERY_TYPE;
	taboo_file = "";
	taboo_type = QUERY_TYPE;
	assembler_program = ASSEMBLER_PROGRAM;
	gene_finding_program = GENE_FINDING_PROGRAM;
	spliced_alignment_program = SPLICED_ALIGNMENT_PROGRAM;
	species = DEFAULT_SPECIES;
	start_k = START_K;
	end_k = END_K;
	step_k = STEP_K;
	clean_round = CLEAN_ROUND;
	contig_limit = CONTIG_LIMIT;
	over_write = OVER_WRITE;
	check_gene_assembled = CHECK_GENE_ASSEMBLED;
	preprocessing_only = PREPROCESSING_ONLY;
	min_score = MIN_SCORE;
	min_coverage = MIN_COVERAGE;
	query_contig_min = QUERY_CONTIG_MIN;
	min_contig_lgth = MIN_CONTIG_LGTH;
	max_contig_lgth = MAX_CONTIG_LGTH;
	masking_round = MASKING_ROUND;
	end_search_length = END_SEARCH_LENGTH;
	end_search_round = END_SEARCH_ROUND;
	extra_rounds = EXTRA_ROUNDS;
	bool k_format = true;
	this->rank = rank;
	this->mpiSize = mpiSize;
	unsigned int insert_size = 0;
	string left_read = "";
	string right_read = "";
	// Parse command line options.
	bool default_e = true;
	bool default_i = true;
	preprocessed_exist = false;
	usage = "SRAssembler v" + VERSION + " usage:\n";
	usage.append("    SRAssembler [options] -q query_file -p parameter_file -l library_file\n");
	usage.append("\n");
	usage.append("-q: Required unless pre-processing only; FASTA-formatted query file.\n");
	usage.append("-t: Query file type; options: 'protein', 'dna'. [Default: " + QUERY_TYPE + "]\n");
	usage.append("-p: Required unless pre-processing only; SRAssembler parameter config file.\n");
	usage.append("-o: SRAssembler output directory. [Default: " + OUT_DIR + "]\n");
	usage.append("-l: Required if the -1 option is not used; sequencing reads library file.\n");
	usage.append("-r: Directory in which to store or from which to retrieve processed reads.\n");
	usage.append("    [Default: output_directory/" + READS_DATA + "]\n");
	usage.append("-k: Specifies the k-mer set to be used by the assembler. [Default: " + int2str(START_K) + ":" + int2str(STEP_K) + ":" + int2str(END_K) + "]\n");
	usage.append("-n: Maximum number of rounds for chromosome walking. [Default: " + int2str(NUM_ROUNDS) + "]\n");
	usage.append("-f: Forgo spliced alignment check after intermediate assembly rounds.\n");
	usage.append("    SRAssembler is forced to continue for the -n specified number of rounds.\n");
	usage.append("-y: Disable SRAssembler resumption from previous checkpoint (will overwrite\n");
	usage.append("    existing output).\n");
	usage.append("\n");
	usage.append("-h: Print extended usage synopsis.");

	extended_usage = "SRAssembler v" + VERSION + " usage:\n";
	extended_usage.append("    SRAssembler [options] -q query_file -p parameter_file (-l library_file | \n");
	extended_usage.append("                -1 reads_file1 [-2 reads_file2 -Z insertsize])\n");
	extended_usage.append("\n");
	extended_usage.append("-q: Required unless pre-processing only; FASTA-formatted query file.\n");
	extended_usage.append("-t: Query file type; options: 'protein', 'dna'. [Default: " + QUERY_TYPE + "]\n");
	extended_usage.append("-p: Required unless pre-processing only; SRAssembler parameter config file.\n");
	extended_usage.append("-o: SRAssembler output directory. [Default: " + OUT_DIR + "]\n");
	extended_usage.append("-T: Directory for temporary file storage. [Default: " + TMP_LOC + "]\n");
	extended_usage.append("\n");
	extended_usage.append("-l: Required if the -1 option is not used; sequencing reads library file.\n");
	extended_usage.append("-1: Required if the -l option is not used; use this option to specify the\n");
	extended_usage.append("    single-end reads file or the left-end reads file for paired-end reads.\n");
	extended_usage.append("-2: Right-end reads file for paired-end reads.\n");
	extended_usage.append("-Z: Insert size of the paired-end reads. Required if the -2 option is used.\n");
	extended_usage.append("-R: Number of reads per pre-preprocessed reads file. [Default: " + int2str(READS_PER_FILE) + "]\n");
	extended_usage.append("-r: Directory in which to store or from which to retrieve processed reads.\n");
	extended_usage.append("    [Default: output_directory/" + READS_DATA + "]\n");
	extended_usage.append("-P: Run the read pre-processing step only, then terminate SRAssembler.\n");
	extended_usage.append("\n");
	extended_usage.append("-A: Contig assembler program; options: 0=>SOAPdenovo2, 1=>ABySS. [Default: " + int2str(ASSEMBLER_PROGRAM) + "]\n");
	extended_usage.append("-k: Specifies the k-mer set to be used by the assembler. [Default: " + int2str(START_K) + ":" + int2str(STEP_K) + ":" + int2str(END_K) + "]\n");
	extended_usage.append("    Format: start_k:interval:end_k. Start_k and end_k must be odd integers, and\n");
	extended_usage.append("    interval must be an even integer, as in the following example: '15:10:45'\n");
	extended_usage.append("    specifies that k-mer values 15, 25, 35, 45 will be tested.\n");
	extended_usage.append("-S: Spliced alignment program; 0=>GeneSeqer, 1=>GenomeThreader. [Default: " + int2str(SPLICED_ALIGNMENT_PROGRAM) + "]\n");
	extended_usage.append("-s: Species model for spliced alignment; options: 'human', 'mouse', 'rat', \n");
	extended_usage.append("    'chicken', 'drosophila', 'nematode', 'aspergillus', 'fission_yeast', \n");
	extended_usage.append("    'arabidopsis', 'maize', 'rice', 'medicago'. [Default: " + DEFAULT_SPECIES + "]\n");
	extended_usage.append("-G: Ab initio gene finding program; options: 0=>None, 1=>Snap. [Default: " + int2str(GENE_FINDING_PROGRAM) + "]\n");
	extended_usage.append("\n");
	extended_usage.append("-i: Minimum contig length for chromosome walking. [Default: " + int2str(QUERY_CONTIG_MIN) + "]\n");
	extended_usage.append("-m: Minimum contig length to accept as a hit. [Default: " + int2str(MIN_CONTIG_LGTH) + "]\n");
	extended_usage.append("-M: Maximum contig length to accept as a hit. [Default: " + int2str(MAX_CONTIG_LGTH) + "]\n");
	extended_usage.append("-e: Minimum spliced alignment score for hits. [Default: " + double2str(MIN_SCORE) + "]\n");
	extended_usage.append("-c: Minimum spliced alignment coverage score for hits. [Default: " + double2str(MIN_COVERAGE) + "]\n");
	extended_usage.append("\n");
	extended_usage.append("-n: Maximum number of rounds for chromosome walking. [Default: " + int2str(NUM_ROUNDS) + "]\n");
	extended_usage.append("-E: Extra rounds of recursion performed once a hit contig is found. [Default: " + int2str(EXTRA_ROUNDS) + "]\n");
	extended_usage.append("-a: The round in which to start read assembly. [Default: " + int2str(ASSEMBLY_ROUND) + "]\n");
	extended_usage.append("-b: The frequency with which to periodically remove unrelated contigs and reads.\n");
	extended_usage.append("    For example, '-b 3' specifies that SRAssembler will purge after two rounds\n");
	extended_usage.append("    of not doing so. [Default: " + int2str(CLEAN_ROUND) + "]\n");
	extended_usage.append("-d: The minimum number of assembled contigs to automatically trigger removal of\n");
	extended_usage.append("    unrelated contigs and reads. If set to '0', do not remove unrelated contigs\n");
	extended_usage.append("    and reads except as scheduled by '-b' option. [Default: " + int2str(CONTIG_LIMIT) + "]\n");
	extended_usage.append("-j: FASTA-formatted file containing sequences used to taboo reads and prevent\n");
	extended_usage.append("    them from being used for assembly.\n");
	extended_usage.append("-J: Taboo file type; options: 'protein', 'dna'. [Default: " + QUERY_TYPE + "]\n");
	extended_usage.append("\n");
	extended_usage.append("-f: Forgo spliced alignment check after intermediate assembly rounds.\n");
	extended_usage.append("    SRAssembler is forced to continue for the -n specified number of rounds.\n");
	extended_usage.append("-y: Disable SRAssembler resumption from previous checkpoint (will overwrite\n");
	extended_usage.append("    existing output).\n");
	extended_usage.append("-x: The round in which to start masking the center of query contigs when\n");
	extended_usage.append("    searching for new reads. If set to '0', do not mask center of query contigs.\n");
	extended_usage.append("    [Default: " + int2str(END_SEARCH_ROUND) + "]\n");
	extended_usage.append("-X: Length of contig ends to leave unmasked as queries to find new reads.\n");
	extended_usage.append("    If set to '0', do not mask center of query contigs. [Default: " + int2str(END_SEARCH_LENGTH) + "]\n");
	extended_usage.append("-z: The round in which to start masking low-complexity regions of contigs before\n");
	extended_usage.append("    searching for reads. If set to '0', do not mask low-complexity regions of\n");
	extended_usage.append("    query contigs. [Default: " + int2str(MASKING_ROUND) + "]\n");
	extended_usage.append("\n");
	extended_usage.append("-v: Verbose output (typically only used for debugging).\n");
	extended_usage.append("-h: Print this extended usage synopsis.");


	char c;
	while((c = getopt(argc, argv, "1:2:a:A:b:c:d:e:E:fG:hi:j:J:k:l:m:M:n:o:p:Pq:r:R:s:S:t:T:vx:X:yz:Z:")) != -1) {
		switch (c){
			case '1':
				left_read = optarg;
				break;
			case '2':
				right_read = optarg;
				break;
			case 'a':
				assembly_round = str2int(optarg);
				if (assembly_round == 0)
					assembly_round = 9999;
				break;
			case 'A':
				assembler_program = str2int(optarg);
				break;
			case 'b':
				clean_round = str2int(optarg);
				break;
			case 'c':
				min_coverage = str2double(optarg);
				break;
			case 'd':
				contig_limit = str2int(optarg);
				if (contig_limit == 0) {
					ignore_contig_explosion = true;
				}
				break;
			case 'e':
				min_score = str2double(optarg);
				break;
			case 'E':
				extra_rounds = str2int(optarg);
				break;
			case 'f':
				check_gene_assembled = false;
				break;
			case 'G':
				gene_finding_program = str2int(optarg);
				break;
			case 'h':
				show_extended_usage();
				return -1;
				break;
			case 'i':
				query_contig_min = str2int(optarg);
				break;
			case 'j':
				taboo_file = optarg;
				break;
			case 'J':
				taboo_type = optarg;
				break;
			case 'k': {
				vector<string> tokens;
				tokenize(optarg, tokens, ":");
				if (tokens.size() != 3)
					k_format = false;
				start_k = str2int(tokens[0]);
				step_k = str2int(tokens[1]);
				end_k = str2int(tokens[2]);
				if (start_k <= 0 || start_k > 150 || start_k % 2 == 0) {
					k_format = false;
					break;
				}
				if (end_k <= 0 || end_k > 150 || end_k % 2 == 0 || end_k < start_k) {
					k_format = false;
					break;
				}
				if (step_k <= 0 || step_k > 150 || step_k % 2 == 1) {
					k_format = false;
					break;
				}
				int k_value = start_k;
				while (k_value < end_k)
					k_value += step_k;
				if (k_value > step_k)
					end_k = k_value;
				break;
			}
			case 'l':
				library_file = optarg;
				break;
			case 'm':
				min_contig_lgth = str2int(optarg);
				break;
			case 'M':
				max_contig_lgth = str2int(optarg);
				break;
			case 'n':
				num_rounds = str2int(optarg);
				break;
			case 'o':
				out_dir = optarg;
				break;
			case 'p':
				param_file = optarg;
				break;
			case 'P':
			//TODO turn this into an optional argument for setting the vmatch prefixlength before preprocessing.
				preprocessing_only = true;
				break;
			case 'q':
				probe_file = optarg;
				break;
			case 'r':
				data_dir = optarg;
				break;
			case 'R':
				reads_per_file = str2int(optarg);
				break;
			case 's':
				species = optarg;
				break;
			case 'S':
				spliced_alignment_program = str2int(optarg);
				break;
			case 't':
				probe_type = optarg;
				break;
			case 'T':
				tmp_loc = optarg;
				break;
			case 'v':
				verbose = true;
				break;
			case 'x':
				end_search_round = str2int(optarg);
				break;
			case 'X':
				end_search_length = str2int(optarg);
				break;
			case 'y':
				over_write = true;
				break;
			case 'z':
				masking_round = str2int(optarg);
				break;
			case 'Z':
				insert_size = str2int(optarg);
				break;
			case '?':
				string msg = "input error! unknown option : -";
				msg.append(1,optopt);
				logger->error(msg);
				return -1;
		}
	}
	if (argc == 1){
		show_usage();
		return -1;
	}
	// Establish the output directories.
	aux_dir = out_dir + "/aux";
	if (data_dir == "")
		data_dir = out_dir + "/" + READS_DATA;
	preprocessed_exist = file_exists(data_dir);
	intermediate_dir = out_dir + "/intermediates";
	log_file = out_dir + "/msg.log";
	spliced_alignment_output_file = out_dir + "/output.aln";
	gene_finding_output_file = out_dir + "/output.ano";
	gene_finding_output_protein_file = out_dir + "/snap.predicted.prot";
	final_contigs_file = out_dir + "/all_contigs.fasta";
	summary_file = out_dir + "/summary.html";
	hit_contigs_file = out_dir + "/hit_contigs.fasta";

	int level = Logger::LEVEL_INFO;
	if (verbose)
		level = Logger::LEVEL_DEBUG;
	logger = Logger::getInstance(level, log_file);

	// Check for missing or illegal options.

	// Create/confirm the existence of the output directory.
	run_shell_command("mkdir -p " + out_dir);
	if (!file_exists(out_dir)){
		logger->error("output directory: " + out_dir + " cannot be created!");
		return -1;
	}

	if (reads_per_file <= 1000){
		logger->error("-R value is not valid or too small");
		return -1;
	}

	if (library_file != ""){
		if (!file_exists(library_file)){
			logger->error("library file: " + library_file + " does not exist!");
			return -1;
		}
		if (!read_library_file())
			return -1;
	// If not using a library configuration file, library attributes are identified from the command line.
	} else {
		if (left_read == "" && !preprocessed_exist){
			logger->error("-1 or library file is required");
			return -1;
		}
		if (!file_exists(left_read) && !preprocessed_exist){
				logger->error("file: " + left_read + " does not exist!");
				return -1;
		}
		Library lib(0, this->data_dir, this->aux_dir, this->logger);
		lib.set_format(FORMAT_FASTQ);
		lib.set_left_read(left_read);
		lib.set_library_name(get_file_base_name(lib.get_left_read()) + "library");
		if (right_read != ""){
			if (insert_size < 1) {
				logger->error("Paired end reads require an insert size");
				return -1;
			}
			lib.set_paired_end(true);
			lib.set_right_read(right_read);
			lib.set_insert_size(insert_size);
			if (!file_exists(right_read) && !preprocessed_exist){
				logger->error("file: " + right_read + " does not exist!");
				return -1;
			}
		}
		this->libraries.push_back(lib);
	}

	// If we're only preprocessing, we don't need the rest of the checks.
	if (preprocessing_only){
		return 0;
	}

	// Confirm the existence of the probe file.
	if (!file_exists(probe_file)){
		logger->error("query file: " + probe_file + " does not exist!");
		return -1;
	}

	if (taboo_file != ""){
		if (!file_exists(taboo_file)){
			logger->error("Taboo file: " + library_file + " does not exist!");
			return -1;
		}
	}

	if (param_file != ""){
		if (!file_exists(param_file)){
			logger->error("Parameter file : " + param_file + " does not exist!");
			return -1;
		} else {
		parameters_dict = read_param_file();
		}
	} else {
		logger->error("Parameter file required!");
		return -1;
	}

	if (num_rounds <= 0){
		logger->error("-n must be larger than 0");
		return -1;
	}

	if (probe_type != TYPE_PROTEIN && probe_type != TYPE_DNA){
		logger->error("-t must be 'protein' or 'dna'");
		return -1;
	} else {
		if (default_e && probe_type == TYPE_DNA)
			mismatch_allowed = MISMATCH_ALLOWED_DNA;
		if (default_i && probe_type == TYPE_DNA)
			init_match_length = INIT_MATCH_LENGTH_DNA;
	}
	if (species != "human" && species != "mouse" && species != "rat" && species != "chicken" && species != "drosophila" && species != "nematode" && species != "fission_yeast" && species != "aspergillus" && species != "arabidopsis" && species != "maize" && species != "rice" && species != "medicago"){
		logger->error("-s valid species is required");
		return -1;
	}
	if (!k_format){
		logger->error("-k format error. The format is : start_k:interval:end_k. The start_k and end_k must be odd values, and the interval must be an even value.");
		return -1;
	}

	return 0;
}

boost::unordered_map<std::string,Params> SRAssembler::read_param_file() {
	ifstream param_file(this->param_file.c_str());
	string line;
	boost::unordered_map<std::string,Params> parameters_dict;
	Params params;
	string program_name;
	bool found_program = false;
	while (getline(param_file, line)){
		line = trim(line);
		if (line.length() == 0) continue;
		if (line.substr(0,1) == "#") continue;
		if (line.substr(0,1) == "[") {
			if (found_program) {
				parameters_dict.insert(make_pair(program_name, params));
				params.clear();
			}
			// This should capture the string between the brackets
			program_name = line.substr(1, line.length() - 2);
			found_program = true;
			continue;
		}
		if (found_program){
			vector<string> tokens;
			tokenize(line, tokens, "=");
			if (tokens.size() == 2){
				string param = trim(tokens[0]);
				string value = trim(tokens[1]);
				params.insert(Params::value_type(param, value));
			}
			// On/off flag parameters are allowed.
			else if (tokens.size() == 1){
				string param = trim(tokens[0]);
				string value = "";
				params.insert(Params::value_type(param, value));
			}
		}
	}
	// Capture the last set of parameters.
	if (found_program) {
		parameters_dict.insert(make_pair(program_name, params));
		params.clear();
	}
	return parameters_dict;
}

Params SRAssembler::get_parameters(string program_name) {
	return this->parameters_dict[program_name];
}

bool SRAssembler::read_library_file() {
	ifstream lib_file(this->library_file.c_str());
	string line;
	Library* lib = NULL;
	while (getline(lib_file, line)){
		line = trim(line);
		if (line == "[LIBRARY]") {
			if (lib != NULL){
				if (lib->get_left_read() == ""){
					logger->error("r1 file is expected in library config file!");
					return false;
				}
				// If the user didn't name their library, name it after the left read file.
				if (lib->get_library_name() == "") {
					lib->set_library_name(get_file_base_name(lib->get_left_read()) + "library");
				}
				this->libraries.push_back(*lib);
			}
			lib = new Library(this->libraries.size(), this->data_dir, this->aux_dir, this->logger);
		} else {
			vector<string> tokens;
			tokenize(line, tokens, "=");
			if (tokens.size() == 2){
				string param = trim(tokens[0]);
				string value = trim(tokens[1]);
				if (param == "library_name" && lib != NULL) {
					lib->set_library_name(value);
				}
				if (param == "insert_size" && lib != NULL) {
					lib->set_insert_size(str2int(value));
				}
				if (param == "direction" && lib != NULL) {
					lib->set_reversed(str2int(value));
				}
				if (param == "r1" && lib != NULL) {
					lib->set_left_read(value);
					if (!file_exists(lib->get_left_read()) && !preprocessed_exist) {
						logger->error("r1 file in config file: " + lib->get_left_read() + " does not exist!");
						return false;
					}
				}
				if (param == "r2" && lib != NULL) {
					lib->set_right_read(value);
					if (!file_exists(lib->get_right_read()) && !preprocessed_exist) {
						logger->error("r2 file in config file: " + lib->get_right_read() + " does not exist!");
						return false;
					}
					lib->set_paired_end(true);
				}
				if (param == "format" && lib != NULL) {
					if (value == "fastq") {
						lib->set_format(FORMAT_FASTQ);
					}
					else if (value == "fasta") {
						lib->set_format(FORMAT_FASTA);
					} else {
						logger->error("format in library config file should be 'fastq' or 'fasta'!");
						return false;
					}
				}
			}
		}
	}
	if (lib != NULL){
		if (lib->get_left_read() == ""){
			logger->error("r1 file is expected in library config file!");
			return false;
		}
		// If the user didn't name their library, name it after the left read file.
		if (lib->get_library_name() == "") {
			lib->set_library_name(get_file_base_name(lib->get_left_read()) + "library");
		}
		this->libraries.push_back(*lib);
	}
	if (this->libraries.size() == 0){
		logger->error(" No [LIBRARY] section found in library config file!");
		return false;
	}
	return true;

}

int SRAssembler::get_file_count(string search_pattern){
	string cmd = "ls -U1dq " + search_pattern + " | wc -l";
	return str2int(run_shell_command_with_return(cmd));
}

void SRAssembler::preprocess_read_chunk(int lib_idx, int read_chunk){
	Library lib = this->libraries[lib_idx];
	logger->running("preprocessing lib " + int2str(lib_idx + 1) + ", reads file (" + int2str(read_chunk) + "/" + int2str(lib.get_num_chunks()) + ")");
	string suffix = int2str(read_chunk) + ".fasta";
	string left_read_file = lib.get_split_read_prefix(lib.get_left_read()) + suffix;
	string right_read_file;
	if (lib.get_paired_end())
		right_read_file = lib.get_split_read_prefix(lib.get_right_read()) + suffix;
	// Create the Vmatch mkvtree indexes.
	Aligner* aligner = get_aligner(0);  // Round 0 means DNA Aligner
	aligner->create_index(lib.get_read_chunk_index_name(read_chunk, LEFT_READ), "dna", left_read_file);
	if (lib.get_paired_end())
		aligner->create_index(lib.get_read_chunk_index_name(read_chunk, RIGHT_READ), "dna", right_read_file);
}

void SRAssembler::send_code(const int& to, const int& action, const int& value1, const int& value2, const int& value3){
	mpi_code code;
	code.action = action;
	code.value1 = value1;
	code.value2 = value2;
	code.value3 = value3;
	mpi_send(get_mpi_code_value(code), to);
}
void SRAssembler::broadcast_code(const int& action, const int& value1, const int& value2, const int& value3){
	mpi_code code;
	code.action = action;
	code.value1 = value1;
	code.value2 = value2;
	code.value3 = value3;
	mpi_bcast(get_mpi_code_value(code));
}

string SRAssembler:: get_contigs_index_name(int round){
   return aux_dir + "/round" + int2str(round);
}


string SRAssembler:: get_query_fasta_file_name(int round){
/* It appears that this function is deceptively named.
 * It does not just return a string naming a fasta file that contains all of the matched reads.
 * This function is also responsible for assembling the contents of that file.
 */
 	if (round > 1){
		// If we have passed the round to start assembling, use the assembled contig files as the query.
		if (assembly_round < round)
			return get_contig_file_name(round-1);
		// If we are still waiting to assemble, use the collected reads as queries.
		else {
			string joined_file = aux_dir + "/matched_reads_joined.fasta";
			return joined_file;
		}
	}
	return probe_file;
}

void SRAssembler::mask_contigs(int round){
	string cmd;
	string contig_file;
	contig_file = get_contig_file_name(round);
	string masked_file = contig_file + ".masked";
	cmd = "cat " + contig_file;
	if (masking_round > 0 && round >= masking_round) {
		// NCBI's dustmasker identifies regions of low complexity.
		cmd += "| dustmasker -outfmt fasta ";
		// Sed command replaces lowercase atcg NOT in the header with capital Ns.
		cmd += "| sed --line-length=0 -e '/^[^>]/s/[atcg]/N/g' ";
		// AWK command converts FASTA to single-line.
		cmd += "| awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }'";
	}
	if (end_search_length > 0 && end_search_round > 0 && round >= end_search_round) {
		// A more complex AWK command will mask any base farther than end_search_length away from the ends of the contigs.
		cmd += "| awk -v keeplength=" + int2str(end_search_length);
		// Set the field separator to an empty string so that AWK can act on every character.
		cmd += " 'BEGIN{FS=OFS=\"\"} ";
		// If the line isn't a header line, change the middle characters into 'n's, then print the current (modified) line.
		cmd += "!/^>/{for (i = keeplength + 1; i <= NF - keeplength; i++) $i=\"n\"} {print $0}' ";
	}
	cmd += "> " + masked_file;
	logger->debug(cmd);
	logger->fragile_run_shell_command(cmd);
}

string SRAssembler:: get_query_fasta_file_name_masked(int round){
	if (round == 0) {
		return taboo_file;
	} else if (round > 1){
		// If we have passed the round to start assembling, use the assembled contig files as the query.
		if (assembly_round < round) {
			// Depending on masking options, this file may not actually be masked in any way.
			string masked_fasta = get_contig_file_name(round-1) + ".masked";
			return masked_fasta;
		// If we are still waiting to assemble, use the collected reads as queries.
		} else {
			string joined_file = aux_dir + "/matched_reads_joined.fasta";
			return joined_file;
		}
	}
	return probe_file;
}

string SRAssembler:: get_contig_file_name(int round){
	return intermediate_dir + "/contigs_" + "r" + int2str(round) + ".fasta";
}

string SRAssembler:: get_matched_reads_file_name(int round){
	return aux_dir + "/found_reads_" + "r" + int2str(round) + ".list";
}

int SRAssembler::do_alignment(int round, int lib_idx, int read_chunk) {
	Library lib = this->libraries[lib_idx];
	Aligner* aligner = get_aligner(round);
	string program_name = aligner->get_program_name();
	string criteria;
	string program;
	if (round == 0) {
		criteria = taboo_type + "_taboo";
	} else if (round == 1) {
		criteria = get_type(1) + "_init";
	} else {
		criteria = "extend_contig";
	}
	program = program_name + "_" + criteria;
	logger->running("Aligning using " + program_name + " criteria " + criteria +": Round " + int2str(round) + ", Lib " + int2str(lib_idx+1) + " of " + int2str(this->libraries.size()) + ", Reads chunk " + int2str(read_chunk) + " of " + int2str(lib.get_num_chunks()));
	Params params = this->get_parameters(program);
	int new_read_count;

	// Reads as queries are necessary when searching against a protein.
	if (round == 1 && probe_type == "protein") {
		aligner->do_alignment(aux_dir + "/qindex", "protein", get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_chunk, LEFT_READ), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		// If library is paired end, the found reads from the second half will be appended to the vmatch_output_file.
		if (lib.get_paired_end())
			aligner->do_alignment(aux_dir + "/qindex", "protein", get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_chunk, RIGHT_READ), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		new_read_count = aligner->parse_output(get_vmatch_output_filename(round, lib_idx, read_chunk), found_reads, lib_idx, read_chunk, lib.get_read_chunk_index_name(read_chunk, LEFT_READ), lib.get_read_chunk_index_name(read_chunk, RIGHT_READ), lib.get_matched_left_reads_filename(round, read_chunk), lib.get_matched_right_reads_filename(round, read_chunk));
		return new_read_count;
	} else if (round == 0 && taboo_type == "protein") {
		aligner->do_alignment(aux_dir + "/tabooindex", "protein", get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_chunk, LEFT_READ), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		if (lib.get_paired_end()) {
			aligner->do_alignment(aux_dir + "/tabooindex", "protein", get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_chunk, RIGHT_READ), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		}
		new_read_count = aligner->ignore_output(get_vmatch_output_filename(round, lib_idx, read_chunk), found_reads, lib_idx, read_chunk);
	/* If we have a dna probe, in round 1 get_query_fasta_file_name_masked() will return the probe_file.
	 * If round is 0, the taboo_file will be used to find reads to exclude.
	 * After round 1 we use the masked contig file as the query.
	 * If masking is disabled the contigs will not be changed but will still be in the masked contig file.
	 * These alignments are much faster than reads against protein because the reads are indexed.
	 */
	} else {
		aligner->do_alignment(lib.get_read_chunk_index_name(read_chunk, LEFT_READ), "dna", get_match_length(round), get_mismatch_allowed(round), get_query_fasta_file_name_masked(round), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		if (lib.get_paired_end()) {
			aligner->do_alignment(lib.get_read_chunk_index_name(read_chunk, RIGHT_READ), "dna", get_match_length(round), get_mismatch_allowed(round), get_query_fasta_file_name_masked(round), params, get_vmatch_output_filename(round, lib_idx, read_chunk));
		}
		if (round == 0) {
			new_read_count = aligner->ignore_output(get_vmatch_output_filename(round, lib_idx, read_chunk), found_reads, lib_idx, read_chunk);
		} else {
			new_read_count = aligner->parse_output(get_vmatch_output_filename(round, lib_idx, read_chunk), found_reads, lib_idx, read_chunk, lib.get_read_chunk_index_name(read_chunk, LEFT_READ), lib.get_read_chunk_index_name(read_chunk, RIGHT_READ), lib.get_matched_left_reads_filename(round, read_chunk), lib.get_matched_right_reads_filename(round, read_chunk));
		}
	}
	return new_read_count;
}

void SRAssembler::do_assembly(int round, int k, int threads) {
	logger->running("Doing assembly: round = " + int2str(round) + ", k = " + int2str(k));
	Assembler* assembler = get_assembler();
	assembler->do_assembly(k, this->libraries, aux_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round), threads, parameters_dict);
}

void SRAssembler::do_spliced_alignment() {
	logger->running("Doing the final spliced alignment ...");
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	string program_name = spliced_aligner->get_program_name();
	Params params = this->get_parameters(program_name);
	spliced_aligner->do_spliced_alignment(this->final_contigs_file, this->probe_type, this->probe_file, this->species, params, this->spliced_alignment_output_file);
	spliced_aligner->get_hit_contigs(min_score, min_coverage, min_contig_lgth, this->final_contigs_file, this->hit_contigs_file, this->spliced_alignment_output_file, this->best_hits);
	logger->running("Done with final spliced alignment.");
}

string_map SRAssembler::do_spliced_alignment(int round) {
	logger->running("Now running the spliced alignment program on the contigs of round " + int2str(round) + " ...");
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	string program_name = spliced_aligner->get_program_name();
	Params params = this->get_parameters(program_name);
	string contig_file = get_contig_file_name(round);
	string output_file = aux_dir + "/query-vs-contig_" + "r" + int2str(round) + ".aln";
	string hit_file = aux_dir + "/hit_contigs_" + "r" + int2str(round) + ".fasta";
	spliced_aligner->do_spliced_alignment(contig_file, this->probe_type, this->probe_file, this->species, params, output_file);
	string_map query_map = spliced_aligner->get_aligned_contigs(min_score, min_coverage, min_contig_lgth, contig_file, hit_file, output_file, round, best_hits);
	// Auxiliary files are removed here.
	spliced_aligner->clean_files(contig_file);
	logger->running("Done with spliced alignment for round " + int2str(round) + ".");
	return query_map;
}

int SRAssembler::do_spliced_alignment(int round, int k) {
	int best_spliced_length;
	logger->debug("Now running the spliced alignment program for round " + int2str(round) + ", k-mer " + int2str(k) + " ...");
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	string program_name = spliced_aligner->get_program_name();
	Params params = this->get_parameters(program_name);
	string contig_file = this->get_assembly_file_name(round, k);
	string output_file = this->get_spliced_alignment_file_name(round, k);
	spliced_aligner->do_spliced_alignment(contig_file, this->probe_type, this->probe_file, this->species, params, output_file);
	best_spliced_length = spliced_aligner->get_longest_match(round, k, min_score, query_contig_min, output_file, best_hits);
	// Auxiliary files are removed here.
	spliced_aligner->clean_files(contig_file);
	logger->debug("Done with spliced alignment for round " + int2str(round) + ", k-mer " + int2str(k) + ".");
	return best_spliced_length;
}

void SRAssembler::do_gene_finding() {
	GeneFinder* gene_finder = get_gene_finder();
	if (gene_finder == NULL)
		return;
	if (!gene_finder->is_available())
		return;
	logger->running("Now running the ab initio gene finding program ...");
	string program_name = gene_finder->get_program_name();
	Params params = this->get_parameters(program_name);
	gene_finder->do_gene_finding(this->hit_contigs_file, this->species, params, this->gene_finding_output_file, this->gene_finding_output_protein_file);
	logger->running("Done with gene finding.");
}

string SRAssembler::get_assembly_file_name(int round, int k){
	return get_assembler()->get_output_contig_file_name(aux_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round));
}

string SRAssembler::get_assembled_scaffold_file_name(int round, int k){
	return get_assembler()->get_output_scaffold_file_name(aux_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round));
}

string SRAssembler::get_spliced_alignment_file_name(int round, int k){
	return aux_dir + "/query-vs-contig_" + "k" + int2str(k) + "_" + "r" + int2str(round) + ".aln";
}

Aligner* SRAssembler::get_aligner(int round) {
	int aligner_type = Aligner::DNA_ALIGNER;
	if (round == 1 && probe_type == "protein")
		aligner_type = Aligner::PROTEIN_ALIGNER;
	Aligner* aligner = Aligner::getInstance(aligner_type, logger->get_log_level(), logger->get_log_file());
	return aligner;
}

Assembler* SRAssembler::get_assembler(){
	return Assembler::getInstance(this->assembler_program, logger->get_log_level(), logger->get_log_file());
}

SplicedAligner* SRAssembler::get_spliced_aligner(){
	return SplicedAligner::getInstance(this->spliced_alignment_program, logger->get_log_level(), logger->get_log_file());
}

GeneFinder* SRAssembler::get_gene_finder(){
	return GeneFinder::getInstance(this->gene_finding_program, logger->get_log_level(), logger->get_log_file());
}

Logger* SRAssembler::get_logger(){
	return logger;
}

void SRAssembler::create_index(int round) {
	// In case of pre-processing only, there is no probe to index.
	if (probe_file == "") {
		return;
	}
	Aligner* aligner = get_aligner(round);
	// This index is used during cleaning rounds, to find contigs that match the probe_file.
	aligner->create_index(aux_dir + "/qindex", probe_type, probe_file);
	// This index is used before walking to remove reads matching a protein taboo_file.
	if (taboo_file != "" && taboo_type == "protein") {
		aligner->create_index(aux_dir + "/tabooindex", taboo_type, taboo_file);
	}
}

string SRAssembler:: get_type(int round){
   return (round == 1)? probe_type: TYPE_DNA;
}

int SRAssembler::get_match_length(int round){
	return (round == 1)? init_match_length: recur_match_length;
}

int SRAssembler::get_mismatch_allowed(int round) {
	return (round == 1)? mismatch_allowed: 0;
}

string SRAssembler::get_vmatch_output_filename(int round, int lib_idx, int read_chunk){
	return tmp_dir + "/vmatch_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "chunk" + int2str(read_chunk);
}

void SRAssembler::merge_mapped_files(int round){
	for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
		Library lib = this->libraries[lib_idx];
		logger->debug("Now merging component matching reads files ...");
		string left_files = aux_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx + 1) + "_chunk*";
		string cmd = "cat " + left_files + " >> " + lib.get_matched_left_reads_filename();
		logger->debug(cmd);
		logger->fragile_run_shell_command(cmd);
		// Auxiliary files are removed here.
		cmd = "rm -f " + left_files;
		logger->debug(cmd);
		logger->fragile_run_shell_command(cmd);

		if (lib.get_paired_end()) {
			string right_files = aux_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx + 1) + "_chunk*";
			cmd = "cat " + right_files + " >> " + lib.get_matched_right_reads_filename();
			logger->debug(cmd);
			logger->fragile_run_shell_command(cmd);
			// Auxiliary files are removed here.
			cmd = "rm -f " + right_files;
			logger->debug(cmd);
			logger->fragile_run_shell_command(cmd);
		}
		logger->fragile_run_shell_command("cp " + lib.get_matched_left_reads_filename() + " " + lib.get_matched_left_reads_filename(round));
		if (lib.get_paired_end())
			logger->fragile_run_shell_command("cp " + lib.get_matched_right_reads_filename() + " " + lib.get_matched_right_reads_filename(round));
	}
	// Auxiliary files are removed here.
	string cmd = "rm -f " + get_contigs_index_name(round) + ".*";;
	logger->debug(cmd);
	logger->fragile_run_shell_command(cmd);
}

long SRAssembler::get_total_read_count(int round){
	long count = 0;
	for (unsigned short int lib_idx=0; lib_idx< this->libraries.size();lib_idx++) {
		Library lib = this->libraries[lib_idx];
		count += get_read_count(lib.get_matched_left_reads_filename(), FORMAT_FASTA);
		if (lib.get_paired_end()) {
			count += get_read_count(lib.get_matched_right_reads_filename(), FORMAT_FASTA);
		}
	}

	return count;
}

Assembly_stats SRAssembler::get_assembly_stats(int round, int k) {
	logger->debug("Counting longest contig length:" + get_assembly_file_name(round, k));
	ifstream contig_file(get_assembly_file_name(round, k).c_str());
	string line;
	vector<int> lens;
	Assembly_stats stats;
	stats.longest_contig = 0;
	stats.total_contig = 0;
	stats.n50 = 0;
	stats.n90 = 0;
	int total_length = 0;
	unsigned int contig_length = 0;
	while (getline(contig_file, line)) {
		if (line.substr(0,1) == ">") {
			if (contig_length > query_contig_min) {
				lens.push_back(contig_length);
				total_length += contig_length;
				stats.total_contig++;
				if (stats.longest_contig < contig_length)
					stats.longest_contig = contig_length;
			}
			contig_length = 0;
		} else
			contig_length += line.length();
	}
	if (contig_length > query_contig_min) {
		total_length += contig_length;
		stats.total_contig++;
		lens.push_back(contig_length);
		if (stats.longest_contig < contig_length)
			stats.longest_contig = contig_length;
	}
	contig_file.close();
	sort(lens.begin(), lens.end());
	int cumu_len = 0;
	for (short i=lens.size()-1;i>=0;i--){
		cumu_len += lens[i];
		if (cumu_len >= total_length/2.0 && stats.n50 == 0)
			stats.n50 = lens[i];
		if (cumu_len >= (double)total_length*0.9){
			stats.n90 = lens[i];
			break;
		}
	}
	return stats;
}

void SRAssembler::save_found_reads(int round){
	string matched_file = get_matched_reads_file_name(round);
	ofstream matched_file_stream;
	// Open the file in append mode.
	matched_file_stream.open(matched_file.c_str(), ios_base::app);
	for (unordered_set<string>::iterator it = found_reads.begin();it != found_reads.end(); ++it) {
		matched_file_stream << string(*it) + '\n';
	}
	matched_file_stream.close();
}

void SRAssembler::load_found_reads(int round){
	string matched_file = get_matched_reads_file_name(round);
	logger->info("Loading the matched reads of round " + int2str(round) + " (" + int2str(rank) + "/" + int2str(mpiSize-1) + ")");
	ifstream matched_file_stream(matched_file.c_str());
	string seq_id;
	int read_count;
	while (getline(matched_file_stream, seq_id)){
		found_reads.insert(seq_id);
		read_count++;
	}
	matched_file_stream.close();
	logger->info("Found matched read (pairs) from round " + int2str(round) + ": " + int2str(read_count));
}

SRAssembler* SRAssembler::getInstance(int pid){
	if (pid == 0){
		_srassembler = new SRAssemblerMaster();
		return _srassembler;
	} else {
		_srassembler = new SRAssemblerSlave();
		return _srassembler;
	}
	return NULL;
}

/* Unlike other DNA vs DNA searches in SRAssembler, we are indexing the contig queries, not the reads.
 * This is so that we can set vmatch to do a complete match, mapping the full read length to the contig.
 */
void SRAssembler::remove_unmapped_reads(unsigned int lib_idx, int round){
	string cmd;
	string contig_file = get_contig_file_name(round);
	// Index for good contigs was created in remove_no_hit_contigs().
	string contig_index = aux_dir + "/cindex";
	Aligner* aligner = get_aligner(round);
	Library lib = this->libraries[lib_idx];

	// Index current matched reads for extraction purposes.
	string left_matched_reads = lib.get_matched_left_reads_filename();
	string left_reads_index = aux_dir + "/left_lib" + int2str(lib_idx + 1) + "_index";
	aligner->create_index(left_reads_index, "dna", left_matched_reads);
	string right_matched_reads;
	string right_reads_index;
	if (lib.get_paired_end()) {
		right_matched_reads = lib.get_matched_right_reads_filename();
		right_reads_index = aux_dir + "/right_lib" + int2str(lib_idx + 1) + "_index";
		aligner->create_index(right_reads_index, "dna", right_matched_reads);
	}

	// Use the found reads as queries against the contigs to identify matching reads.
	string program_name = aligner->get_program_name();
	program_name += "_reads_vs_contigs";
	Params params = get_parameters(program_name);
	string vmatch_outfile = aux_dir + "/reads_vs_contigs.lib" + int2str(lib_idx+1) + ".round" + int2str(round) + ".vmatch";
	aligner->do_alignment(contig_index, "reads", 0, 2, left_matched_reads, params, vmatch_outfile);
	if (lib.get_paired_end()) {
		aligner->do_alignment(contig_index, "reads", 0, 2, right_matched_reads, params, vmatch_outfile);
	}

	// Use vseqselect to collect matching reads.
	logger->debug("Removing reads in library " + int2str(lib_idx + 1) + " without hits against contigs in round " + int2str(round));
	cmd = "vseqselect -seqnum " + vmatch_outfile + " " + left_reads_index + " | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + left_matched_reads;
	logger->debug(cmd);
	logger->safe_run_shell_command(cmd);
	cmd = "cp " + left_matched_reads + " " + lib.get_matched_left_reads_filename(round);
	logger->safe_run_shell_command(cmd);
	if (lib.get_paired_end()) {
		cmd = "vseqselect -seqnum " + vmatch_outfile + " " + right_reads_index + " | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + right_matched_reads;
		logger->debug(cmd);
		logger->safe_run_shell_command(cmd);
		cmd = "cp " + right_matched_reads + " " + lib.get_matched_right_reads_filename(round);
		logger->safe_run_shell_command(cmd);
	}
	// Auxiliary files are removed here.
	cmd = "rm -f " + vmatch_outfile + " " + left_reads_index + "*";
	logger->fragile_run_shell_command(cmd);
	if (lib.get_paired_end()) {
		cmd = "rm -f " + right_reads_index + "*";
		logger->fragile_run_shell_command(cmd);
	}
}

void finalized(){
	mpi_finalize();
}

int main(int argc, char * argv[] ) {
	long int start_time = time(0);
	int rank, mpiSize;

	SRAssembler* instance = NULL;
	try {
		mpi_init(argc,argv);
		mpiSize=mpi_get_size();
		rank=mpi_get_rank();

		instance = SRAssembler::getInstance(rank);
		int ret = instance->init(argc, argv, rank, mpiSize);
		if (ret == -1) {
			throw -1;
		}
		instance->do_preprocessing();
		instance->remove_taboo_reads();
		instance->do_walking();
	} catch (int e) {
		mpi_code code;
		code.action = ACTION_EXIT;
		code.value1 = 0;
		code.value2 = 0;
		mpi_bcast(get_mpi_code_value(code));
		finalized();
		return -1;
	}
	finalized();
	if (rank == 0) {
		string str = "Execution time: " + int2str(time(0) - start_time) + " seconds";
		instance->get_logger()->info(str);
	}
	return 0;
}
