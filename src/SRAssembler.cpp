/*
 * SRAssembler.cpp
 *
 *  Created on: Oct 12, 2011
 *      Author: Hsien-chao Chou
 */

#include <ctime>
#include "SRAssembler.h"
#include "SRAssemblerMaster.h"
#include "SRAssemblerSlave.h"
#include "SRAssemblerSlave.h"
#include <vector>


SRAssembler* SRAssembler::_srassembler = NULL;

SRAssembler::SRAssembler() {
	// TODO Auto-generated constructor stub
	//cerr << "SRAssembler constructed." << endl ;
}

SRAssembler::~SRAssembler() {
	// TODO Auto-generated destructor stub
	cerr << "SRAssembler destructed" << endl ;
}

int SRAssembler::init(int argc, char * argv[], int rank, int mpiSize) {
//cerr << "I am an instance and my PPID is " + int2str(getppid()) + ".\n";
	//set the default values
	init_match_length = INIT_MATCH_LENGTH_PROTEIN;
	recur_match_length = RECUR_MATCH_LENGTH;
	mismatch_allowed = MISMATCH_ALLOWED_PROTEIN;
	num_rounds = NUM_ROUNDS;
	verbose = VERBOSE;
	assembly_round = ASSEMBLY_ROUND;
	reads_per_file = READS_PER_FILE;
	fastq_format = FASTQ_INTERLEAVED;
	out_dir = OUT_DIR;
	data_dir = "";
	type = QUERY_TYPE;
	assembler_program = ASSEMBLER_PROGRAM;
	gene_finding_program = GENE_FINDING_PROGRAM;
	spliced_alignment_program = SPLICED_ALIGNMENT_PROGRAM;
	species = DEFAULT_SPECIES;
	start_k = START_K;
	end_k = END_K;
	step_k = STEP_K;
	clean_round = CLEAN_ROUND;
	over_write = OVER_WRITE;
	check_gene_assembled = CHECK_GENE_ASSEMBLED;
	preprocessing_only = PREPROCESSING_ONLY;
	min_score = MIN_SCORE;
	min_coverage = MIN_COVERAGE;
	ini_contig_size = INI_CONTIG_SIZE;
	min_contig_lgth = MIN_CONTIG_LGTH;
	max_contig_lgth = MAX_CONTIG_LGTH;
	bool k_format = true;
	this->rank = rank;
	this->mpiSize = mpiSize;
	unsigned int insert_size = INSERT_SIZE;
	string left_read = "";
	string right_read = "";
	//parse command line options
	bool default_e = true;
	bool default_i = true;
	preprocessed_exist = false;
	usage = "Usage:\n";
	usage.append("\n");
	usage.append("    SRAssembler [options] -q query_file -p parameter_file\n");
	usage.append("                [-l library_file or -1 reads_file1 -2 reads_file2]\n");
	usage.append("\n");
	usage.append("-q: Required; FASTA-formatted query file.\n");
	usage.append("-t: Query file type; options: 'protein', 'cdna' [Default: " + QUERY_TYPE + "].\n");
	usage.append("-p: Required; SRAssembler parameter configuration file.\n\n");

	usage.append("-l: Required if the -1 option is not used; reads library file.\n");
	usage.append("-1: Required if the -l option is not used; use this option to specify the single-end reads file\n");
	usage.append("    or the left-end reads file for paired-end reads.\n");
	usage.append("-2: Right-end reads file for paired-end reads.\n");
	usage.append("-z: Insert size of the paired-end reads [Default: " + int2str(INSERT_SIZE) + "].\n");
	usage.append("-r: Directory in which to store or from which to retrieve the pre-processed reads [Default: output directory/" + READS_DATA + "].\n");
	usage.append("-o: SRAssembler output directory [Default: current directory].\n\n");

	usage.append("-P: Run the read pre-processing step only, then terminate SRAssembler.\n");
	usage.append("-x: Number of reads per pre-preprocessed reads file [Default: " + int2str(READS_PER_FILE) + "].\n");
	usage.append("-A: Assembler program choice; options: 0=>SOAPdenovo2, 1=>ABySS [Default: " + int2str(ASSEMBLER_PROGRAM) + "].\n");
	usage.append("-k: Specifies the k-mer set to be used by the assembler; format: start_k:interval:end_k.\n");
	usage.append("    Start_k and end_k must be odd integers, and interval must be an even integer, similar to the following example:\n");
	usage.append("    '15:10:45' specifies that k-mer values 15, 25, 35, 45 will be tested.[Default: " + int2str(START_K) + ":" + int2str(STEP_K) + ":" + int2str(END_K) + "].\n");
	usage.append("-S: Spliced alignment program; options: 0=>GeneSeqer, 1=>GenomeThreader,\n");
	usage.append("    2=>Exonerate [Default: " + int2str(SPLICED_ALIGNMENT_PROGRAM) + "].\n");
	usage.append("-s: Species model for spliced alignment; options (for GenomeThreader and GeneSeqer):\n");
	usage.append("    'human', 'mouse', 'rat', 'chicken', 'drosophila', 'nematode', 'fission_yeast', 'aspergillus', 'arabidopsis',\n");
	usage.append("    'maize', 'rice', 'medicago' [DEFAULT: " + DEFAULT_SPECIES + "].\n");
	usage.append("-G: Ab initio gene finding program; options: 0=>None, 1=>Snap [Default: " + int2str(GENE_FINDING_PROGRAM) + "].\n\n");

	usage.append("-i: Initial contig size for chromosome walking [Default: " + int2str(INI_CONTIG_SIZE) + "].\n");
	usage.append("-m: Minimum contig length to be reported [Default: " + int2str(MIN_CONTIG_LGTH) + "].\n");
	usage.append("-M: Maximum contig length to be reported [Default: " + int2str(MAX_CONTIG_LGTH) + "].\n");
	usage.append("-e: Minimum spliced alignment score for hits [Default: " + double2str(MIN_SCORE) + "].\n");
	usage.append("-c: Minimum spliced alignment coverage score for hits [Default: " + double2str(MIN_COVERAGE) + "].\n\n");

	usage.append("-n: Maximum number of rounds for chromosome walking [Default: " + int2str(NUM_ROUNDS) + "].\n");
	usage.append("-a: The number of the round in which to start read assembly [Default: " + int2str(ASSEMBLY_ROUND) + "].\n");
	usage.append("-b: The number of the round in which to periodically remove unrelated contigs and reads. For example,\n");
	usage.append("    “-b 3” specifies that SRAssembler will remove unrelated contigs and reads after assembly rounds 3, 6, 9, ... [Default: " + int2str(CLEAN_ROUND) + "].\n");
	usage.append("-w: Forgo spliced alignment check after intermediate assembly rounds [SRAssembler will continue for the -n specified number of rounds].\n");
	usage.append("-y: Disable SRAssembler resumption from previous checkpoint [will overwrite existing output].\n\n");

	usage.append("-v: Verbose output.\n");
	usage.append("-h: Print this usage synopsis.");


	char c;
	while((c = getopt(argc, argv, "q:t:p:l:1:2:z:r:o:Px:A:k:S:s:G:i:m:M:e:c:n:a:b:wyvh")) != -1) {
		 switch (c){
			 case 'q':
				 query_file = optarg;
				 break;
			 case 't':
				 type = optarg;
				 break;
			 case 'p':
			 	 param_file = optarg;
			 	 break;
			 case 'l':
				 library_file = optarg;
				 break;
			 case '1':
				 left_read = optarg;
				 break;
			 case '2':
				 right_read = optarg;
				 break;
			 case 'z':
				 insert_size = str2int(optarg);
				 break;
			 case 'r':
			 	 data_dir = optarg;
			 	 break;
			 case 'o':
				 out_dir = optarg;
				 break;
			 case 'P':
			 	 preprocessing_only = 1;
			 	 break;
			 case 'x':
				 reads_per_file = str2int(optarg);
				 break;
			 case 'A':
			 	 assembler_program = str2int(optarg);
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
			 case 'S':
				 spliced_alignment_program = str2int(optarg);
				 break;
			 case 's':
				 species = optarg;
				 break;
			 case 'G':
			 	 gene_finding_program = str2int(optarg);
			 	 break;
			 case 'i':
				 ini_contig_size = str2int(optarg);
				 break;
			 case 'm':
				 min_contig_lgth = str2int(optarg);
				 break;
			 case 'M':
				 max_contig_lgth = str2int(optarg);
				 break;
			 case 'e':
				 min_score = str2double(optarg);
				 break;
			 case 'c':
				 min_coverage = str2double(optarg);
				 break;
			 case 'n':
				 num_rounds = str2int(optarg);
				 break;
			 case 'a':
				 assembly_round = str2int(optarg);
				 if (assembly_round == 0)
					 assembly_round = 9999;
				 break;
			 case 'b':
				 clean_round = str2int(optarg);
				 break;
			 case 'w':
				 check_gene_assembled = 0;
				 break;
			 case 'y':
			 	 over_write = 1;
			 	 break;
			 case 'v':
				 verbose = 1;
				 break;
			 case 'h':
				 show_usage();
				 return -1;
			 	 break;
			 case '?':
				 string msg = "input error! unknown option : -";
				 msg.append(1,optopt);
				 print_message(msg);         //unknown options
				 return -1;
		  }
	}
	if (argc == 1){
		show_usage();
		return -1;
	}
	tmp_dir = out_dir + "/tmp";
	//mem_dir = "/dev/shm/SRAssembler" + int2str(getppid());
	if (data_dir == "")
		data_dir = out_dir + "/" + READS_DATA;
	preprocessed_exist = file_exists(data_dir);
	results_dir = out_dir + "/output";
	intermediate_dir = results_dir + "/intermediates";
	log_file = results_dir + "/msg.log";
	spliced_alignment_output_file = results_dir + "/output.aln";
	gene_finding_output_file = results_dir + "/output.ano";
	gene_finding_output_protein_file = results_dir + "/snap.predicted.prot";
	//spliced_alignment_gsq_file = results_dir + "/output.gsq";
	final_contigs_file = results_dir + "/all_contigs.fasta";
	summary_file = results_dir + "/summary.html";
	hit_contigs_file = results_dir + "/hit_contigs.fasta";

	int level = Logger::LEVEL_INFO;
	if (verbose)
		level = Logger::LEVEL_DEBUG;
	logger = Logger::getInstance(level, log_file);


	if (param_file != ""){
		if (!file_exists(param_file)){
			print_message("Parameter file : " + param_file + " does not exist!");
			return -1;
		//} else {
		//parameter_dict = read_param_file();
		}
	} else {
		print_message("Parameter file required!");
		return -1;
	}

	if (library_file != ""){
		if (!file_exists(library_file)){
			print_message("library file: " + library_file + " does not exist!");
			return -1;
		}
		if (!read_library_file())
			return -1;
	} else {
		if (left_read == "" && !preprocessed_exist){
			print_message("-1 or -P or library file is required");
			return -1;
		}
		Library lib(0, this->data_dir, this->tmp_dir,this->logger);
		lib.set_format(FORMAT_FASTQ);
		lib.set_left_read(left_read);
		if (right_read != ""){
			lib.set_paired_end(true);
			lib.set_right_read(right_read);
			lib.set_insert_size(insert_size);
			if (!file_exists(right_read) && !preprocessed_exist){
				print_message("file: " + right_read + " does not exist!");
				return -1;
			}
		}
		this->libraries.push_back(lib);
	}

	//check missing options
	if (num_rounds <= 0){
		print_message("-n must be larger than 0");
		return -1;
	}
	if (reads_per_file <= 1000){
		print_message("-x value is not valid or too small");
		return -1;
	}

	if (type != TYPE_PROTEIN && type != TYPE_CDNA){
		print_message("-t must be 'protein' or 'cdna'");
		return -1;
	} else {
		if (default_e && type == TYPE_CDNA)
			mismatch_allowed = MISMATCH_ALLOWED_CDNA;
		if (default_i && type == TYPE_CDNA)
			init_match_length = INIT_MATCH_LENGTH_CDNA;
	}
	if (species != "human" && species != "mouse" && species != "rat" && species != "chicken" && species != "drosophila" && species != "nematode" && species != "fission_yeast" && species != "aspergillus" && species != "arabidopsis" && species != "maize" && species != "rice" && species != "medicago"){
		print_message("-s is required");
		return -1;
	}
	if (!k_format){
		print_message("-k format error. The format is : start_k:interval:end_k. The start_k and end_k must be odd value, and The interval must be even value.");
		return -1;
	}
	//check the existence of read files
	if (!file_exists(query_file)){
		print_message("file: " + query_file + " does not exist!");
		return -1;
	}

	if (!file_exists(out_dir)){
		print_message("output directory: " + out_dir + " does not exist!");
		return -1;
	}

	return 0;
}

//TODO this should read through one time and store all of the potential parameters so that we aren't always rereading the file
//boost::unordered_map<std::string,Params> SRAssembler::read_param_file() {
	//ifstream param_file(this->param_file.c_str());
	//string line;
	//Params params;
	//string program_name;
	//while (getline(param_file, line)){
		//line = trim(line);
		//if (line.length() == 0) continue;
		//if (line.substr(0,1) == "[") {
			//program_name = line.substr(1, line.length() - 2)
			//found_program = true;
			//continue;
		//}
		//if (line.substr(0,1) == "#") continue;
		//if (found_program && line.substr(0,1) == "[") break;
		//if (found_program){
			//vector<string> tokens;
			//tokenize(line, tokens, "=");
			//if (tokens.size() == 2){
				//string param = trim(tokens[0]);
				//string value = trim(tokens[1]);
				//params.insert(Params::value_type(param, value));
			//}
		//}
	//}
	//return params;
//}

//Params SRAssembler::get_parameters(string program_name) {
//}
Params SRAssembler::read_param_file(string program_name) {
	ifstream param_file(this->param_file.c_str());
	string line;
	Params params;
	bool found_program = false;
	while (getline(param_file, line)){
		line = trim(line);
		if (line.length() == 0) continue;
		if (line == ("[" + program_name + "]")) {
			found_program = true;
			continue;
		}
		if (line.substr(0,1) == "#") continue;
		if (found_program && line.substr(0,1) == "[") break;
		if (found_program){
			vector<string> tokens;
			tokenize(line, tokens, "=");
			if (tokens.size() == 2){
				string param = trim(tokens[0]);
				string value = trim(tokens[1]);
				params.insert(Params::value_type(param, value));
			}
		}
	}
	return params;
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
					print_message("r1 file is expected in config file!");
					return false;
				}
				this->libraries.push_back(*lib);
			}
			lib = new Library(this->libraries.size(), this->data_dir, this->tmp_dir, this->logger);
		} else {
			vector<string> tokens;
			tokenize(line, tokens, "=");
			if (tokens.size() == 2){
				string param = trim(tokens[0]);
				string value = trim(tokens[1]);
				if (param == "r1" && lib != NULL) {
					lib->set_left_read(value);
					if (!file_exists(lib->get_left_read()) && !preprocessed_exist) {
						print_message("r1 file in config file: " + lib->get_left_read() + " does not exist!");
						return false;
					}
				}
				if (param == "r2" && lib != NULL) {
					lib->set_right_read(value);
					if (!file_exists(lib->get_right_read()) && !preprocessed_exist) {
						print_message("r2 file in config file: " + lib->get_right_read() + " does not exist!");
						return false;
					}
					lib->set_paired_end(true);
				}

				if (param == "insert_size" && lib != NULL) {
					lib->set_insert_size(str2int(value));
				}
				if (param == "direction" && lib != NULL) {
					lib->set_reversed(str2int(value));
				}
				if (param == "format" && lib != NULL) {
					if (value == "fastq") {
						lib->set_format(FORMAT_FASTQ);
					}
					else if (value == "fasta") {
						lib->set_format(FORMAT_FASTA);
					} else {
						print_message("format in config file should be 'fastq' or 'fasta'!");
						return false;
					}
				}
			}
		}
	}
	if (lib != NULL){
		if (lib->get_left_read() == ""){
			print_message("r1 file is expected in config file!");
			return false;
		}
		this->libraries.push_back(*lib);
	}
	if (this->libraries.size() == 0){
		print_message(" No [LIBRARY] section found in config file!");
		return false;
	}
	return true;

}

int SRAssembler::get_file_count(string search_pattern){
	string cmd = "ls -l " + search_pattern + " | wc -l";
	logger->debug(cmd);
	return str2int(run_shell_command_with_return(cmd));
}

int SRAssembler::count_preprocessed_reads(int lib_idx){
	// This uses a system call to count the lines in all the fasta files in the split reads directory
	string cmd = "wc -l " + data_dir + "/lib" + int2str(lib_idx+1) + "/*part*.fasta | tail -n 1 | cut -d' ' -f3";
	logger->debug(cmd);
	return str2int(run_shell_command_with_return(cmd)) / 2;
}

void SRAssembler::preprocess_read_part(int lib_idx, int read_part){
	Library lib = this->libraries[lib_idx];
	logger->info("preprocessing lib " + int2str(lib_idx + 1) + ", reads file (" + int2str(read_part) + "/" + int2str(lib.get_num_parts()) + ")");
	string suffix = int2str(read_part) + ".fasta";
	string left_read_file = lib.get_split_read_prefix(lib.get_left_read()) + suffix;
	string right_read_file;
	if (lib.get_paired_end())
		right_read_file = lib.get_split_read_prefix(lib.get_right_read()) + suffix;
	// Create the Vmatch mkvtree indexes
	Aligner* aligner = get_aligner(0);  // Round 0 means DNA Aligner
	// Creating both indices in one function is not maximally efficient when using a lot of processors, but it's not terrible.
	aligner->create_index(lib.get_read_part_index_name(read_part, LEFT_READ), "dna", left_read_file);
	if (lib.get_paired_end())
		aligner->create_index(lib.get_read_part_index_name(read_part, RIGHT_READ), "dna", right_read_file);
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
//TODO This may become unnecessary
string SRAssembler:: get_contigs_index_name(int round){
   return tmp_dir + "/round" + int2str(round);
}

/* It appears that this function is deceptively named.
 * It does not just return a string naming a fasta file that contains all of the matched reads.
 * This function is also responsible for assembling the contents of that file.
 */
 //TODO break into multiple functions, this is ridiculous
string SRAssembler:: get_query_fasta_file_name(int round){
	if (round > 1){
		if (assembly_round < round)
			return get_contig_file_name(round-1);
		else {
			string joined_file = tmp_dir + "/matched_reads_joined.fasta";
			string cmd;
			// For each library, append the matched reads to the matched_reads_joined.fasta
			for (unsigned i=0;i<this->libraries.size();i++){
				string left_file = tmp_dir + "/matched_reads_left_" + "lib" + int2str(i+1) + ".fasta";
				string right_file = tmp_dir + "/matched_reads_right_" + "lib" + int2str(i+1) + ".fasta";
				if (libraries[i].get_paired_end()){
					cmd = "cat " + left_file + " " + right_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				} else {
					cmd = "cat " + left_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				}
			}
			return joined_file;
		}
	}
	return query_file;
}

void SRAssembler::mask_contigs(int round){
	string cmd;
	string contig_file;
	contig_file = get_contig_file_name(round);
	string masked_file = contig_file + ".masked";
	// Sed command replaces lowercase letters NOT in the header with capital Ns. AWK command converts FASTA to single-line.
	cmd = "dustmasker -in " + contig_file + " -outfmt fasta -out - | sed '/^[^>]/s/[a-z]/N/g' | awk '!/^>/ { printf \"%s\", $0; n = \"\\n\" } /^>/ { print n $0} END { printf n }' > " + masked_file;
	run_shell_command(cmd);
}

string SRAssembler:: get_query_fasta_file_name_masked(int round){
	if (round > 1){
		if (assembly_round < round) {
			string masked_fasta = get_contig_file_name(round-1) + ".masked";
			return masked_fasta;
		} else {
			string joined_file = tmp_dir + "/matched_reads_joined.fasta";
			string cmd;
			for (unsigned i=0;i<this->libraries.size();i++){
				string left_file = tmp_dir + "/matched_reads_left_" + "lib" + int2str(i+1) + ".fasta";
				string right_file = tmp_dir + "/matched_reads_right_" + "lib" + int2str(i+1) + ".fasta";
				if (libraries[i].get_paired_end()){
					cmd = "cat " + left_file + " " + right_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				} else {
					cmd = "cat " + left_file + " >> " + joined_file;
					logger->debug(cmd);
					run_shell_command(cmd);
				}
			}
			return joined_file;
		}
	}
	return query_file;
}

string SRAssembler:: get_contig_file_name(int round){
	return intermediate_dir + "/contigs_" + "r" + int2str(round) + ".fasta";
}

string SRAssembler:: get_matched_reads_file_name(int round){
	//return tmp_dir + "/matched_reads_" + "r" + int2str(round) + "_" + "rank" + int2str(this->rank) + ".list";
	//return tmp_dir + "/matched_reads_" + "r" + int2str(round) + "_" + "rank" + "0" + ".list";
	return tmp_dir + "/found_reads_" + "r" + int2str(round) + ".list";
}

int SRAssembler::do_alignment(int round, int lib_idx, int read_part) {
	Library lib = this->libraries[lib_idx];
	logger->info("Aligning: round = " + int2str(round) + " Lib (" + int2str(lib_idx+1) + "/" + int2str(this->libraries.size()) + "), Reads (" + int2str(read_part) + "/" + int2str(lib.get_num_parts()) + ")");
	Aligner* aligner = get_aligner(round);
	string program_name = aligner->get_program_name();
	if (round == 1) {
		program_name += "_" + get_type(1) + "_init";
	} else {
		program_name += "_extend_contig";
	}
	logger->debug("... using Vmatch criteria: " + program_name);
	//TODO we read this parameter file A LOT. We should import the parameters for each program_name once.
	Params params = this->read_param_file(program_name);
	int new_read_count;
	// Reads as queries are necessary when searching against a protein.
	if (round == 1 && type == "protein") {
		aligner->do_alignment(get_contigs_index_name(round), get_type(round), get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_part, LEFT_READ), params, get_vmatch_output_filename(round, lib_idx, read_part));
		if (lib.get_paired_end())
			aligner->do_alignment(get_contigs_index_name(round), get_type(round), get_match_length(round), get_mismatch_allowed(round), lib.get_split_file_name(read_part, RIGHT_READ), params, get_vmatch_output_filename(round, lib_idx, read_part));
		new_read_count = aligner->parse_output(get_vmatch_output_filename(round, lib_idx, read_part), found_reads, lib_idx, read_part, lib.get_read_part_index_name(read_part, LEFT_READ), lib.get_read_part_index_name(read_part, RIGHT_READ), lib.get_matched_left_reads_filename(round, read_part), lib.get_matched_right_reads_filename(round, read_part));
		//save_found_reads(round);
		return new_read_count;
	// After round 1 we use the masked contig file as the query
	} else {
		aligner->do_alignment(lib.get_read_part_index_name(read_part, LEFT_READ), get_type(round), get_match_length(round), get_mismatch_allowed(round), get_query_fasta_file_name_masked(round), params, get_vmatch_output_filename(round, lib_idx, read_part));
		if (lib.get_paired_end())
			aligner->do_alignment(lib.get_read_part_index_name(read_part, RIGHT_READ), get_type(round), get_match_length(round), get_mismatch_allowed(round), get_query_fasta_file_name_masked(round), params, get_vmatch_output_filename(round, lib_idx, read_part));
		new_read_count = aligner->parse_output(get_vmatch_output_filename(round, lib_idx, read_part), found_reads, lib_idx, read_part, lib.get_read_part_index_name(read_part, LEFT_READ), lib.get_read_part_index_name(read_part, RIGHT_READ), lib.get_matched_left_reads_filename(round, read_part), lib.get_matched_right_reads_filename(round, read_part));
		//save_found_reads(round);
		return new_read_count;
	}
}

void SRAssembler::do_assembly(int round, int k) {
	logger->info("Doing assembly: round = " + int2str(round) + ", k = " + int2str(k));
	Assembler* assembler = get_assembler();
	assembler->do_assembly(k, this->libraries, tmp_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round));
}

void SRAssembler::do_spliced_alignment() {
	logger->info("Doing the final spliced alignment ...");
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	string program_name = spliced_aligner->get_program_name();
	Params params = this->read_param_file(program_name);
	spliced_aligner->do_spliced_alignment(this->final_contigs_file, type, this->query_file, this->species, params, this->spliced_alignment_output_file);
	spliced_aligner->get_hit_contigs(min_score, min_coverage, min_contig_lgth, this->final_contigs_file, this->hit_contigs_file, this->spliced_alignment_output_file);
	logger->info("Done.");
}

string_map SRAssembler::do_spliced_alignment(int round) {
	logger->info("Now running the spliced alignment program ...");
	SplicedAligner* spliced_aligner = get_spliced_aligner();
	string program_name = spliced_aligner->get_program_name();
	Params params = this->read_param_file(program_name);
	string contig_file = get_contig_file_name(round);
	string output_file = tmp_dir + "/query-vs-contig_" + "r" + int2str(round) + ".aln";
	string hit_file = tmp_dir + "/hit_contigs_" + "r" + int2str(round) + ".fasta";
	spliced_aligner->do_spliced_alignment(contig_file, type, this->query_file, this->species, params, output_file);
	string_map query_map = spliced_aligner->get_aligned_contigs(min_score, min_coverage, min_contig_lgth, contig_file, hit_file, output_file);
	//RM HERE
	spliced_aligner->clean_files(contig_file);
	logger->info("Done.");
	return query_map;
}

void SRAssembler::do_gene_finding() {
	GeneFinder* gene_finder = get_gene_finder();
	if (gene_finder == NULL)
		return;
	if (!gene_finder->is_available())
		return;
	logger->info("Now running the ab initio gene finding program ...");
	string program_name = gene_finder->get_program_name();
	Params params = this->read_param_file(program_name);
	gene_finder->do_gene_finding(this->hit_contigs_file, this->species, params, this->gene_finding_output_file, this->gene_finding_output_protein_file);
	logger->info("Done.");
}

string SRAssembler::get_assembly_file_name(int round, int k){
	return get_assembler()->get_output_contig_file_name(tmp_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round));
}

string SRAssembler::get_assembled_scaf_file_name(int round, int k){
	return get_assembler()->get_output_scaffold_file_name(tmp_dir + "/assembly_" + "k" + int2str(k) + "_" + "r" + int2str(round));
}

Aligner* SRAssembler::get_aligner(int round) {
	int aligner_type = Aligner::DNA_ALIGNER;
	if (round == 1 && type == "protein")
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
	Aligner* aligner = get_aligner(round);
	aligner->create_index(get_contigs_index_name(round), get_type(round), get_query_fasta_file_name_masked(round));
	aligner->create_index(tmp_dir + "/qindex", type, query_file);
}

string SRAssembler:: get_type(int round){
   return (round == 1)? type: TYPE_CDNA;
}

int SRAssembler::get_match_length(int round){
	return (round == 1)? init_match_length: recur_match_length;
}

int SRAssembler::get_mismatch_allowed(int round) {
	return (round == 1)? mismatch_allowed: 0;
}

string SRAssembler::get_vmatch_output_filename(int round, int lib_idx, int read_part){
	return mem_dir + "/vmatch_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx+1) + "_" + "part" + int2str(read_part);
}

void SRAssembler::merge_mapped_files(int round){
	for (unsigned int lib_idx=0;lib_idx<this->libraries.size();lib_idx++){
		Library lib = this->libraries[lib_idx];
		logger->debug("Now merging component matching reads files ...");
		string left_files = tmp_dir + "/matched_reads_left_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx + 1) + "_part*";
		string cmd = "cat " + left_files + " >> " + lib.get_matched_left_reads_filename();
		logger->debug(cmd);
		run_shell_command(cmd);
		//RM HERE
		cmd = "rm -f " + left_files;
		logger->debug(cmd);
		run_shell_command(cmd);

		if (lib.get_paired_end()) {
			string right_files = tmp_dir + "/matched_reads_right_" + "r" + int2str(round) + "_" + "lib" + int2str(lib_idx + 1) + "_part*";
			cmd = "cat " + right_files + " >> " + lib.get_matched_right_reads_filename();
			logger->debug(cmd);
			run_shell_command(cmd);
			//RM HERE
			cmd = "rm -f " + right_files;
			logger->debug(cmd);
			run_shell_command(cmd);
		}
		if (round > 1) {
			run_shell_command("cp " + lib.get_matched_left_reads_filename() + " " + lib.get_matched_left_reads_filename(round));
			if (lib.get_paired_end())
				run_shell_command("cp " + lib.get_matched_right_reads_filename() + " " + lib.get_matched_right_reads_filename(round));
		}
	}
	//RM HERE
	string cmd = "rm -f " + get_contigs_index_name(round) + ".*";;
	logger->debug(cmd);
	run_shell_command(cmd);
	logger->debug("done.");
}

int SRAssembler::get_total_read_count(int round){
	int count = 0;
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
			if (contig_length > ini_contig_size) {
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
	if (contig_length > ini_contig_size) {
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
	// Open the file in append mode
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
