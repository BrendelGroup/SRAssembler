/*
 * SRAssemblerMaster.h
 *
 *  Created on: Oct 13, 2011
 *      Author: hchou
 */

#ifndef SRASSEMBLERMASTER_H_
#define SRASSEMBLERMASTER_H_

#include "SRAssembler.h"
#include <vector>

class SRAssemblerMaster: public SRAssembler {
public:
	SRAssemblerMaster();
	int init(int argc, char * argv[], int rank, int mpiSize);
	void show_usage();
	void print_message(const string&);
	int do_assembly(int round);
	void do_preprocessing();
	void do_walking();
	virtual ~SRAssemblerMaster();
private:
	void process_long_contigs(int round, int k);
	void prepare_final_contigs_file(int round);
	void remove_hit_contigs(vector<string> &contig_list, int round);
	void remove_no_hit_contigs(int round);
	int get_start_round();
	void load_long_contigs();
	void remove_contigs_no_hits(int round);
	void remove_unmapped_reads(int round);
	void clean_tmp_files(int round);
	void create_folders();
	void save_query_list();
	void load_query_list();
	void load_saved_contigs();
	void output_header();
	void output_libraries();
	void output_summary_header();
	void output_summary(int);
	void output_spliced_alignment();
	void output_gene_finding();
	void get_query_list();
	string get_saved_contig_file_name();
	int contig_number;
	unsigned best_k;
	string output_content;
	string summary_header;
	string summary_n50;
	string summary_n90;
	string summary_max;
	string summary_total;
	string summary_best;
	//vector<string> long_contigs;
	vector<string> query_list;
};

#endif /* SRASSEMBLERMaster_H_ */
