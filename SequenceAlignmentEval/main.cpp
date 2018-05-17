#include <cstdlib>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "Evaluation.h"

using namespace std;
using namespace boost::program_options;
using namespace boost::filesystem;

int main(int argc,char ** argv)
{
	options_description desc;
	string genbank;
	string fasta;
	string output;
	desc.add_options()
		("help,h","print the current message")
		("reference,r",value<string>(&genbank),"the genbank file being used as a reference")
		("dna,d",value<string>(&fasta),"the sequence alignment result in FASTA format which is being evaluated")
		("output,o",value<string>(&output)->default_value("evaluation.txt"),"where the evaluation report will be output");
	variables_map vm;
	store(parse_command_line(argc,argv,desc),vm);
	notify(vm);

	if(1 == argc || vm.count("help")) {
		cout<<desc;
		return EXIT_SUCCESS;
	}

	if(1 != vm.count("reference")) {
		cerr<<"the genbank file must be given"<<endl;
		return EXIT_FAILURE;
	}

	if(1 != vm.count("dna")) {
		cerr<<"the DNA sequence being evaluated must be given"<<endl;
		return EXIT_FAILURE;
	}

	path outputfile(output);
	if(exists(outputfile)) {
		cerr<<"the output path already exists, won't override it!"<<endl;
		return EXIT_FAILURE;
	}

	Evaluation eval(genbank);
	bool retval = eval.evaluate(fasta,output);
	if(false == retval) {
		cerr<<"evaluation failed"<<endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}

