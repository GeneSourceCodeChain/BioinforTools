#include <cmath>
#include <stdexcept>
#include <fstream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include "Evaluation.h"

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

Evaluation::Evaluation(string genbank)
{
	//load stop codon position from genbank which is supposed to be correct
    	if(false == exists(genbank)) throw runtime_error(genbank + " doesn't exists");
	stopCodonPos = getStopCodonPos(genbank);
}

Evaluation::~Evaluation()
{
}

bool Evaluation::evaluate(string fasta,string output)
{
	//evaluate the possibility of a sequence alignment result
	if(false == exists(fasta)) return false;
	std::ofstream out(output);
	if(false == out.is_open()) return false;
	string sequence = getDNASequence(fasta);
	long ORFMaxLen;
	ORFList ORFs;
	boost::tie(ORFMaxLen,ORFs) = getORF(sequence);

	auto foreground = count(sequence,ORFs,IsLongPeptideChain());
	auto background = count(sequence,ORFs,IsShortPeptideChain());

	for(long index = 3 ; index < ORFMaxLen ; index += 3) {
		//for every orf length
		double avg_log_odds = 0;
		long ORF_count = 0;	//how many orfs of length index
		long positive = 0;
		long gene_count = 0;
		long pos_and_real = 0;
		for(ORFList::nth_index<1>::type::const_iterator it = ORFs.get<1>().lower_bound(index) ; it != ORFs.get<1>().upper_bound(index) ; it++) {
			if(it->length == index) {
				double P_prob = 0;
				double Q_prob = 0;
				P_prob = logsum(P_prob,get<0>(foreground)[sequence[it->start_pos]]);
				P_prob = logsum(P_prob,get<1>(foreground)[sequence[it->start_pos]][sequence[it->start_pos + 1]]);
				P_prob = logsum(P_prob,get<2>(foreground)[sequence[it->start_pos]][sequence[it->start_pos + 1]][sequence[it->start_pos + 2]]);
				Q_prob = logsum(Q_prob,get<0>(background)[sequence[it->start_pos]]);
				Q_prob = logsum(Q_prob,get<1>(background)[sequence[it->start_pos]][sequence[it->start_pos + 1]]);
				Q_prob = logsum(Q_prob,get<2>(background)[sequence[it->start_pos]][sequence[it->start_pos + 1]][sequence[it->start_pos + 2]]);
				for(int ii = it->start_pos ; ii < it->start_pos + it->length - 3 ; ii++) {
					P_prob = logsum(P_prob,get<3>(foreground)[sequence[ii]][sequence[ii + 1]][sequence[ii + 2]][sequence[ii + 3]]);
					Q_prob = logsum(Q_prob,get<3>(background)[sequence[ii]][sequence[ii + 1]][sequence[ii + 2]][sequence[ii + 3]]);
				}
				double log_odds = log(P_prob / Q_prob);
				avg_log_odds += log_odds;
				ORF_count++;
				if(log_odds > 0) positive++;
				long stop_pos = it->start_pos + it->length - 1;
				//search the stop_pos in stopCodonPos
				if(stopCodonPos.end() != stopCodonPos.find(stop_pos)) {
					gene_count++;
					if(log_odds > 0) pos_and_real++;
				}
			}
		}//end for lower(index)->upper(index)
		avg_log_odds /= ORF_count;
		if(gene_count > 0) {
			out<<"ORF Len: "<<index
				<<", gene/ORF: "<<gene_count<<"/"<<ORF_count
				<<", real/pos: "<<pos_and_real<<"/"<<positive
				<<", avg log odds: "<<avg_log_odds<<endl;
		}
	}
	return true;
}

set<long> Evaluation::getStopCodonPos(string genbank)
{
	regex expr("([0-9]+)\\.\\.([0-9]+)");
	set<long> stopCodonPos;
	std::ifstream in(genbank);
	string str;
	while(false == in.eof()) {
		in >> str; trim(str);
		if(str == "CDS") {
			in >> str; trim(str);
			match_results<string::const_iterator> what;
			if(regex_match(str,what,expr)) {
				string from(what[1].begin(),what[1].end());
				string to(what[2].begin(),what[2].end());
				stopCodonPos.insert(lexical_cast<long>(to) - 1);
			}
		}
	}
	return stopCodonPos;
}

string Evaluation::getDNASequence(string fasta)
{
	std::ifstream in(fasta);
	string sequence;
	string line;
	while(false == in.eof()) {
		getline(in,line);
		if(line[0] != '>') sequence += line;
	}
	regex expr("[^GCAT]");
	string normalized_sequence;
	//replace all non CGAT alphabet with T
	normalized_sequence = regex_replace(sequence,expr,"T",boost::match_default|boost::format_all);
	
	return normalized_sequence;
}

boost::tuple<long,ORFList> Evaluation::getORF(string sequence)
{
	long ORFMaxLen = 0;
	ORFList ORFs;
	long ORF_starts[3] = {0,1,2};
	for(long index = 2 ; index < sequence.size() ; index++) {
		long codon_start = index - 2;
		string codon = sequence.substr(codon_start,3);
		if(codon == "TAA" || codon == "TAG" || codon == "TGA") {
			//if the stop codon is met
			long ORF_start = ORF_starts[codon_start % 3];
			long length = index - ORF_start + 1;
#ifndef NDEBUG
			assert(length % 3 == 0);
#endif
			ORFMaxLen = max(ORFMaxLen,length);
			ORFs.insert(ORF(ORF_start,length));
			//start from the next nucleotide
			ORF_starts[codon_start % 3] = index + 1;
		}
	}
	return boost::make_tuple(ORFMaxLen,ORFs);
}

double Evaluation::logsum(double lx,double ly)
{
	//log(a+b)=log(a(1+b/a))=log(a) + log(1 + b/a)
	//=log(a) + log(1+exp(log(b/a))) = log(a) + log(1 + exp(log(b) - log(a)))
	if(isnan(lx)) return ly;
	if(isnan(ly)) return lx;
	return (lx > ly)?(lx + log(1 + exp(ly - lx))):(ly + log(1 + exp(lx - ly)));
}

