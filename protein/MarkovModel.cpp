#include <stdexcept>
#include <ifstream>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include "MarkovModel.h"

using namespace std;
using namespace boost;
using namespace boost::algorithm;
using namespace boost::filesystem;

MarkovModel::MarkovModel(string genbank,string fasta)
{
    if(false == exists(genbank)) throw runtime_error(genbank + " doesn't exists");
	if(false == exists(fasta)) throw runtime_error(fasta + "doesn't exists");
	stopCodonPos = getStopCodonPos(genbank);
	sequence = getDNASequence(fasta);
	boost::tie(ORFMaxLen,ORFs) = getORF(sequence);
}

MarkovModel::~MarkovModel()
{
}

auto MarkovModel::count(unary_function<long,bool> & predicate)
{
	long singlesum = 0;
	long doublesum = 0;
	long triplesum = 0;
	long quadsum = 0;
	map<char,long> singlecount;
	map<char,map<char,long> > doublecount;
	map<char,map<char,map<char,long> > > triplecount;
	map<char,map<char,map<char,map<char,long> > > > quadcount;
	for(auto & ORF : ORFs) {
		if(predicate(ORF.second)) {
			//count how many singles doubles triples quad can be contained in the sequence
			singlesum += ORF.second;
			doublesum += ORF.second - 1;
			triplesum += ORF.second - 2;
			quadsum += ORF.second - 3;
			//count the occurrence of each type of nucleotide
			for(int index = ORF.first ; index < ORF.first + ORF.second ; index++) singlecount[sequence[index]]++;
			for(int index = ORF.first ; index < ORF.first + ORF.second - 1 ; index++) doublecount[sequence[index]][sequence[index + 1]]++;
			for(int index = ORF.first ; index < ORF.first + ORF.second - 2 ; index++) triplecount[sequence[index]][sequence[index + 1]][sequence[index + 2]]++;
			for(int index = ORF.first ; index < ORF.first + ORF.second - 3 ; index++) quadcount[sequence[index]][sequence[index + 1]][sequence[index + 2]][sequence[index + 3]]++;
		}
	}
	//normalize
	for(auto & s : singlecount) s.second = log(static_cast<double>(s.second) / singlesum);
	for(auto & s : doublecount)
		for(auto & d : s.second) d.second = log(static_cast<double>(d.second) / doublesum);
	for(auto & s : triplecount)
		for(auto & d : s.second)
			for(auto & t : d.second) t.second = log(static_cast<double>(t.second) / triplesum);
	for(auto & s : quadcount)
		for(auto & d : s.second)
			for(auto & t : d.second)
				for(auto & q : t.second) q.second = log(static_cast<double>(q.second) / quadsum);
	return boost::make_tuple(singlesum,doublesum,triplesum,quadsum,singlecount,doublecount,triplecount,quadcount);
}

set<long> MarkovModel::getStopCodonPos(string genbank)
{
	regex expr("([0-9]+)\\.\\.([0-9]+)");
	set<long> stopCodonPos;
	ifstream in(genbank);
	string str;
	while(false == in.eof()) {
		in >> str; trim(str);
		if(str == "CDS") {
			in >> str; trim(str);
			cmatch what;
			if(regex_match(str,what,expr)) {
				string from(what[1].begin(),what[1].end());
				string to(what[2].begin(),what[2].end());
				stopCodonPos.insert(lexical_cast<long>(to) - 1);
			}
		}
	}
	return stopCodonPos;
}

string MarkovModel::getDNASequence(string fasta)
{
	ifstream in(fasta);
	string sequence;
	string line;
	while(false == in.eof()) {
		getline(in,line);
		if(line[0] != '>') sequence += line;
	}
	regex expr("[^GCAT]");
	string normalized_sequence;
	//replace all non CGAT alphabet with T
	regex_replace(normalized_sequence,sequence.begin(),sequence.end(),expr,"T",boost::match_default|boost::format_all);
	
	return normalized_sequence;
}

boost::tuple<long,vector<pair<long,long> > > MarkovModel::getORF(string sequence)
{
	long ORFMaxLen = 0;
	vector<pair<long,long> > ORFs;
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
			ORFs.push_back(make_pair(ORF_start,length));
			//start from the next nucleotide
			ORF_starts[codon_start % 3] = index + 1;
		}
	}
	return boost::make_tuple(ORFMaxLen,ORFs);
}

