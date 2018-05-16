#ifndef MARKOVMODEL_H
#define MARKOVMODEL_H

#include <string>
#include <set>
#include <vector>
#include <functional>
#include <boost/tuple/tuple.hpp>

using namespace std;

class MarkovModel {
	set<long> stopCodonPos;
	string sequence;
	long ORFMaxLen;
	vector<pair<long,long> > ORFs;
public:
    MarkovModel(string genbank,string fasta);
    virtual ~MarkovModel();
protected:
	set<long> getStopCodonPos(string genbank);
	string getDNASequence(string fasta);
	boost::tuple<long,vector<pair<long,long> > > getORF(string sequence);
	auto count(unary_function<long,bool> & predicate);
};

#endif

