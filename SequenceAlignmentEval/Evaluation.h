#ifndef EVALUATION_H
#define EVALUATION_H

#include <string>
#include <set>
#include <vector>
#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

using namespace std;
using namespace boost::multi_index;

struct ORF {
	long start_pos;
	long length;
	ORF(long s,long l):start_pos(s),length(l){}
};

typedef multi_index_container<
	ORF,
	indexed_by<
		ordered_unique<
			member<ORF,long,&ORF::start_pos>
		>,
		ordered_non_unique<
			member<ORF,long,&ORF::length>
		>
	>
> ORFList;

class Evaluation {
	set<long> stopCodonPos;
public:
	Evaluation(string genbank);
    	virtual ~Evaluation();
	bool evaluate(string fasta,string output);
protected:
	set<long> getStopCodonPos(string genbank);
	string getDNASequence(string fasta);
	boost::tuple<long,ORFList> getORF(string sequence);
	template<typename T> boost::tuple<
		map<char,double>,
		map<char,map<char,double> >,
		map<char,map<char,map<char,double> > >,
		map<char,map<char,map<char,map<char,double> > > >
	> count(const string & sequence,const ORFList & ORFs,const T & predicate);
	double logsum(double lx,double ly);
};

template<typename T> 
boost::tuple<
	map<char,double>,
	map<char,map<char,double> >,
	map<char,map<char,map<char,double> > >,
	map<char,map<char,map<char,map<char,double> > > > 
> Evaluation::count(const string & sequence, const ORFList & ORFs, const T & predicate)
{
	long singlesum = 0;
	long doublesum = 0;
	long triplesum = 0;
	long quadsum = 0;
	map<char,double> singlecount;
	map<char,map<char,double> > doublecount;
	map<char,map<char,map<char,double> > > triplecount;
	map<char,map<char,map<char,map<char,double> > > > quadcount;
	for(auto & orf : ORFs) {
		if(predicate(orf.length)) {
			//count how many singles doubles triples quad can be contained in the sequence
			singlesum += orf.length;
			doublesum += orf.length - 1;
			triplesum += orf.length - 2;
			quadsum += orf.length - 3;
			//count the occurrence of each type of nucleotide
			for(int index = orf.start_pos ; index < orf.start_pos + orf.length ; index++) singlecount[sequence[index]]+=1;
			for(int index = orf.start_pos ; index < orf.start_pos + orf.length - 1 ; index++) doublecount[sequence[index]][sequence[index + 1]]+=1;
			for(int index = orf.start_pos ; index < orf.start_pos + orf.length - 2 ; index++) triplecount[sequence[index]][sequence[index + 1]][sequence[index + 2]]+=1;
			for(int index = orf.start_pos ; index < orf.start_pos + orf.length - 3 ; index++) quadcount[sequence[index]][sequence[index + 1]][sequence[index + 2]][sequence[index + 3]]+=1;
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
	return boost::make_tuple(singlecount,doublecount,triplecount,quadcount);
}

class IsLongPeptideChain : public unary_function<long,bool> {
public:
	bool operator()(const long & length) const {
		return length >= 1400;
	}
};

class IsShortPeptideChain : public unary_function<long,bool> {
public:
	bool operator()(const long & length) const {
		return length <= 50;
	}
};

#endif

