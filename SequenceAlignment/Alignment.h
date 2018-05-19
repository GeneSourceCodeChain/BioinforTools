#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <string>
#include <map>
#include <functional>
#include <vector>
#include <stack>
#include <utility>
#include <limits>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
namespace ublas = boost::numeric::ublas;

class Alignment {
public:
	Alignment();
	virtual ~Alignment();
	template<typename T> vector<pair<char,char> > alignment(string a,string b,T score_policy);
	void print_alignment(const vector<pair<char,char> > & a);
};

template<typename T> 
vector<pair<char,char> > Alignment::alignment(string a, string b, T policy)
{
	boost::to_upper(a);
	boost::to_upper(b);
	//align with Smith-Waterman algorithm
	//1)create score table
	ublas::matrix<pair<int,int> > tracetable(a.size() + 1,b.size() + 1,make_pair(0,0));
	for(int h = 1 ; h < tracetable.size1() ; h++) tracetable(h,0) = make_pair(h - 1,0);
	for(int w = 1 ; w < tracetable.size2() ; w++) tracetable(0,w) = make_pair(0,w - 1);
	int maxscore = numeric_limits<int>::min();
	pair<int,int> maxpos;
	{
		//to free score matrix right after its usage
		ublas::matrix<int> score = ublas::zero_matrix<int>(a.size() + 1,b.size() + 1);
		for(int h = 1 ; h < score.size1() ; h++)
			for(int w = 1 ; w < score.size2() ; w++) {
				int AForwardBForward = score(h - 1,w - 1) + policy(a[h - 1],b[w - 1]);		//both string A and B consume one alphabet
				int AForwardBStop = score(h - 1, w) + policy(a[h - 1],'-');			//string A consumes one alphabet
				int AStopBForward = score(h, w - 1) + policy('-',b[w - 1]);			//string B consumes one alphabet
				pair<int,int> step(h,w);
				int s = (AForwardBForward > 0)?(step = make_pair(h - 1,w - 1),AForwardBForward):0;
				s = (AForwardBStop > s)?(step = make_pair(h - 1, w),AForwardBStop):s;
				s = (AStopBForward > s)?(step = make_pair(h, w - 1),AStopBForward):s;
				score(h,w) = s;									//save score for latter calculation
				tracetable(h,w) = step;							//save trace to trace back easily
				maxscore = (s > maxscore)?(maxpos = make_pair(h,w),s):maxscore;
			}
	}
	//2)back trace
	stack<pair<int,int> > backtrace;
	pair<int,int> pos = maxpos;
	do {
		backtrace.push(pos);
		pos = tracetable(pos.first,pos.second);
	} while(false == (pos.first == backtrace.top().first && pos.second == backtrace.top().second));
	//3) convert trace into readable format
	vector<pair<char,char> > trace;
	pair<int,int> prev_step(-1,-1);
	while(backtrace.size()) {
		pair<int,int> step = backtrace.top();
		backtrace.pop();
		if(prev_step.first == -1 || prev_step.second == -1) {
			//without gap
			trace.push_back(make_pair(a[step.first - 1],b[step.second - 1]));
		} else {
			//with gap
#ifndef NDEBUG
			assert(false == (step.first == prev_step.first && step.second == prev_step.second));
#endif
			if(step.first == prev_step.first) {
#ifndef NDEBUG
				assert(0 <= step.second - 1 && step.second - 1 < b.size());
#endif
				trace.push_back(make_pair('-',b[step.second - 1]));
			} else if(step.second == prev_step.second) {
#ifndef NDEBUG
				assert(0 <= step.first - 1 && step.first - 1 < a.size());
#endif
				trace.push_back(make_pair(a[step.first - 1],'-'));
			} else {
#ifndef NDEBUG
				assert(0 <= step.first - 1 && step.first - 1 < a.size());
				assert(0 <= step.second - 1 && step.second - 1 < b.size());
#endif
				trace.push_back(make_pair(a[step.first - 1],b[step.second - 1]));
			}
		}
		prev_step = step;
	}
	return trace;
}

class BasicPolicy : public binary_function<char, char, int> {
public:
	int operator()(const char & c1,const char & c2) {
#ifndef NDEBUG
		assert('A' <= c1 && c1 <= 'Z' || c1 == '-');
		assert('A' <= c2 && c2 <= 'Z' || c2 == '-');
#endif
		return c1 == c2 ? 2 : -1;
	}
};

class PenaltyPolicy : public binary_function<char, char, int> {
	ublas::matrix<int> BLOSUM62;
	map<char,int> index;
	int gap_penalty;
public:
	PenaltyPolicy(int gap_penalty = -4);
	int operator()(const char & c1,const char & c2) {
#ifndef NDEBUG
		assert('A' <= c1 && c1 <= 'Z' || c1 == '-');
		assert('A' <= c2 && c2 <= 'Z' || c2 == '-');
#endif
		if(c1 == '-' || c2 == '-') return gap_penalty;
		int row = index[c1]; assert(row != -1); //can't input unrealistic protein
		int col = index[c2]; assert(col != -1); //can't input unrealistic protein
		return BLOSUM62(row,col);
	}
};

#endif

