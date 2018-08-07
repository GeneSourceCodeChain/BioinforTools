#include <string>
#include "Model.hpp"

#define BATCHSIZE 51
#define TRAINSIZE 10000

using namespace std;

class Pathogenicity : public Model {
	string setup_sa_model(string input = "sa_sequence",bool isTrainable);
	string setup_ss_model(string input = "ss_sequence",bool isTrainable);
	string residual_unit(string input,string output,int in_size,int out_size,int kernel,bool isResidual = true,bool isTrainable = true);
protected:
	virtual void setupMainNet(ModelUtil & network,bool isTrainable = true);
public:
	Pathogenicity(const string & name);
	virtual ~Pathogenicity();
};
