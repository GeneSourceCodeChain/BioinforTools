#include <stdexcept>
#include "Model.hpp"

using namespace std;

Model::Model(const string & name)
:model_name(name),workspace(nullptr)
{
	LoadModule("","libcaffe2_detectron_ops_gpu.so");
	//setup structure
	ModelUtil network(init,train);
	setupMainNet(network,true);
	setupSaveNet();
	setupLoadNet();
	//choose device to use
	init.mutable_device_option()->set_device_type(CUDA);
	train.mutable_device_option()->set_device_type(CUDA);
	save.mutable_device_option()->set_device_type(CUDA);
	load.mutable_device_option()->set_device_type(CUDA);
	//create train workspace
	workspace.RunNetOnce(init);
	//instantiate network
	train_net = CreateNet(train,workspace);
	save_net = CreateNet(save,workspace);
	load_net = CreateNet(load,workspace);
}

void Model::train()
{
	//do the iteration
	for(int i = 0 ; ; i++) {
		cout<<"iter:"<<i<<endl;
		train_net->Run();
		if(i % 1000 == 0) {
			cout<<"saving params"<<endl;
			remove_all(model_name + "_params");
			save_net->Run();
		}
	}
}

bool Model::exportModel(const string & dir)
{
	if(false == exists(dir) || false == is_directory(dir)) {
		cerr<<"invalid directory!"<<endl;
		return false;
	}
	NetDef init,predict;
	ModelUtil network(init,predict);
	setupMainNet(network,false);
	network.predict.WriteText((path(dir) / (model_name + "_predict.pbtxt")).string());
	load.WriteText((path(dir) / (model_name + "_load.pbtxt")).string());
	return true;
}

void Model::setupMainNet(ModelUtil & network,bool isTrainable)
{
}

void Model::setupSaveNet()
{
	NetUtil InitNet(init);
	NetUtil SaveNet(save);
	vector<string> params;
	for(auto & op : InitNet.net.op()) {
		if(op.type() == "CreateDB") continue;
		for(auto & output: op.output())
			params.push_back(output);
	}
	SaveNet.AddSaveOp(params,"lmdb",model_name + "_params");
}

void Model::setupLoadNet()
{
	NetUtil InitNet(init);
	NetUtil LoadNet(load);
	vector<string> params;
	for(auto & op : InitNet.net.op()) {
		if(op.type() == "CreateDB") continue;
		for(auto & output : op.output())
			params.push_back(output);
	}
	LoadNet.AddLoadOp(params,"lmdb",model_name + "_params");
}

