#ifndef MODEL_HPP
#define MODEL_HPP

#include <cstdlib>
#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <caffe2/core/init.h>
#include <caffe2/core/context_gpu.h>
#include <caffe2/core/operator.h>
#include <caffe2/core/module.h>
#include <caffe2/util/blob.h>
#include <caffe2/util/model.h>
#include <caffe2/util/net.h>
#include <cvplot/cvplot.h>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace boost::filesystem;
using namespace cv;
using namespace caffe2;
using namespace cvplot;

class Model {
protected:
	const string model_name;
	//network structure
	NetDef init;
	NetDef train;
	NetDef save;
	NetDef load;
	//network instance
	unique_ptr<NetBase> train_net;
	unique_ptr<NetBase> save_net;
	unique_ptr<NetBase> load_net;
	//workspace
	Workspace workspace;
	//structure setup function
	virtual void setupMainNet(ModelUtil & network,bool isTrainable = true) = 0;
	virtual void setupSaveNet();
	virtual void setupLoadNet();
public:
	Model(const string & name);
	virtual void train();
	virtual bool exportModel(const string & dir);
};

#endif

