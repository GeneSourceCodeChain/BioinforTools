#include <boost/lexical_cast.hpp>
#include "Pathogenicity.hpp"

using namespace boost;

Pathogenicity::Pathogenicity(const string & name)
:Model(name)
{
}

Pathogenicity::~Pathogenicity()
{
}

void Pathogenicity::setupMainNet(ModelUtil & network,bool isTrainable)
{
	network.init.AddCreateDbOp("db","lmdb","dataset");
	network.predict.AddInput("db");
	//dimensions of input tensors are all 51(batch)x1(channel)x20(w)
	network.AddTensorProtosDbInputOp(
		"db",{
			"orig_seq","snp_seq","conservation_full",
			"conservation_primates","conservation_otherspecies",
			"label"
		},
		BATCHSIZE
	);
	network.AddSumOp({"conservation_full","conservation_primates","conservation_otherspecies"},"cons_counts");
	network.AddStopGradientOp("cons_counts");
	string struc = setup_ss_model("cons_counts",isTrainable);
	string solv = setup_sa_model("cons_counts",isTrainable);
	network.AddConv1DOps("orig_seq","1d_conv_orig",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("snp_seq","1d_conv_snp",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("conservation_full","1d_conv_mammals",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("conservation_primates","1d_conv_primates",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("conservation_otherspecies","1d_conv_vertebrates",1,40,{1},{0,0},{1},!isTrainable);
	network.AddSumOp({"1d_conv_orig","1d_conv_mammals","1d_conv_primates","1d_conv_vertebrates",struc,solv},"merge_orig_conserv");
	network.AddSumOp({"1d_conv_snp","1d_conv_mammals","1d_conv_primates","1d_conv_vertebrates",struc,solv},"merge_snp_conserv");
	residual_unit("merge_orig_conserv","orig_residual",40,40,5,false,isTrainable);
	residual_unit("merge_snp_conserv","snp_residual",40,40,5,false,isTrainable);
	//snp_orig_conv_merge dimension is batch x 40(channel) x (20 + 20)(w)
	//skip_snp_orig_conv_merge dimension is batch x 40(channel) x (20 + 20)(w)
	network.AddConcatOp({"orig_residual","snp_residual"},"snp_orig_conv_merge",2);
	network.AddConcatOp({"orig_residual","snp_residual"},"skip_snp_orig_conv_merge",2);
	network.AddConv1DOps("snp_orig_conv_merge","1d_conv_reduce",40,40,{1},{2,2},{5},!isTrainable);
	network.AddConv1DOps("skip_snp_orig_conv_merge","1d_skip_reduce",40,40,{1},{2,2},{5},!isTrainable);
	string conv_output = "1d_conv_reduce";
	string skip_output = "1d_skip_reduce";
	for(int i = 0 ; i < 3 ; i++) {
		//conv branch
		for(int j = 0 ; j < 2 ; j++) {
			residual_unit(conv_output,"residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1),40,40,5,true,isTrainable);
			conv_output = "residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1);
		}
		network.AddConv1DOps(conv_output,"denseforskip_" + lexical_cast<string>(i + 1),40,40,{1},{0,0},{1},!isTrainable);
		//sum conv branch with skip branch
		network.AddSumOp({skip_output,"denseforskip_" + lexical_cast<string>(i + 1)},"skip_" + lexical_cast<string>(i + 1));
		skip_output = "skip_" + lexical_cast<string>(i + 1);
	}
	residual_unit(skip_output,"residual_final",40,40,1,true,isTrainable);
	//conv_final dimension is batch x 2 x (20 + 20) (w)
	network.AddConv1DOps("residual_final","conv_final",40,2,{1},{0,0},{1},!isTrainable);
	network.AddSigmoidOp("conv_final","conv_final");
	//output dimension is batch x 2 x 1
	network.AddGlobalMaxPool1DOp("conv_final","raw_output");
	network.AddReshapeOp("raw_output","output",{BATCHSIZE,2});
	//classification
	network.AddSoftmaxOp("output","softmax",1);
	if(isTrainable) {
		//loss
		network.AddLabelCrossEntropyOp("softmax","label","xent");
		network.AddAveragedLossOp("xent","loss");
		network.AddAccuracyOp("softmax","label","accuracy");
		//backpropagation
		network.AddConstantFillWithOp(1.0,"loss","loss_grad");
		network.predict.AddGradientOps();
		network.AddIterOps();
		network.AddLearningRateOp("iter","lr",-0.01,0.9,100*round(static_cast<float>(TRAINSIZE)/BATCHSIZE));
		string optimizer = "adam";
		network.AddOptimizerOps(optimizer);
	}
}

string Pathogenicity::setup_sa_model(string input,bool isTrainable)
{
	//create network;
	ModelUtil network(init,predict);
	//conv branch (batchsize x 40 x 20)
	network.AddConv1DOps(input,"sa_1d_conv_sequence",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("sa_1d_conv_sequence","sa_1d_conv_down",40,40,{1},{0,0},{1},!isTrainable);
	//skip branch (batchsize x 40 x 20)
	network.AddConv1DOps(input,"sa_skip_1d_conv_sequence",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("sa_skip_1d_conv_sequence","sa_1d_skip_down",40,40,{1},{0,0},{1},!isTrainable);
	//setup residual units
	string conv_output = "sa_1d_conv_down";
	string skip_output = "sa_1d_skip_down";
	for(int i = 0 ; i < 3 ; i++) {
		//conv branch
		for(int j = 0 ; j < 2 ; j++) {
			residual_unit(conv_output,"sa_residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1),40,40,5,true,isTrainable);
			conv_output = "sa_residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1);
		}
		network.AddConv1DOps(conv_output,"sa_denseforskip_" + lexical_cast<string>(i + 1),40,40,{1},{0,0},{1},!isTrainable);
		//sum conv branch with skip branch
		network.AddSumOp({skip_output,"sa_denseforskip_" + lexical_cast<string>(i + 1)},"sa_skip_" + lexical_cast<string>(i + 1));
		//update skip branch to the latest output tensor
		skip_output = "sa_denseforskip_" + lexical_cast<string>(i + 1);
	}
	residual_unit(skip_output,"sa_residual_final",40,40,1,true,isTrainable);
	return "sa_residual_final";
}

string Pathogenicity::setup_ss_model(string input,bool isTrainable)
{
	//create network
	ModelUtil network(init,predict);
	//conv branch (batchsize x 40 x 20)
	network.AddConv1DOps(input,"ss_1d_conv_sequence",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("ss_1d_conv_sequence","ss_1d_conv_down",40,40,{1},{0,0},{1},!isTrainable);
	//skip branch (batchsize x 40 x 20)
	network.AddConv1DOps(input,"ss_skip_1d_conv_sequence",1,40,{1},{0,0},{1},!isTrainable);
	network.AddConv1DOps("ss_skip_1d_conv_sequence","ss_1d_skip_down",40,40,{1},{0,0},{1},!isTrainable);
	//setup residual units;
	string conv_output = "ss_1d_conv_down";
	string skip_output = "ss_1d_skip_down";
	for(int i = 0 ; i < 3 ; i++) {
		//conv branch
		for(int j = 0 ; j < 2 ; j++) {
			residual_unit(conv_output,"ss_residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1),40,40,5,true,isTrainable);
			conv_output = "ss_residual_" + lexical_cast<string>(i + 1) + "_" + lexical_cast<string>(j + 1);
		}
		network.AddConv1DOps(conv_output,"ss_denseforskip_" + lexical_cast<string>(i + 1),40,40,{1},{0,0},{1},!isTrainable);
		//sum conv branch with skip branch
		network.AddSumOp({skip_output,"ss_denseforskip_" + lexical_cast<string>(i + 1)},"ss_skip_" + lexical_cast<string>(i + 1));
		//update skip branch to the latest output tensor
		skip_output = "ss_denseforskip_" + lexical_cast<string>(i + 1);
	}
	residual_unit(skip_output,"ss_residual_final",40,40,1,true,isTrainable);
	return "ss_residual_final";
}

string Pathogenicity::residual_unit(string input,string output,int in_size,int out_size,int kernel,bool isResidual, bool isTrainable) {
	//input dimension (batchsize x channel x w)
	//output dimension (batchsize x channel x w)
	ModelUtil network(init,predict);
	network.AddSpatialBNOps(input,"sa_BatchNormalization_" + output + "_1",in_size,1e-5f,0.9,!isTrainable);
	network.AddReluOp("sa_BatchNormalization_" + output + "_1","sa_relu_" + output + "_1");
	network.AddConv1DOps("sa_relu_" + output + "_1","sa_conv1d_" + output + "_1",in_size,out_size,{1},{kernel / 2,kernel / 2},{kernel},!isTrainable);
	network.AddSpatialBNOps("sa_conv1d_" + output + "_1","sa_BatchNormalization_" + output + "_2",out_size,1e-5f,0.9,!isTrainable);
	network.AddReluOp("sa_BatchNormalization_" + output + "_2","sa_relu_" + output + "_2");
	network.AddConv1DOps("sa_relu_" + output + "_2","sa_conv1d_" + output + "_2",out_size,out_size,{1},{kernel / 2,kernel / 2},{kernel},!isTrainable);
	if(isResidual)
		network.AddSumOp({"sa_conv1d_" + output + "_2",input},output);
	else
		network.AddCopyOp("sa_conv1d_" + output + "_2",output);
	return output;
}

