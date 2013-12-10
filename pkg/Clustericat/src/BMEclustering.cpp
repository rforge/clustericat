/*
 * rcpp_reading.cpp
 *
 *  Created on: 2 oct. 2012
 *      Author: matthieu
 */

#include "BMEclustering.h"
#include "conversion.h"
#include "eigen/Eigen/Dense"
#include "BME/src/datafile.h"
#include "BME/src/datafile.cpp"
#include "BME/src/model.h"
#include "BME/src/model.cpp"
#include "BME/src/MCMCAlgo.h"
#include "BME/src/MCMCAlgo.cpp"
#include "BME/src/paramblock.h"
#include "BME/src/parameters.h"
#include "BME/src/paramblock.cpp"
#include "BME/src/parameters.cpp"
#include <vector>
#include <time.h>

using namespace Rcpp ;
using namespace Eigen ;

RcppExport SEXP BMEclustering(SEXP mat, SEXP moda, SEXP nb_cluster, SEXP partition_initiale, SEXP nb_init,SEXP stop_criterion){
	srand(time(0));

	MatrixXi data=convertMatrix<MatrixXi,NumericMatrix>(mat);
	VectorXi modalite=convertvector<VectorXi,NumericVector>(moda);
	datafile dataF(data,modalite);

	NumericMatrix red=convertMatrix<NumericMatrix,MatrixXi>(dataF.Get_mat_datafile());
	VectorXi partition_vbles=convertvector<VectorXi,NumericVector>(partition_initiale);
	MatrixXi m;
	int g=as<int>(nb_cluster);
	//int borne=as<int>(nbiter);
	m.resize(g,partition_vbles.rows());
	for (int k=0;k<g;k++){
		m.row(k)=partition_vbles;
	}
	NumericMatrix test=convertMatrix<NumericMatrix,MatrixXi>(m);
	int nbinit=as<int>(nb_init);
	int borne=as<int>(stop_criterion);
	MCMCAlgo ref(dataF,m,m.rows(),2,1,2,6,borne);
	for (int ini=0;ini<nbinit;ini++){
		MCMCAlgo test(dataF,m,m.rows(),2,1,2,6,borne);
		if (test.Get_best_bic()>ref.Get_best_bic()){
			ref=test;
		}
	}
	//sauvegarde des caractéristiques du modèle
	NumericMatrix model=convertMatrix<NumericMatrix,MatrixXi>(ref.Get_best_omega());
	double bic=ref.Get_best_bic();
	double likelihood=ref.Get_best_like();
	NumericMatrix probapost=convertMatrix<NumericMatrix,MatrixXd>(ref.Sauv_probapost());
	NumericVector localise=convertvector<NumericVector,VectorXi>(dataF.Get_localise());

	vector< vector< NumericVector > > tau;
	vector< vector< vector< NumericVector > > > delta;
	vector< vector< vector< NumericVector > > > alpha;
	vector< vector< double > > rho;
	tau.resize(g);
	rho.resize(g);
	delta.resize(g);
	alpha.resize(g);
	for (int k=0;k<g;k++){
		tau[k].resize(ref.Get_best_omega().row(k).maxCoeff()+1);
		rho[k].resize(ref.Get_best_omega().row(k).maxCoeff()+1);
		delta[k].resize(ref.Get_best_omega().row(k).maxCoeff()+1);
		alpha[k].resize(ref.Get_best_omega().row(k).maxCoeff()+1);
		for (int b=0;b<=ref.Get_best_omega().row(k).maxCoeff();b++){
		//for (int b=0;b<2;b++){
			//vectblock[k][b]=ref.Get_rho(k,b);
			//vector < VectorXd > passe=ref.Get_alpha(k,b);
			if ((ref.Get_best_omega().row(k).array()==b).count()>0){
				vector<VectorXd> passe=ref.Get_alpha(k,b);
				alpha[k][b].resize((ref.Get_best_omega().row(k).array()==b).count());
				for (int loc=0;loc<((ref.Get_best_omega().row(k).array()==b).count());loc++){
					alpha[k][b][loc]=convertvector<NumericVector,VectorXd>(passe[loc]);
				}
			}
			if ((ref.Get_best_omega().row(k).array()==b).count()>1){
				MatrixXi passe=ref.Get_delta(k,b);
				delta[k][b].resize(passe.cols());
				tau[k][b]=convertvector<NumericVector,VectorXd>(ref.Get_tau(k,b));
				for (int loc=0;loc<passe.cols();loc++){
					delta[k][b][loc]=convertvector<NumericVector,VectorXi>(passe.col(loc));
				}
				//delta[k][b]=convertMatrix<NumericMatrix,MatrixXi>(ref.Get_delta(k,b));
				rho[k][b]=ref.Get_rho(k,b);
			}
		}
	}

	List param=List::create(Rcpp::Named("tau")=tau,Rcpp::Named("rho")=rho,Rcpp::Named("delta")=delta,Rcpp::Named("alpha")=alpha,Rcpp::Named("proportions")=convertvector<NumericVector,VectorXd>(ref.Get_propor()));
    List desc_model = List::create(Rcpp::Named("sigma")=model,Rcpp::Named("bic")=bic,Rcpp::Named("likelihood")=likelihood,Rcpp::Named("probapost")=probapost,Rcpp::Named("partition")=localise,Rcpp::Named("nbcluster")=nb_cluster,Rcpp::Named("parameters")=param);

    return desc_model;
}



