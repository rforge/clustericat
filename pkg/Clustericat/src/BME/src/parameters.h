/*
 * parameters.h
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include "../../eigen/Eigen/Dense"
#include "paramblock.h"
#include <vector>
using namespace std;

class parameters {
protected:
	vector<vector< param_block > > m_param;
	VectorXd m_propor;
	vector<MatrixXd> m_proba_block;
	MatrixXd m_proba;
	double m_loglikelihood,m_bic;
	int m_nbparam;

public:
	parameters();
	parameters(datafile, model);
	parameters(datafile, model,model,parameters,int);
	parameters(datafile, model,model,parameters,int,int);
	virtual ~parameters();
	VectorXd uniforme(int dim);
	void Optimise_gamma(int,int,datafile, model ,int);
	void Probapost(model,const MatrixXi & );
	void Likelihood(const VectorXd & eff);
	void Compte_nbparam(datafile,model);
	void Estimation(int,int,int, datafile , model );
	void Estimation_discrete(int , datafile , model);
	void Estimation_continue(int, datafile, model);
	void Affiche_param();
	void Mstep(datafile , model);
	double Get_bic();
	double Get_like();
	MatrixXd Sauve_proba_post();
	int Get_nb_param();
	vector<VectorXd> Get_alpha(int,int);
	VectorXd Get_tau(int,int);
	double Get_rho(int,int);
	MatrixXi Get_delta(int,int);
	VectorXd Get_propor();

};

#endif /* PARAMETERS_H_ */
