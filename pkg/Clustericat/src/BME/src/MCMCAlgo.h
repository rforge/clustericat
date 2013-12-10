/*
 * MCMCAlgo.h
 *
 *  Created on: 12 juil. 2012
 *      Author: matthieu
 */

#ifndef MCMCALGO_H_
#define MCMCALGO_H_
#include "../../eigen/Eigen/Dense"
#include "datafile.h"
#include "model.h"
#include "parameters.h"
#include <vector>
//#include "paramblock.h"

class MCMCAlgo {
protected:
	model best_model;
	parameters best_parametres;

public:
	MCMCAlgo(datafile dat,MatrixXi mo, int nbc,int nbinit,int nbiter,int iterdiscret,int itercond,int iterglobal);
	MCMCAlgo();
	virtual ~MCMCAlgo();
	double Get_best_bic();
	double Get_best_like();
	MatrixXi Get_best_omega();
	MatrixXd Sauv_probapost();
	int Get_nb_param();
	vector<VectorXd> Get_alpha(int,int);
	VectorXd Get_tau(int,int);
	double Get_rho(int,int);
	MatrixXi Get_delta(int,int);
	VectorXd Get_propor();
};

#endif /* MCMCALGO_H_ */
