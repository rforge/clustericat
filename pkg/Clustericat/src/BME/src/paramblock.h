/*
 * paramblock.h
 *
 *  Created on: 11 juil. 2012
 *      Author: matthieu
 */

#ifndef PARAMBLOCK_H_
#define PARAMBLOCK_H_
#include "../../eigen/Eigen/Dense"
#include <vector>
#include "datafile.h"
#include "model.h"
using namespace std;
using namespace Eigen;

class param_block {
protected:
	vector<VectorXd> m_alpha;
	VectorXd m_beta;
	VectorXi m_ind_st;
	MatrixXi m_gamma;
	double m_epsilon;

public:
	param_block();
//	param_block(int , int , datafile , model);
	param_block(int , int , datafile , model, VectorXd ,int);
	param_block(int , int , datafile, model ,MatrixXi);
	void Genere_gamma2(VectorXi modalite);
	void Genere_gamma(int , int , datafile , model,VectorXd);
	void Genere_gamma3(int , int , datafile , model,VectorXd);
	double Proba_bloc_gamma(const MatrixXi &, const VectorXd &  eff, const MatrixXi &g1);
	void Qui_dpt(const VectorXi & who, const MatrixXi & mat, const MatrixXi & gamma);
	virtual ~param_block();
	VectorXd uniform(int dim);
	VectorXd proba_indpt(const VectorXi & who, const MatrixXi & mat);
	VectorXd proba_dpt(const VectorXi & who,const MatrixXi & mat);
	param_block Optimise_gamma(int k, int b, datafile dat, model mod, int iter,const VectorXd & poids,const VectorXd & eff);
	double Estime_theta_cond(int k, int b, datafile dat, model mod,const VectorXd & poids,const MatrixXi & gam,const VectorXd & eff);
	param_block actu_param(const VectorXi & who,const MatrixXi & mat,const VectorXd & tik,const VectorXd & poids,const VectorXd & eff,param_block);
	VectorXd proba_dpt(const VectorXi & who,const MatrixXi & mat,param_block parambl);
	VectorXd proba_indpt(const VectorXi & who, const MatrixXi & mat,param_block parambl);
	void Mstep(const VectorXi & , const MatrixXi &, const VectorXd &, const VectorXd &, const VectorXd &);
	MatrixXi Genere_possible_gam(int,int);
	int isNaN();
	vector<VectorXd> Get_alpha();
	VectorXd Get_tau();
	double Get_rho();
	MatrixXi Get_delta();
};

#endif /* PARAMBLOCK_H_ */
