/*
 * model.h
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#ifndef MODEL_H_
#define MODEL_H_
#include "../../eigen/Eigen/Dense"
#include <vector>
using namespace std;
using namespace Eigen;

class model {
protected:
	MatrixXi m_omega;
	vector<vector<VectorXi> > m_block;
	int m_g,nb_param;


public:
	model();
	model(VectorXi mo, int nbc);
	model(MatrixXi mo, int nbc);
	virtual ~model();
	const MatrixXi & Get_model();
	const VectorXi & Get_var_block(int,int);
	int Get_nbcl();
	int Get_nb_param();

};

#endif /* MODEL_H_ */
