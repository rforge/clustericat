/*
 * datafile.h
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#ifndef DATAFILE_H_
#define DATAFILE_H_
#include "../../eigen/Eigen/Dense"
using namespace std;
using namespace Eigen;

class datafile {
protected:
	int m_n,m_d,m_unique;
	MatrixXi m_mat;
	VectorXd m_effectif ;
	VectorXi m_modalities,m_localise;

public:
	datafile(MatrixXi data);
	datafile(MatrixXi data, VectorXi moda);
	virtual ~datafile();
	const MatrixXi & Get_mat_datafile();
	const VectorXd & Get_eff_datafile();
	VectorXd  Get_eff_datafile(VectorXd);
	const VectorXi Get_modalite(VectorXi);
	const VectorXi & Get_localise();
};
#endif /* DATAFILE_H_ */
