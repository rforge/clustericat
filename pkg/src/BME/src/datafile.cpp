/*
 * datafile.cpp
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#include "datafile.h"
//#include <vector>
#include "../../eigen/Eigen/Dense"
using namespace Eigen;
//using namespace Rcpp;



datafile::datafile(MatrixXi  mat){
	// TODO Auto-generated constructor stub
	m_n=mat.rows();
	m_d=mat.cols();
	m_localise.resize(m_n);
	MatrixXi mat2;
	VectorXd eff;
	mat2.resize(m_n,m_d);
	mat2.Zero(m_n,m_d);
	mat2.row(0)=mat.row(0);
	eff.resize(m_n);
	eff.Zero(m_n);
	eff(0)=1;
	m_localise(0)=0;
	m_unique=1;
	int vu=0;
	for (int i=1;i<m_n;i++){
		vu=0;
		for (int j=0;j<(m_unique);j++){
			if ((mat2.row(j)==mat.row(i))){
				vu=1;
				eff(j)+=1;
				m_localise(i)=j;
				break;
			}
		}
		if (vu==0){
			mat2.row(m_unique)=mat.row(i);
			eff(m_unique)=1;
			m_localise(i)=m_unique;
			m_unique+=1;
		}
	}
	m_mat.resize(m_unique,m_d);
	m_mat=mat2.block(0,0,m_unique,m_d);
	m_effectif.resize(m_unique);
	m_effectif=eff.head(m_unique);
	m_modalities.resize(m_d);
	m_modalities=m_mat.colwise().maxCoeff();
}


datafile::datafile(MatrixXi  mat, VectorXi moda){
	// TODO Auto-generated constructor stub
	m_n=mat.rows();
	m_d=mat.cols();
	m_localise.resize(m_n);
	MatrixXi mat2;
	VectorXd eff;
	mat2.resize(m_n,m_d);
	mat2.Zero(m_n,m_d);
	mat2.row(0)=mat.row(0);
	eff.resize(m_n);
	eff.Zero(m_n);
	eff(0)=1;
	m_localise(0)=0;
	m_unique=1;
	int vu=0;
	for (int i=1;i<m_n;i++){
		vu=0;
		for (int j=0;j<(m_unique);j++){
			if ((mat2.row(j)==mat.row(i))){
				vu=1;
				eff(j)+=1;
				m_localise(i)=j;
				break;
			}
		}
		if (vu==0){
			mat2.row(m_unique)=mat.row(i);
			eff(m_unique)=1;
			m_localise(i)=m_unique;
			m_unique+=1;
		}
	}
	m_mat.resize(m_unique,m_d);
	m_mat=mat2.block(0,0,m_unique,m_d);
	m_effectif.resize(m_unique);
	m_effectif=eff.head(m_unique);
	m_modalities.resize(m_d);
	m_modalities=moda;
}


datafile::~datafile() {
	// TODO Auto-generated destructor stub
}

const MatrixXi & datafile::Get_mat_datafile(){
	return m_mat;
}

const VectorXd & datafile::Get_eff_datafile(){
	return m_effectif;
}

VectorXd  datafile::Get_eff_datafile(VectorXd poids){
	VectorXd eff=m_effectif;
	eff=eff.array()*poids.array();
	return eff;
}

const VectorXi  datafile::Get_modalite(VectorXi who){
	VectorXi moda;
	moda.resize(who.rows());
	for (int j=0;j<who.rows();j++){
		moda(j)=m_modalities(who(j));
	}
	return(moda);
}

const VectorXi & datafile::Get_localise(){
	return m_localise;
}
