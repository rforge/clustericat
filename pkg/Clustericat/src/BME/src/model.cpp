/*
 * model.cpp
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#include <iostream>
#include "model.h"
#include "../../eigen/Eigen/Dense"
using namespace std;
using namespace Eigen;

model::model(){

}



model::model(VectorXi mo, int nbc) {
	// TODO Auto-generated constructor stub
	m_g=nbc;
	m_omega.resize(m_g,mo.rows());
	for (int k=0;k<m_g;k++){
		m_omega.row(k)=mo;
	}
	m_block.resize(m_g);
	for (int k=0;k<m_g;k++){
		m_block[k].resize(m_omega.rowwise().maxCoeff()(k)+1);
		for (int b = 0; b < (m_omega.rowwise().maxCoeff()(k)+1); ++b) {
			if ((m_omega.row(k).array()==b).count()>0){
				m_block[k][b].resize((m_omega.row(k).array()==b).count());
				int pl=0,var=0;
				while(pl<((m_omega.row(k).array()==b).count())){
					if (m_omega(k,var)==b){
						m_block[k][b](pl)=var;
						pl++;
					}
					var++;
				}
			}
		}
	}

}

model::model(MatrixXi mo, int nbc) {
	// TODO Auto-generated constructor stub
	m_g=nbc;
	m_omega=mo;
	m_block.resize(m_g);
	for (int k=0;k<m_g;k++){
		m_block[k].resize(m_omega.rowwise().maxCoeff()(k)+1);
		for (int b = 0; b < (m_omega.rowwise().maxCoeff()(k)+1); ++b) {
			if ((m_omega.row(k).array()==b).count()>0){
				m_block[k][b].resize((m_omega.row(k).array()==b).count());
				int pl=0,var=0;
				while(pl<((m_omega.row(k).array()==b).count())){
					if (m_omega(k,var)==b){
						m_block[k][b](pl)=var;
						pl++;
					}
					var++;
				}
			}
		}
	}

}
model::~model() {
	// TODO Auto-generated destructor stub
}

const MatrixXi & model::Get_model(){
	return(m_omega);
}

const VectorXi & model::Get_var_block(int k,int b){
	return(m_block[k][b]);
}


int model::Get_nb_param(){
	return(nb_param);
}

int model::Get_nbcl(){
	return(m_g);
}
