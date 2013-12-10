/*
 * parameters.cpp
 *
 *  Created on: 10 juil. 2012
 *      Author: matthieu
 */

#include "parameters.h"
#include "datafile.h"
#include "paramblock.h"
#include "model.h"
#include "../../eigen/Eigen/Dense"
#include <vector>
#include <cmath>
#include "math.h"



parameters::parameters(){

}

parameters::parameters(datafile dat, model nv_mod, model ref_mod,parameters ref_param,int compo){
	const MatrixXi & omega=nv_mod.Get_model(),ref_omega=ref_mod.Get_model(),mat=dat.Get_mat_datafile();
	const int g=omega.rows(),unique=mat.rows();
	m_proba=ref_param.m_proba;
	m_proba_block.resize(g);
	m_param.resize(g);
	for (int k=0;k<g;k++){
		if (k!=compo){
			m_param[k].resize(omega.rowwise().maxCoeff()(k)+1);
			m_param[k]=ref_param.m_param[k];
			m_proba_block[k].resize(unique,omega.rowwise().maxCoeff()(k)+1);
			m_proba_block[k]=ref_param.m_proba_block[k];
			for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
				if ((omega.row(k).array()==b).any()){
					m_param[k][b]=ref_param.m_param[k][b];
				}
			}
		}else{
			m_param[k].resize(omega.rowwise().maxCoeff()(k)+1);
			m_proba_block[k].resize(unique,omega.rowwise().maxCoeff()(k)+1);
			for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
				if ((omega.row(k).array()==b).any()){
					if ((((omega.row(k).array()==b)==(ref_omega.row(k).array()==b)).prod())==1){
						m_param[k][b]=ref_param.m_param[k][b];
					}else{
						m_param[k][b]=param_block(k,b,dat,nv_mod,m_proba.col(k).array()/m_proba.rowwise().sum().array(),1);
						if ((omega.row(k).array()==b).count()>1){
							int prem=0;
							while(omega(k,prem)!=b){prem++;}
							if (mat.col(prem).maxCoeff()>5){
									m_param[k][b]=m_param[k][b].Optimise_gamma(k,b,dat,nv_mod,5,m_proba.col(k).array()/m_proba.rowwise().sum().array(),dat.Get_eff_datafile());
							}
						}
					}
				}
			}
		}

	}
	m_propor=uniforme(g);
	Probapost( nv_mod , mat );
	Compte_nbparam(dat,nv_mod);
	Likelihood(dat.Get_eff_datafile());
	Estimation(1,0,6,dat,nv_mod);
}


parameters::parameters(datafile dat, model nv_mod, model ref_mod,parameters ref_param,int compo,int iter){
	const MatrixXi & omega=nv_mod.Get_model(),ref_omega=ref_mod.Get_model(),mat=dat.Get_mat_datafile();
	const int g=omega.rows(),unique=mat.rows();
	m_proba=ref_param.m_proba;
	m_proba_block.resize(g);
	m_param.resize(g);
	for (int k=0;k<g;k++){
		if (k!=compo){
			m_param[k].resize(omega.rowwise().maxCoeff()(k)+1);
			m_param[k]=ref_param.m_param[k];
			m_proba_block[k].resize(unique,omega.rowwise().maxCoeff()(k)+1);
			m_proba_block[k]=ref_param.m_proba_block[k];
			for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
				if ((omega.row(k).array()==b).any()){
					m_param[k][b]=ref_param.m_param[k][b];
				}
			}
		}else{
			m_param[k].resize(omega.rowwise().maxCoeff()(k)+1);
			m_proba_block[k].resize(unique,omega.rowwise().maxCoeff()(k)+1);
			for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
				if ((omega.row(k).array()==b).any()){
					if ((((omega.row(k).array()==b)==(ref_omega.row(k).array()==b)).prod())==1){
						m_param[k][b]=ref_param.m_param[k][b];
					}else{
						m_param[k][b]=param_block(k,b,dat,nv_mod,m_proba.col(k).array()/m_proba.rowwise().sum().array(),1);
						if ((omega.row(k).array()==b).count()>1){
							int prem=0;
							while(omega(k,prem)!=b){prem++;}
							if (mat.col(prem).maxCoeff()>5){
									m_param[k][b]=m_param[k][b].Optimise_gamma(k,b,dat,nv_mod,5,m_proba.col(k).array()/m_proba.rowwise().sum().array(),dat.Get_eff_datafile());
							}
						}
					}
				}
			}
		}

	}
	m_propor=uniforme(g);
	Probapost( nv_mod , mat );
	Compte_nbparam(dat,nv_mod);
	Likelihood(dat.Get_eff_datafile());
	Estimation(1,0,iter,dat,nv_mod);
}


parameters::parameters(datafile dat, model mod){
	// TODO Auto-generated constructor stub
	const MatrixXi & omega=mod.Get_model(),mat=dat.Get_mat_datafile();
	const int g=omega.rows(),unique=mat.rows();
	m_proba=MatrixXd::Ones(unique,g);
	m_proba_block.resize(g);
	m_param.resize(g);

	for (int k=0;k<g;k++){
		m_param[k].resize(omega.rowwise().maxCoeff()(k)+1);
		m_proba_block[k].resize(unique,omega.rowwise().maxCoeff()(k)+1);
		for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
			if ((omega.row(k).array()==b).any()){
				m_param[k][b]=param_block(k,b,dat,mod,VectorXd::Ones(mat.rows()),1);
			}
		}
	}
	m_propor=uniforme(g);
	Probapost( mod , mat );
	Compte_nbparam(dat,mod);
	Likelihood(dat.Get_eff_datafile());
	Estimation(1,0,6,dat,mod);
}

parameters::~parameters() {
	// TODO Auto-generated destructor stub
	for(int k=0;k<m_proba.cols();k++){
		m_param[k].resize(0);
	}
	m_param.resize(0);
	m_proba_block.resize(0);
}

VectorXd parameters::uniforme(int dim){
	VectorXd rep(dim);
	for (int i=0;i<dim;i++){
		rep(i)= rand()/(double)RAND_MAX ;
	}
	rep=rep.array()/rep.sum();
	return(rep);
}

void parameters::Estimation(int alter, int int_discrete,int int_continue, datafile dat, model mode){
	for (int it=0;it<alter;it++){
		//Estimation_discrete(int_discrete, dat, mode);
		Probapost(mode, dat.Get_mat_datafile());
		Likelihood(dat.Get_eff_datafile());
		Estimation_continue(int_continue,dat,mode);
	}
	Probapost(mode, dat.Get_mat_datafile());
	Likelihood(dat.Get_eff_datafile());

}

void parameters::Estimation_discrete(int nbtent, datafile dat, model mod){
	const MatrixXi & omega=mod.Get_model();
	const int g=omega.rows();
	for (int k=0;k<g;k++){
		for (int b=0;b<(omega.row(k).maxCoeff()+1);b++){
			if ((omega.row(k).array()==b).count()>1){
				m_param[k][b]=m_param[k][b].Optimise_gamma(k,b,dat,mod,nbtent, m_proba.col(k).array()/m_proba.rowwise().sum().array(),dat.Get_eff_datafile());
			}
		}
	}

}



void parameters::Estimation_continue(int nbiter,datafile dat, model mod){
	Probapost(mod, dat.Get_mat_datafile());
	for (int it=0;it<nbiter;it++){
		Mstep(dat, mod);
		Probapost(mod, dat.Get_mat_datafile());
	}
	Likelihood(dat.Get_eff_datafile());
}



void parameters::Mstep(datafile dat, model mod){
	const MatrixXi & omega=mod.Get_model(),mat=dat.Get_mat_datafile();
	const VectorXd & eff=dat.Get_eff_datafile();
	for (int k=0;k<m_proba.cols();k++){
		m_propor(k)= (eff.array()*(m_proba.col(k)).array()/m_proba.rowwise().sum().array()).sum() / eff.sum();
		for (int b=0;b<mat.cols();b++){
			if ((omega.row(k).array()==b).any()){
				const VectorXi & who=mod.Get_var_block(k,b);
				m_param[k][b].Mstep(who,mat,m_proba_block[k].col(b),m_proba.col(k).array()/m_proba.rowwise().sum().array(),eff);
			}
		}
	}
}

void parameters::Compte_nbparam(datafile dat,model mod){
	const MatrixXi & omega=mod.Get_model();
	int g=omega.rows();
	m_nbparam=g-1;
	for (int k=0;k<g;k++){
		for (int b=0;b<(omega.rowwise().maxCoeff()(k)+1) ;b++){
			if ((omega.row(k).array()==b).any()){
				const VectorXi & who=mod.Get_var_block(k,b);
				const VectorXi modalite=dat.Get_modalite(who);
				if (who.rows()==1){
					m_nbparam+=modalite(0)-1;
				}else{
					for (int h=0;h<modalite.rows();h++){
						m_nbparam+=modalite(h)-1;
					}
					m_nbparam+=modalite(0);
					if ((modalite(1)==2)&&(who.rows()==2)){m_nbparam--;}
				}
			}
		}
	}
}

void parameters::Likelihood(const VectorXd & eff){
	 VectorXd loc;
	 loc.resize(m_proba.rows());
	 loc=(m_proba.rowwise().sum().array().log());
	 m_loglikelihood=eff.transpose()*loc;
	 m_bic=m_loglikelihood - 0.5*m_nbparam*log(eff.sum());
}

void parameters::Probapost(model mod, const MatrixXi & mat){
	const MatrixXi & omega=mod.Get_model();
	m_proba=MatrixXd::Ones(mat.rows(),omega.rows());
	for (int k = 0; k < omega.rows(); ++k) {
		m_proba_block[k]=m_proba_block[k].Ones(mat.rows(),(omega.rowwise().maxCoeff()(k)+1));
		for (int b = 0; b < (omega.rowwise().maxCoeff()(k)+1); ++b){

			if ((omega.row(k).array()==b).any()){
				const VectorXi & who=mod.Get_var_block(k,b);
				if (who.rows()>1){
					m_proba_block[k].col(b)=m_param[k][b].proba_indpt(who,mat)+m_param[k][b].proba_dpt(who,mat);
				}else{
					m_proba_block[k].col(b)=m_param[k][b].proba_indpt(who,mat);
				}
			}
		}
		m_proba.col(k)=m_propor(k)* (m_proba_block[k].rowwise().prod()).array();
	}
}

double parameters::Get_bic(){
	return(m_bic);
}

double parameters::Get_like(){
	return(m_loglikelihood);
}


MatrixXd parameters::Sauve_proba_post(){
	return(m_proba);
}

int parameters::Get_nb_param(){
	return(m_nbparam);
}

vector<VectorXd> parameters::Get_alpha(int k,int b){
	return(m_param[k][b].Get_alpha());
}

VectorXd parameters::Get_tau(int k,int b){
	return(m_param[k][b].Get_tau());
}

double parameters::Get_rho(int k,int b){

	return(m_param[k][b].Get_rho());
}

MatrixXi parameters::Get_delta(int k,int b){
	return(m_param[k][b].Get_delta());
}

VectorXd parameters::Get_propor(){
	return(m_propor);
}
