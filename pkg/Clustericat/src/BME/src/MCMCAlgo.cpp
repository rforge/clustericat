/*
 * MCMCAlgo.cpp
 *
 *  Created on: 12 juil. 2012
 *      Author: matthieu
 */

#include "MCMCAlgo.h"
#include "../../eigen/Eigen/Dense"
#include <vector>
#include "datafile.h"
#include "model.h"
//#include "paramblock.h"
#include "parameters.h"
#include "math.h"

MCMCAlgo::MCMCAlgo(){

}


MCMCAlgo::MCMCAlgo(datafile dat,MatrixXi mo, int nbc,int nbinit,int nbiter,int iterdiscret,int itercond,int iterglobal) {
	// TODO Auto-generated constructor stub
	model best_mod(mo,nbc);
	model actuel_mod=best_mod;
	parameters best_param(dat,best_mod);
	best_param.Estimation(nbiter,iterdiscret,itercond,dat,best_mod);
	parameters actuel_param=best_param;
	int compt=0,k=0,b1=0,b2=0,bnv=0,vu=0;


	while(compt<iterglobal){
		const MatrixXi & actuel_omega=actuel_mod.Get_model();
		MatrixXi candidats_omega=actuel_omega;
		compt++;
		vector<parameters> candidats_param;
		vector<model> candidats_mod;
		VectorXd proba_tirage;
		k=(rand()%actuel_omega.rows());
		b1= actuel_omega(k,rand()%actuel_omega.cols());
		b2= actuel_omega(k,rand()%actuel_omega.cols());
		bnv=0;
		vu=0;
		while((actuel_omega.row(k).array()==bnv).any()){bnv++;}
		if (bnv==actuel_omega.cols()){
			proba_tirage=VectorXd::Zero(2);
			candidats_param.resize(1);
			candidats_mod.resize(1);
			while(actuel_omega(k,vu)!=b1){vu++;}
			candidats_omega=actuel_omega;
			candidats_omega(k,vu)=b2;
			candidats_mod[0]=model(candidats_omega,nbc);
			candidats_param[0]=parameters(dat, candidats_mod[0],actuel_mod,actuel_param,k);
			proba_tirage(0)=candidats_param[0].Get_bic();
			proba_tirage(1)=actuel_param.Get_bic();
			if ((rand()/(double)RAND_MAX) < 1/(1+exp(proba_tirage(1)-proba_tirage(0)))){
				actuel_mod=candidats_mod[0];
				actuel_param=candidats_param[0];
				if (actuel_param.Get_bic()>best_param.Get_bic()){
					if (actuel_param.Get_bic()>(best_param.Get_bic()+0.001)){
						compt=0;
					}
					best_mod=actuel_mod;
					best_param=actuel_param;
				}
			}

		}else{
			int cb=(actuel_omega.row(k).array()==b1).count(),var=0,indice=0;
			proba_tirage=VectorXd::Zero(2*cb + 2);
			candidats_param.resize(2*cb +1);
			candidats_mod.resize(2*cb +1);
			VectorXi qui=VectorXi::Zero(cb);
			vu=0;
			while (vu<cb){
				if (actuel_omega(k,var)==b1){
					candidats_omega=actuel_omega;
					candidats_omega(k,var)=b2;
					candidats_mod[indice]=model(candidats_omega,nbc);
					candidats_param[indice]=parameters(dat, candidats_mod[indice],actuel_mod,actuel_param,k);//,1+int(10*compt/iterglobal));
					proba_tirage(indice)=candidats_param[indice].Get_bic();
					indice++;
					candidats_omega=actuel_omega;
					candidats_omega(k,var)=bnv;
					candidats_mod[indice]=model(candidats_omega,nbc);
					candidats_param[indice]=parameters(dat,candidats_mod[indice],actuel_mod,actuel_param,k);//,1+int(10*compt/iterglobal));
					proba_tirage(indice)=candidats_param[indice].Get_bic();
					indice++;
					vu++;
				}
				var++;
			}
			candidats_omega=actuel_omega;
			vu=0;
			var=0;
			while (vu<cb){
				if (actuel_omega(k,var)==b1){
					candidats_omega(k,var)=bnv;
					vu++;
				}
				var++;
			}
			candidats_mod[indice]=model(candidats_omega,nbc);
			candidats_param[indice]=parameters(dat, candidats_mod[indice],actuel_mod,actuel_param,k);//,1+int(10*compt/iterglobal));
			proba_tirage(indice)=candidats_param[indice].Get_bic();
			indice++;
			actuel_param.Estimation_continue(1+int(10*compt/iterglobal),dat,actuel_mod);
			proba_tirage(indice)=actuel_param.Get_bic();
			proba_tirage=(proba_tirage.array()-proba_tirage.maxCoeff()).array().exp();
			for (int loc=1;loc<proba_tirage.rows();loc++){proba_tirage(loc)+=proba_tirage(loc-1);}
			proba_tirage=proba_tirage.array()/proba_tirage(proba_tirage.rows()-1);
			double alea=(rand()/(double)RAND_MAX);
			if ((proba_tirage.array()<alea).count()<indice){
				actuel_mod=candidats_mod[(proba_tirage.array()<alea).count()];
				actuel_param=candidats_param[(proba_tirage.array()<alea).count()];
				if (actuel_param.Get_bic()>best_param.Get_bic()){
					if (actuel_param.Get_bic()>(best_param.Get_bic()+0.001)){
						compt=0;
					}
					best_mod=actuel_mod;
					best_param=actuel_param;
				}
			}
		}
		candidats_mod.resize(0);
		candidats_param.resize(0);
	}


	best_model=best_mod;
	best_parametres=best_param;

}

MCMCAlgo::~MCMCAlgo() {
	// TODO Auto-generated destructor stub
}


double MCMCAlgo::Get_best_bic(){
	return(best_parametres.Get_bic());
}
double MCMCAlgo::Get_best_like(){
	return(best_parametres.Get_like());
}


MatrixXi MCMCAlgo::Get_best_omega(){
	return(best_model.Get_model());
}


MatrixXd MCMCAlgo::Sauv_probapost(){
	return(best_parametres.Sauve_proba_post());
}

int MCMCAlgo::Get_nb_param(){
	return(best_parametres.Get_nb_param());
}


vector<VectorXd> MCMCAlgo::Get_alpha(int k,int b){
	return(best_parametres.Get_alpha(k,b));
}

VectorXd MCMCAlgo::Get_tau(int k,int b){
	return(best_parametres.Get_tau(k,b));
}

double MCMCAlgo::Get_rho(int k,int b){
	return(best_parametres.Get_rho(k,b));
}

MatrixXi MCMCAlgo::Get_delta(int k,int b){
	return(best_parametres.Get_delta(k,b));
}

VectorXd MCMCAlgo::Get_propor(){
	return(best_parametres.Get_propor());
}
