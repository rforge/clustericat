/*
 * paramblock.cpp
 *
 *  Created on: 11 juil. 2012
 *      Author: matthieu
 */

#include "paramblock.h"
#include "datafile.h"
#include "model.h"
#include "../../eigen/Eigen/Dense"
#include <vector>
using namespace std;
using namespace Eigen;

VectorXd param_block::uniform(int dim){
	VectorXd rep(dim);
	for (int i=0;i<dim;i++){
		rep(i)= rand()/(double)RAND_MAX ;
	}
	rep=rep.array()/rep.sum();
	return(rep);
}

param_block::param_block() {
	// TODO Auto-generated constructor stub

}

/*
param_block::param_block(int k, int b, datafile dat, model mod) {
	// TODO Auto-generated constructor stub
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	if (who.rows()==1){
		m_alpha.resize(1);
		m_alpha[0]=uniform(modalite(0));
		m_epsilon=0;
	}else{
		m_alpha.resize(who.rows());
		for (int j=0;j<who.rows();j++){
			m_alpha[j]=uniform(modalite(j));
		}
		m_beta=uniform(modalite(0));
		m_epsilon=rand()/(double)RAND_MAX;
		//Genere_gamma(modalite);
		Genere_gamma2( k,  b,  dat,  mod);
		Qui_dpt(who,dat.Get_mat_datafile(),m_gamma);
	}

}
*/

param_block::param_block(int k, int b, datafile dat, model mod, VectorXd poids,int test) {
	// TODO Auto-generated constructor stub
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	if (who.rows()==1){
		m_alpha.resize(1);
		m_alpha[0]=uniform(modalite(0));
		m_epsilon=0;
	}else{
		m_alpha.resize(who.rows());
		for (int j=0;j<who.rows();j++){
			m_alpha[j]=uniform(modalite(j));
		}
		m_beta=uniform(modalite(0));
		m_epsilon=rand()/(double)RAND_MAX;
		//Genere_gamma(modalite);
		Genere_gamma( k,  b,  dat,  mod, poids);

		Qui_dpt(who,dat.Get_mat_datafile(),m_gamma);
	}

}

param_block::param_block(int k, int b, datafile dat, model mod, MatrixXi struc) {
	// TODO Auto-generated constructor stub
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	m_alpha.resize(who.rows());
	if (who.rows()==1){
		m_alpha[0]=uniform(modalite(0));
		m_epsilon=0;
	}else{
		for (int j=0;j<who.rows();j++){
			m_alpha[j]=uniform(modalite(j));
		}
		m_beta=uniform(modalite(0));
		m_epsilon=rand()/(double)RAND_MAX;
		m_gamma=struc;
		Qui_dpt(who,dat.Get_mat_datafile(),m_gamma);
	}
}



double param_block::Proba_bloc_gamma(const MatrixXi & mat, const VectorXd &  eff, const MatrixXi & gamma){
	VectorXd beta(gamma.cols());
	beta=uniform(gamma.cols());
	MatrixXd alpha=MatrixXd::Zero(2,gamma.cols());

	alpha.row(0)=uniform(gamma.cols());
	alpha.row(1).head(gamma.row(1).maxCoeff() +1)=uniform(gamma.row(1).maxCoeff() +1);
	double epsilon=rand()/(double)RAND_MAX;

	VectorXd proba_indpt=VectorXd::Constant(mat.rows(),epsilon);
	for (int j=0;j<mat.cols();j++){
		for (int i=0;i<mat.rows();i++){
			proba_indpt(i)*=alpha(j,mat(i,j)-1);
		}
	}

	VectorXd proba=proba_indpt;
	for (int i=0;i<mat.rows();i++){
		if ((mat(i,1)-1)==gamma(1,mat(i,0)-1)){
			proba(i)+=epsilon*beta(mat(i,0)-1);
		}
	}

	for (int it=0;it<2		;it++){
		VectorXd si=proba_indpt.array()/proba.array();
		alpha=alpha.array()*0;
		beta=beta.array()*0;
		for (int j=0;j<mat.cols();j++){
			for (int i=0;i<mat.rows();i++){
				alpha(j,mat(i,j)-1)+= eff(i)*si(i);
			}
			alpha.row(j)=alpha.row(j).array()/alpha.row(j).sum();
		}
		for (int i=0;i<mat.rows();i++){
			beta(mat(i,0)-1)+= eff(i)*(1-si(i));
		}
		beta=beta.array()/beta.sum();
		epsilon=(eff.array()-(eff.array()*si.array()).array()).sum()/eff.sum();

		proba_indpt=VectorXd::Constant(mat.rows(),epsilon);
		for (int j=0;j<mat.cols();j++){
			for (int i=0;i<mat.rows();i++){
				proba_indpt(i)*=alpha(j,mat(i,j)-1);
			}
		}
		proba=proba_indpt;
		for (int i=0;i<mat.rows();i++){
			if ((mat(i,1)-1)==gamma(1,mat(i,0)-1)){
				proba(i)+=epsilon*beta(mat(i,0)-1);
			}
		}
	}
	return((proba.array().log().array() * eff.array()).sum());
//	return(0.1);

}

void param_block::Genere_gamma(int k, int b, datafile dat, model mod, VectorXd poids){
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	const VectorXd & eff=dat.Get_eff_datafile(poids);
	const MatrixXi & data=dat.Get_mat_datafile();
	m_gamma.resize(modalite.rows(),modalite(0));
	m_gamma=m_gamma.array()*0;
	for (int m=0;m<modalite(0);m++){
		m_gamma(0,m)=m;
	}
	for (int j=1;j<(modalite.rows());j++){
		if ((modalite(j-1)<5)){
			const MatrixXi & possible=Genere_possible_gam(modalite(j-1),modalite(j));
			MatrixXi reduit(data.rows(),2);
			reduit.col(0)=data.col(who(j-1));
			reduit.col(1)=data.col(who(j));
			VectorXd opt(possible.rows());
			MatrixXi ga(2,modalite(j-1));
			for (int loc=0;loc<modalite(j-1);loc++){ga(0,loc)=loc;}


			for (int loc=0;loc<possible.rows();loc++){
				ga.row(1)=possible.row(loc);
				opt(loc)=Proba_bloc_gamma(reduit,eff,ga);
			}
			opt=(opt.array()-opt.maxCoeff()).array().exp();
			opt=opt.array()/opt.sum();
			for (int loc=1;loc<possible.rows();loc++){
				opt(loc)+=opt(loc-1);
			}
			double alea=(rand()/(double)RAND_MAX);
			for (int loc=0;loc<m_gamma.cols();loc++){
				m_gamma(j,loc)=possible((opt.array()<alea).count(),m_gamma(j-1,loc));
			}
		}else{
			VectorXi coord(modalite(j-1));
			for (int h=0;h<modalite(j-1);h++){
				if (modalite(j)<=h){
					coord(h)=(rand()%modalite(j));
				}else{
					coord(h)=h;
				}
			}
			for (int loc=0;loc<m_gamma.cols();loc++){
				m_gamma(j,loc)=coord(m_gamma(j-1,loc));
			}
		}
	}
}

void param_block::Genere_gamma3(int k, int b, datafile dat, model mod, VectorXd poids){
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	const VectorXd & eff=(dat.Get_eff_datafile(poids));
	const MatrixXi & data=dat.Get_mat_datafile();
	m_gamma.resize(modalite.rows(),modalite(0));
	for (int m=0;m<modalite(0);m++){
		m_gamma(0,m)=m;
	}

	Vector2d opt;
	MatrixXi g1(2,2),g2(2,2);
	g1(0,0)=0;
	g1(0,1)=1;
	g1(1,0)=0;
	g1(1,1)=1;
	g2(0,0)=0;
	g2(0,1)=1;
	g2(1,0)=1;
	g2(1,1)=0;
	for (int j=1;j<(modalite.rows());j++){
		MatrixXi reduit(data.rows(),2);
		reduit.col(0)=data.col(who(j-1));
		reduit.col(1)=data.col(who(j));
		opt(0)=Proba_bloc_gamma(reduit,eff,g1);
		opt(1)=Proba_bloc_gamma(reduit,eff,g2);

		if ((rand()/(double)RAND_MAX)<(opt.array()-opt.maxCoeff()).array().exp()(0)){
			m_gamma.row(j)=m_gamma.row(j-1);
		}else{
			m_gamma(j,0)=m_gamma(j-1,1);
			m_gamma(j,1)=m_gamma(j-1,0);
		}

	}
}

MatrixXi param_block::Genere_possible_gam(int m1,int m2){
	MatrixXi cand(pow(m2,m1),m1);
	int repete=0,li=0;
	for (int h=0;h<cand.cols();h++){
		li=0;
		while(li<(cand.rows()-1)){
			for (int val=0;val<m2;val++){
				for(int it=0;it<(pow(m2,repete));it++){
					cand(li,h)=val;
					li++;
				}
			}
		}
		repete++;
	}
	VectorXi conserve=VectorXi::Ones(cand.rows());
	for (int ligne=0;ligne<cand.rows();ligne++){
		int h=0;
		while((h<m2)&&(conserve(ligne)!=0)){
			if ((cand.row(ligne).array()!=h).all()){
				conserve(ligne)=0;
			}
			h++;
		}
	}
	MatrixXi poss((conserve.array()>0).count(),m1);
	li=0;
	int ligne=0;
	while(li<conserve.sum()){
		if (conserve(ligne)>0){
			poss.row(li)=cand.row(ligne);
			li++;
		}
		ligne++;
	}
	return(poss);
}


void param_block::Genere_gamma2(VectorXi modalite){
	m_gamma.resize(modalite.rows(),modalite(0));
	for (int h=0;h<m_gamma.cols();h++){
		m_gamma(0,h)=h;
	}
	for (int j=1;j<m_gamma.rows();j++){
		for (int h=0;h<m_gamma.cols();h++){
			if (modalite(j)<h){
				m_gamma(j,h)=(rand()%modalite(j));
			}else{
				m_gamma(j,h)=h;
			}
		}
	}
}

VectorXd param_block::proba_indpt(const VectorXi & who, const MatrixXi & mat){
	VectorXd m_proba_indpt=VectorXd::Constant(mat.rows(),(1-m_epsilon));
	for (int j=0;j<who.rows();j++){
		for (int i = 0; i < mat.rows(); ++i) {
			m_proba_indpt(i)*=m_alpha[j](mat(i,who(j))-1);
		}
	}
	return(m_proba_indpt);
}

VectorXd param_block::proba_indpt(const VectorXi & who, const MatrixXi & mat,param_block parambl){
	VectorXd m_proba_indpt=VectorXd::Constant(mat.rows(),(1-parambl.m_epsilon));
	for (int j=0;j<who.rows();j++){
		for (int i = 0; i < mat.rows(); ++i) {
			m_proba_indpt(i)*=parambl.m_alpha[j](mat(i,who(j))-1);
		}
	}
	return(m_proba_indpt);
}

VectorXd param_block::proba_dpt(const VectorXi & who,const MatrixXi & mat){
	VectorXd rep=VectorXd::Zero(mat.rows());


	for(int i=0;i<m_ind_st.rows();i++){
		rep(m_ind_st(i))= m_beta(mat(m_ind_st(i),who(0))-1);
	}
	rep=rep.array()*(m_epsilon);
	return(rep);
}

VectorXd param_block::proba_dpt(const VectorXi & who,const MatrixXi & mat,param_block parambl){
	VectorXd rep=VectorXd::Zero(mat.rows());
	for(int i=0;i<parambl.m_ind_st.rows();i++){
		rep(parambl.m_ind_st(i))= parambl.m_beta(mat(parambl.m_ind_st(i),who(0))-1);
	}
	rep=rep.array()*(parambl.m_epsilon);
	return(rep);
}

param_block::~param_block() {
	// TODO Auto-generated destructor stub
	m_alpha.resize(0);
}

void param_block::Qui_dpt(const VectorXi & who, const MatrixXi & mat, const MatrixXi & gamma){
	VectorXi qui=VectorXi::Zero(mat.rows());
	int cb=0,dpt=1;
	for (int i=0;i<mat.rows();i++){
		dpt=1;
		for (int j=1;j<who.rows();j++){
			if ((mat(i,who(j))-1)!=gamma(j,mat(i,who(0))-1)){
				dpt=0;
			}
		}
		if (dpt==1){
			qui(cb)=i;
			cb++;
		}
	}
	m_ind_st.resize(cb);
	m_ind_st=qui.head(cb);
}


double param_block::Estime_theta_cond(int k, int b, datafile dat, model mod,const VectorXd & poids,const MatrixXi & gam,const VectorXd & eff){
	param_block parambl(k,b,dat,mod,gam);
	const MatrixXi & mat=dat.Get_mat_datafile();
	const VectorXi & who=mod.Get_var_block(k,b);
	VectorXd tik=VectorXd::Ones(mat.rows());
	VectorXd pro_indptbl=proba_indpt(who,mat,parambl);
	VectorXd pro_dptbl=proba_dpt(who,mat,parambl);

	for (int it=0;it<50;it++){
		tik=pro_indptbl.array() / (pro_dptbl+pro_indptbl).array();
		parambl=actu_param(who,mat,tik,poids,eff,parambl);
		pro_dptbl=proba_dpt(who,mat,parambl);
		pro_indptbl=proba_indpt(who,mat,parambl);
	}
	return((pro_dptbl+pro_indptbl).array().log().array()*(eff.array()*poids.array())).sum();
}


param_block param_block::actu_param(const VectorXi & who,const MatrixXi & mat,const VectorXd & tik,const VectorXd & poids,const VectorXd & eff,param_block parambl){
	for (int j=0;j<who.rows();j++){
		parambl.m_alpha[j]*=0;
		 for (int i=0;i<mat.rows();i++){
			 parambl.m_alpha[j](mat(i,who(j))-1)+=tik(i)*(poids(i)*eff(i));
		 }
		 parambl.m_alpha[j]=parambl.m_alpha[j]/parambl.m_alpha[j].sum();

	}
	parambl.m_beta*=0;
	parambl.m_epsilon=0;
	for (int i=0;i<parambl.m_ind_st.rows();i++){
		parambl.m_beta(mat(parambl.m_ind_st(i),who(0))-1)+=(1-tik(parambl.m_ind_st(i)))*(poids(parambl.m_ind_st(i))*eff(parambl.m_ind_st(i)));
		parambl.m_epsilon+=(1-tik(parambl.m_ind_st(i)))*(poids(parambl.m_ind_st(i))*eff(parambl.m_ind_st(i)));
	}
	if (parambl.m_beta.sum()>0){
		parambl.m_beta=parambl.m_beta/parambl.m_beta.sum();
		parambl.m_epsilon=parambl.m_epsilon/(poids.transpose()*eff);
	}
	return(parambl);
}

param_block param_block::Optimise_gamma(int k, int b, datafile dat, model mod, int iter,const VectorXd & poids,const VectorXd & eff){
	const VectorXi & who=mod.Get_var_block(k,b);
	const VectorXi modalite=dat.Get_modalite(who);
	int j=0,h=0,h2=0,chgt;
	double actuel_likelihood_completed,cand_likelihood_completed;
	actuel_likelihood_completed=Estime_theta_cond( k, b, dat,  mod, poids, m_gamma,eff);

	for (int it=0;it<iter;it++){
		MatrixXi gamma_cand=m_gamma;
		j=1+(rand()%(modalite.rows()-1));
		h=(rand()%modalite(0));
		h2=(rand()%(modalite(0)-1));
		if( h2>=h){h2++;if (h2>modalite(0)){h2=1;}}
	//	gamma_cand.col(h).tail(who.rows()+1-j)=m_gamma.col(h2).tail(who.rows()+1-j);
	//	gamma_cand.col(h2).tail(who.rows()+1-j)=m_gamma.col(h).tail(who.rows()+1-j);
		gamma_cand(j,h)=m_gamma(j,h2);
		gamma_cand(j,h2)=m_gamma(j,h);
		cand_likelihood_completed=Estime_theta_cond( k, b, dat,  mod, poids,gamma_cand,eff);
		if ((rand()/(double)RAND_MAX) < 1/(1+exp(actuel_likelihood_completed-cand_likelihood_completed))){
			m_gamma=gamma_cand;
			chgt=1;
		}

	}
	// Faire une estimation correcte de Theta associé au Gamma !!!!
	param_block ret(k,b,dat,mod,m_gamma);
	return(ret);
}


void param_block::Mstep(const VectorXi & who, const MatrixXi & mat, const VectorXd  & proba_block, const VectorXd & m_proba, const VectorXd & eff){
	if (who.rows()>1){
		VectorXd proba_dep=proba_dpt(who, mat).array()/proba_block.array();
		m_beta*=0;
		m_epsilon=0;
		for (int j=0;j<who.rows();j++){
			m_alpha[j]=VectorXd::Zero(m_alpha[j].rows());
		}
	//	vector<VectorXd> sauv=m_alpha;
		for (int i=0;i<mat.rows();i++){
			 if (proba_block(i)==0){
				 proba_dep(i)=0;
			 }
			for (int j=0;j<who.rows();j++){
				 m_alpha[j](mat(i,who(j))-1)+=(1-proba_dep(i))*(m_proba(i)*eff(i));
			}
		}
		for (int i=0;i<m_ind_st.rows();i++){
			m_beta(mat(m_ind_st(i),who(0))-1)+=(proba_dep(m_ind_st(i)))*(m_proba(m_ind_st(i))*eff(m_ind_st(i)));
			m_epsilon+=(proba_dep(m_ind_st(i)))*(m_proba(m_ind_st(i))*eff(m_ind_st(i)));
		}
		for (int j=0;j<who.rows();j++){
				 m_alpha[j]=m_alpha[j]/m_alpha[j].sum();
		}
		if (m_beta.sum()>0){
			m_beta=m_beta/m_beta.sum();
			m_epsilon=m_epsilon/(m_proba.transpose()*eff);
		}
	}else{
		for (int j=0;j<who.rows();j++){
			m_alpha[j]*=0;
			 for (int i=0;i<mat.rows();i++){
				 m_alpha[j](mat(i,who(j))-1)+=(m_proba(i)*eff(i));
			 }
			 m_alpha[j]=m_alpha[j]/m_alpha[j].sum();
		}
	}
}


int param_block::isNaN(){
	int NaN=0;
	for (int j=0;j<m_alpha.size();j++){
		for (int loc=0; loc<m_alpha[j].rows();loc++){
			if (isnan(m_alpha[j](loc))){
				NaN=1;
			}
		}
	}
	if (isnan(m_epsilon)){
		NaN=1;
	}
	for(int loc=0; loc<m_beta.rows();loc++){
		if (isnan(m_beta(loc))){
			NaN=1;
		}
	}
	return(NaN);
}

vector<VectorXd> param_block::Get_alpha(){
	return(m_alpha);
}

VectorXd param_block::Get_tau(){
	return(m_beta);
}

double param_block::Get_rho(){
	return(m_epsilon);
}

MatrixXi param_block::Get_delta(){
	return(m_gamma);
}
