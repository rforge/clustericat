/*
 * rcpp_reading.h
 *
 *  Created on: 2 oct. 2012
 *      Author: matthieu
 */

#ifndef _BMEcategorical_RCPP_READING_H_
#define _BMEcategorical_RCPP_READING_H_

#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP BMEclustering(SEXP mat,SEXP moda, SEXP nb_cluster, SEXP partition_initiale, SEXP nb_init, SEXP stop_criterion) ;

#endif /* RCPP_READING_H_ */
