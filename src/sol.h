/**
 * @file sol.h
 */

#include "poly.h"
#include "basis.h"
#include <iostream>
#include <armadillo>

// Definiton of functions to compute the solution, with different steps of its optimization

arma::mat solutionref(arma::vec,arma::vec,int,int,Basis);
arma::mat solution1(arma::vec,arma::vec,int,int,Basis);
arma::mat solution2(arma::vec,arma::vec,int,int,Basis);
arma::mat solution3(arma::vec,arma::vec,int,int,Basis);
arma::mat solution4(arma::vec,arma::vec,int,int,Basis);
