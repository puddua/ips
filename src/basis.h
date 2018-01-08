#ifndef BASIS_H
#define BASIS_H

/**
 *@file basis.h
 */

#include "poly.h"
#include <armadillo>
#include <cmath>
#include <time.h>
#include <iostream>

/**
 *Class to compute the truncation of the basis, which define domain of definition of the result 
 */


class Basis{
 public:
  /**
   *Quantic numbers
   */
  int mMax;
  arma::ivec nMax;
  arma::imat n_zMax;
  /**
   *Parameters to define the truncation
   */
  double br;
  double bz;
  
  Basis(double,double,int,double);
  arma::vec rPart(arma::vec,int,int);
  arma::vec zPart(arma::vec,int);
  arma::mat basisFunc(int,int,int,arma::vec,arma::vec);
};

#endif
