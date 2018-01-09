#ifndef POLY_H
#define POLY_H
/**
 *@file poly.h
 */
#include <armadillo>

class Poly{
 public:
  /**
   * Matrix which contains value of the Hermite Polynomial
   */
  arma::mat her;

  /**                                                                                                    
   * Cube which contains value of the Laguerre Polynomial                                   
   */
  arma::cube guerre;

  Poly();
  void calcHermite(int,arma::vec);
  void calcLaguerre(int,int,arma::vec);  
  arma::vec laguerre(int,int);
  arma::vec hermite(int);
};

#endif
