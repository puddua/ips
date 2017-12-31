#ifndef POLY_H
#define POLY_H

#include <armadillo>

class Poly{
 public:
  arma::mat her;
  arma::cube guerre;

  Poly();
  void calcHermite(int,arma::vec);
  void calcLaguerre(int,int,arma::vec);  
  arma::vec laguerre(int,int);
  arma::vec hermite(int);
};

#endif
