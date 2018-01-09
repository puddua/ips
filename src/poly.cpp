/**
 *@file poly.cpp
 */
#include "poly.h"

/**
 * Constructor of the class Poly
 */
Poly::Poly(){
  
}


/**
 * Function to evaluate Hermitian polynom on a vector
 *@param n Degre of the polynom
 *@param z Vector to compute the polynom 
 */

void Poly::calcHermite(int n,arma::vec z){
  int size_z=z.n_elem;
  her=arma::zeros(n,size_z);
  arma::rowvec zz=z.t();
  her.row(0)=arma::ones(size_z).t();
  her.row(1)=2*zz;

  int i;
  for(i=2;i<n;i++){
    her.row(i)=2*zz%her.row(i-1)-2*(i-1)*her.row(i-2);
  }
}

/**
 *Return a row of the Hermitian polynomial
 *@param n Index of the row to return
 *@return Row of the Hermitian polynomials
 */
arma::vec Poly::hermite(int n){
  return her.row(n).t();
}

/**                                                                             
 * Function to evaluate Laguerre polynomials on a vector                           
 *@param n Degre of the polynom 
 *@param m 
 *@param z Vector to compute the polynom                                        

*/

void Poly::calcLaguerre(int m,int n,arma::vec z){
  int size_z=z.n_elem;
  guerre=arma::zeros(size_z,m,n);
  arma::rowvec zz=z.t();
  int i,j;
  guerre.slice(0)=arma::ones(size_z,m);
  for(i=0;i<m;i++){
    guerre.slice(1).col(i)=arma::ones(size_z)-z+i;
  }
  
  for(i=2;i<n;i++){
    for(j=0;j<m;j++){
      guerre.slice(i).col(j)=(2+(j-1-z)/i)%guerre.slice(i-1).col(j)-(1+(j-1)/(double)i)*guerre.slice(i-2).col(j);
    }
  }

}

/**                                                                                            
 *Return a row of the Laguerre polynomial
 *@param m Index of the column                                                                
 *@param n Index of the slice                                                         
 *@return Row of the Hermitian polynomials                                                                        
 */


arma::vec Poly::laguerre(int m,int n){
  return guerre.slice(n).col(m);
}
