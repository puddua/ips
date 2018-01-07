/**
 *file basis.cpp
 */

#include "basis.h"

#define _USE_MATH_DEFINES


/**
 *Necessary function to compute the quantic number m,corresponding to the following formula:
 *\f$ n_z^\textrm{max}(i) \equiv (N+2).Q^\frac{2}{3}+\frac{1}{2}-i.Q \f$
 *@param i
 *@param N
 *@param Q
 @return 
 */
int nmax(int i,int N,double Q){
  return ((double)N+2.0)*std::pow(Q,2.0/3.0)+0.5-(double)i*Q;

}

/**
 *Constructor of the class Basis
 *@param br
 *@param bz
 *@param N
 *@param Q
 */

Basis::Basis(double br,double bz,int N,double Q){
  int i=0;
  this->br=br;
  this->bz=bz;
  while(nmax(i,N,Q)>=1){
    i++;
  }
  mMax=i-1;
  

  
  nMax=arma::zeros<arma::ivec>(mMax);
  int m;
 
  for(m=0;m<mMax;m++){
    nMax(m)=floor(0.5*(mMax-m-1))+1;
  }

  n_zMax=arma::zeros<arma::imat>(mMax,nMax(0));
  for(m=0;m<mMax;m++){
    for(i=0;i<nMax(0);i++){
      n_zMax(m,i)=nmax(m+2*i+1,N,Q);
      if(n_zMax(m,i)<0)
	n_zMax(m,i)=0;
    }
  }

}

/**
 *Factiorial function, which use long int to avoid problems with big numbers
 *@param n
 *@result factioral(n)
 */
long fact(long n) {
  return n > 1?(n * fact(n-1)):1;
}

/**
 * Compute of rPart according to the following formula:
 * \f$ R(r_\perp, m, n)
 \equiv
 \frac{1}{b_{\perp}\sqrt{\pi}}
 \sqrt{\frac{n!}{(n+|m|)!}}
 e^{-\frac{r_{\perp}^2}{2b_{\perp}^2}}
 \left(\frac{r_{\perp}}{b_{\perp}}\right)^{|m|}
 L_n^{|m|}\left(\frac{r_{\perp}^2}{b_{\perp}^2}\right) \f$
 *@param r Vector r
 *@param m Quantic number m
 *@param n Quantic number n
 *@result Result vector
 */

arma::vec Basis::rPart(arma::vec r,int m,int n ){
  Poly pol;
  int i;
	 
  pol.calcLaguerre(m+2,n+2,arma::pow(r,2)/pow(br,2));
  arma::vec res=arma::vec(arma::size(r));
  res=(1/(br*sqrt(M_PI)))*sqrt((double)fact(n)/fact(n+abs(m)))*arma::exp(-arma::pow(r,2)/(2*pow(br,2)))%arma::pow(r/br,abs(m))%pol.laguerre(abs(m),n);
  return res;   
}

/**                                                                                                                                            
 * Compute of zPart according to the following formula: 
 * \f$ Z(z, n_z) \equiv \phi_{n_z}(z)=\frac{1}{\sqrt{b_z}} \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}}e^{-\frac{z^2}{2b_z^2}}H_{n_z}\left(\frac{z}{b_z}\right) \f$
 *@param z Vector z
 * @param nz Quantic number n_z
 *@return Result vector
 */


arma::vec Basis::zPart(arma::vec z,int nz){
  Poly pol;
  pol.calcHermite(nz+2,z/bz);
  arma::vec res=arma::vec(arma::size(z));
  /*std::cout<<"nz!="<<fact(nz)<<std::endl;
  std::cout<<"2^nz="<<pow(2,nz)<<std::endl;
  std::cout<<"cst z= "<<(1/(sqrt(bz)*sqrt(pow(2,nz)*sqrt(M_PI)*fact(nz))))<<std::endl;*/
  //  std::cout<<"gneee?"<</*pow(2,nz)*sqrt(M_PI)*/fact(nz)<<std::endl;
  res=(1/(sqrt(bz*pow(2,nz)*sqrt(M_PI)*fact(nz))))*arma::exp(-arma::pow(z,2)/(2*pow(bz,2)))%pol.hermite(nz);
  return res;
}

/**
 *Compute the wave function \f$ \psi_{m,n,n_z}(r_\perp, \theta, z)=Z(z, n_z)*R(r_\perp, m, n) \f$ 
 *@param m Quantic number m
 *@param n Quantic number n
 *@param n_z Quantic number n_z
 *@param z Vector z
 *@param r Vector r
 *@result Matrix representing the wave function
 */

arma::mat Basis::basisFunc(int m,int n,int n_z,arma::vec z,arma::vec r){
  arma::mat res=zPart(z,n_z)*(rPart(r,m,n).t());
  return res;
}


