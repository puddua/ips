/**                                                               
 * @file main.cpp                                                              
 */


#include <iostream>
#include <armadillo>
#include "main.h"
#include "basis.h"
//#include "hermite.h" 
using namespace std;


/**
   A modifier domi
 *Function to write value of \f$ \Psi \f$ for each degree n in a new file "result.dat"
 *@param sol Solution to write in the file
 *@param z Vector z which discretises the solution
 */

void write_result(arma::mat res,arma::colvec z,arma::colvec r){
  ofstream file;
  file.open("result.dat",ios::out);
  arma::mat plot=res;
  //  plot.insert_cols(0,z);
  //plot.insert_cols(1,r);
  file<<plot;
  file.close();
}

/**                                                                             
 osef domi
 *Function to write value of Hermite polynomials for each degree n in a new file "hermite.dat"                         
 *@param her Hermite polynomial to write in the file                                                                             
 *@param z Vector z which discretises the solution                                                                     
 */

void write_hermit(Poly pol,arma::colvec z){
  ofstream file;
  file.open("hermite.dat",ios::out);
  arma::mat plot=pol.her.t();
  plot.insert_cols(0,z);
  file<<plot;
  file.close();

}

/**
 * osef domi
 */

void write_laguerre(Poly pol,arma::colvec z,int s){
  ofstream file;
  file.open("laguerre.dat",ios::out);
  arma::mat plot=(pol.guerre.slice(s)).t();
  plot.insert_cols(0,z);
  file<<plot;
  file.close();

}



/**
 A modifier domi: faire un mini résumé sur le problème qu'on résout dans le projet
 * The main objective of this application is to solve the quantum harmonic oscillator in 1 dimension
 *For this, we solve the Schrödinger equation given by:
 *\f$ \hat{H}_{(z)} \Psi_n(z) =E_n \Psi_n (z) \f$
 * where the 1D-Hamiltonian and 1D-Momentum operator are defined as:
 *  \f$ \hat{H}_{(z)} \equiv \frac{\hat{p}^2_{(z)}}{2m}+\frac{1}{2}mw^2\hat{z}^2 \f$
 * and 
 *  \f$ \hat{p}_{(z)} \equiv -i\hbar\frac{\partial}{\partial z} \f$
 */
int main(int argc,char **argv){
  
  /*  if(argc!=3){
    cerr<<"Execution requires 2 arguments: look the Readme"<<endl;
    exit(1);
    }*/

  
  int s_z=100;
  int s_r=100;
  arma::colvec r=arma::linspace<arma::colvec>(-1,1,s_z);
  arma::colvec z=arma::linspace<arma::colvec>(-1,1,s_r);

  Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);
  arma::mat res=basis.solution1(z,r,s_z,s_r);
  write_result(res,z,r);
  
  return 0;
}

