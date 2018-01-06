/**                                                               
 * @file main.cpp                                                              
 */


#include <iostream>
#include <armadillo>
#include "main.h"
#include "basis.h"
#include "sol.h"
//#include "hermite.h" 
using namespace std;


/**
 *Function to write the result of the density in a new file "result.dat"
 *@param res Solution to write in the file 
*/

void write_result(arma::mat res){
  ofstream file;
  file.open("result.dat",ios::out);
  arma::mat plot=res;
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
 * The main objective of this application is to compute the local density of a nuclear system in 2 or 3 dimensions
 *For this, we solve the Schrödinger equation given by:
 *\f$ \hat{H}_{(z)} \Psi_n(z) =E_n \Psi_n (z) \f$
 * where the 1D-Hamiltonian and 1D-Momentum operator are defined as:
 *  \f$ \hat{H}_{(z)} \equiv \frac{\hat{p}^2_{(z)}}{2m}+\frac{1}{2}mw^2\hat{z}^2 \f$
 * and 
 *  \f$ \hat{p}_{(z)} \equiv -i\hbar\frac{\partial}{\partial z} \f$
 */
int main(int argc,char **argv){

  if(argc!=3){
    cerr<<"Execution requires 2 arguments: look the Readme"<<endl;
    exit(1);
  }
  
  int s_r=atoi(argv[1]);
  int s_z=atoi(argv[2]);
  arma::colvec r=arma::linspace<arma::colvec>(-1,1,s_z);
  arma::colvec z=arma::linspace<arma::colvec>(-1,1,s_r);

  Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);
  arma::mat res=solution3(z,r,s_z,s_r,basis);
  write_result(res);
  
  return 0;
}

