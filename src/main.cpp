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
 *Function to write a wave function psi in a new file "psi.dat"           
 *@param psi Wave function to write in the file               
 */

void write_psi(arma::mat psi){
  ofstream file;
  file.open("psi.dat",ios::out);
  arma::mat plot=psi;
  file<<plot;
  file.close();
}




/**
 * The main objective of this application is to compute the local density \f$ rho \f$ of a nuclear system in 2 or 3 dimensions:
 *
 */
int main(int argc,char **argv){

  if(argc!=3){
    cerr<<"Execution requires 2 arguments: look the Readme"<<endl;
    exit(1);
  }
  
  //initialisation of vectors
  int s_r=atoi(argv[1]);
  int s_z=atoi(argv[2]);
  arma::colvec r=arma::linspace<arma::colvec>(-10,10,s_z);
  arma::colvec z=arma::linspace<arma::colvec>(-10,10,s_r);

  //initialisation of the basis
  Basis basis(1.935801664793151,      2.829683956491218,     14,     1.3);

  //compute and writing of the result
  arma::mat res=solution3(z,r,s_z,s_r,basis);
  write_result(res);

  /*arma::mat psi;
  psi=basis.basisFunc(basis.mMax/2,basis.nMax(basis.mMax/2),basis.n_zMax(basis.mMax/2,basis.mMax/2),z,r);          
  write_psi(psi);
  */
  return 0;
}

