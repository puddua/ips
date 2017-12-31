/**
 *file basis.cpp
 */

#include "basis.h"

#define _USE_MATH_DEFINES

int nmax(int i,int N,double Q){
  return ((double)N+2.0)*std::pow(Q,2.0/3.0)+0.5-(double)i*Q;

}

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

long fact(long n) {
  /*int i,factorial;
  for (i = 0; i <= n; i++){

    if (i == 0)
      factorial = 1;
    else
      factorial = factorial * i;
  }
  return factorial;*/
  return n > 1?(n * fact(n-1)):1;
}

arma::vec Basis::rPart(arma::vec r,int m,int n ){
  Poly pol;
  int i;
	 
  pol.calcLaguerre(m+2,n+2,arma::pow(r,2)/pow(br,2));
  arma::vec res=arma::vec(arma::size(r));
  // std::cout<<"n!="<<fact(n)<<"n+m!="<<fact(n+abs(m))<<std::endl;
  //std::cout<<"cst r= "<<(1/(br*sqrt(M_PI)))*sqrt((double)fact(n)/fact(n+abs(m)))<<std::endl;
  res=(1/(br*sqrt(M_PI)))*sqrt((double)fact(n)/fact(n+abs(m)))*arma::exp(-arma::pow(r,2)/(2*pow(br,2)))%arma::pow(r/br,abs(m))%pol.laguerre(abs(m),n);
  return res;   
}

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

arma::mat Basis::basisFunc(int m,int n,int n_z,arma::vec z,arma::vec r){
  arma::mat res=zPart(z,n_z)*(rPart(r,m,n).t());
  return res;
}


arma::mat Basis::solutionref(arma::vec zVals,arma::vec rVals,int nbR,int nbZ){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  for (int m = 0; m < mMax; m++){
      for (int n = 0; n < nMax(m); n++){
	  for (int n_z = 0; n_z < n_zMax(m, n); n_z++){
	      j=0;
	      for (int mp = 0; mp < mMax; mp++){
		  for (int np = 0; np < nMax(mp); np++){
		      for (int n_zp = 0; n_zp < n_zMax(mp, np); n_zp++){
			  arma::mat funcA = basisFunc( m,  n,  n_z, zVals, rVals);
			  arma::mat funcB = basisFunc(mp, np, n_zp, zVals, rVals);
			  result+= funcA % funcB * rho(i,j); // mat += mat % mat * double
			  j++;
			}
		    }
		}
	      i++;
	    }
	}
    }
  t2=clock();
  tf=(double)(t2-t1)/CLOCKS_PER_SEC;
  std::cout<<"temps d'execution = "<<tf<<std::endl;
  return result;
}


arma::mat Basis::solution1(arma::vec zVals,arma::vec rVals,int nbR,int nbZ){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  for (int m = 0; m < mMax; m++){
    for (int n = 0; n < nMax(m); n++){
      for (int n_z = 0; n_z < n_zMax(m, n); n_z++){
	j=0;
	//for (int mp = 0; mp < mMax; mp++){
	//for (int np = 0; np < nMax(mp); np++){
	//  for (int n_zp = 0; n_zp < n_zMax(mp, np); n_zp++){
	      arma::mat funcA = basisFunc( m,  n,  n_z, zVals, rVals);
	      arma::mat funcB = basisFunc(m, n, n_z, zVals, rVals);
	      result+= funcA % funcB * rho(i,j); // mat += mat % mat * double
	      j++;
	      /*	    }
	  }
	  }*/
	i++;
      }
    }
  }
  t2=clock();
  tf=(double)(t2-t1)/CLOCKS_PER_SEC;
  std::cout<<"temps d'execution = "<<tf<<std::endl;
  return result;
}

arma::mat Basis::solution2(arma::vec zVals,arma::vec rVals,int nbR,int nbZ){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  for (int m = 0; m < mMax; m++){
    for (int n = 0; n < nMax(m); n++){
      for (int n_z = 0; n_z < n_zMax(m, n); n_z++){
	j=0;
	for (int mp = 0; mp < mMax; mp++){
	  for (int np = 0; np < nMax(mp); np++){
	    for (int n_zp = 0; n_zp < n_zMax(mp, np); n_zp++){
	      arma::mat funcA = basisFunc( m,  n,  n_z, zVals, rVals);
	      arma::mat funcB = basisFunc(mp, np, n_zp, zVals, rVals);
	      result+= funcA % funcB * rho(i,j); // mat += mat % mat * double
	      j++;
	    }
	  }
	}
	i++;
      }
    }
  }
  t2=clock();
  tf=(double)(t2-t1)/CLOCKS_PER_SEC;
  std::cout<<"temps d'execution = "<<tf<<std::endl;
  return result;
}



