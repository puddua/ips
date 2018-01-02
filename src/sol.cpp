#include "sol.h"


arma::mat solutionref(arma::vec zVals,arma::vec rVals,int nbR,int nbZ,Basis basis){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  for (int m = 0; m < basis.mMax; m++){
    for (int n = 0; n < basis.nMax(m); n++){
      for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++){
	j=0;
	for (int mp = 0; mp < basis.mMax; mp++){
	  for (int np = 0; np < basis.nMax(mp); np++){
	    for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++){
	      arma::mat funcA = basis.basisFunc( m,  n,  n_z, zVals, rVals);
	      arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
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


arma::mat solution1(arma::vec zVals,arma::vec rVals,int nbR,int nbZ,Basis basis){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  for (int m = 0; m < basis.mMax; m++){
    for (int n = 0; n < basis.nMax(m); n++){
      for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++){
	j=0;
	arma::mat funcA = basis.basisFunc( m,  n,  n_z, zVals, rVals);
	for (int mp = 0; mp < basis.mMax; mp++){
	  for (int np = 0; np < basis.nMax(mp); np++){
	    for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++){
	      //              arma::mat funcA = basisFunc( m,  n,  n_z, zVals, rVals);
	      arma::mat funcB = basis.basisFunc(mp, np, n_zp, zVals, rVals);
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


arma::mat solution2(arma::vec zVals,arma::vec rVals,int nbR,int nbZ,Basis basis){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  arma::mat tmpr,tmprp;
  for (int m = 0; m < basis.mMax; m++){
    for (int n = 0; n < basis.nMax(m); n++){
      tmpr=basis.rPart(rVals,m,n).t();
      for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++){
	j=0;
	arma::mat funcA = basis.zPart(zVals,n_z)*tmpr;
	for (int mp = 0; mp < basis.mMax; mp++){
	  for (int np = 0; np < basis.nMax(mp); np++){
	    tmprp=basis.rPart(rVals,mp,np).t();
	    for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++){
	      arma::mat funcB = basis.zPart(zVals,n_zp)*tmprp;
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


arma::mat solution3(arma::vec zVals,arma::vec rVals,int nbR,int nbZ,Basis basis){
  arma::mat rho;
  rho.load("rho.arma", arma::arma_ascii);
  int i=0;
  int j;
  int k=0;
  arma::mat result = arma::zeros(nbR, nbZ); // number of points on r- and z- axes
  double t1,t2,tf;
  t1=clock();
  arma::cube conv=arma::zeros(basis.mMax,basis.nMax(0),basis.n_zMax(0,0));
  for (int m = 0; m < basis.mMax; m++)
    for (int n = 0; n < basis.nMax(m); n++)
      for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
	{
	  conv(m,n,n_z)=k;
	  k++;
	}
  arma::mat tmpr,tmprp;
  for (int m = 0; m < basis.mMax; m++){
    for (int n = 0; n < basis.nMax(m); n++){
      tmpr=basis.rPart(rVals,m,n).t();
      for (int n_z = 0; n_z < basis.n_zMax(m, n); n_z++){
	arma::mat funcA = basis.zPart(zVals,n_z)*tmpr;
	for (int mp = 0; mp < basis.mMax; mp++){
	  if(mp==m){
	    for (int np = 0; np < basis.nMax(mp); np++){
	      tmprp=basis.rPart(rVals,mp,np).t();
	      for (int n_zp = 0; n_zp < basis.n_zMax(mp, np); n_zp++){
		arma::mat funcB = basis.zPart(zVals,n_zp)*tmprp;
		j=conv(m,np,n_zp);
		result+= funcA % funcB * rho(i,j); // mat += mat % mat * double
	      }
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
