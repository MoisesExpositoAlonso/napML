#include <stdlib.h>
#include <cstdio>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <math.h>
#include <vector>
#include <list>
#include <string>
#include <iostream>

#include <fstream>
#include <sstream>


// #define ARMA_64BIT_WORD 1
// //// https://stackoverflow.com/questions/40592054/large-matrices-in-rcpparmadillo-via-the-arma-64bit-word-define

// when armadillo is loaded, remove this below
//#include <Rcpp.h>
#include <RcppArmadillo.h>

#include <bigmemory/MatrixAccessor.hpp>
#include <bigmemory/isna.hpp>


// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(bigmemory)]]
// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppArmadillo)]]

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// #include "utils.h"

////////////////////////////////////////////////////////////////////////////////
/// Constants
////////////////////////////////////////////////////////////////////////////////

// #define MIN_NUM = std::numeric_limits<float>::min(); // problem is that it does not know the type
const double MIN_NUM = std::numeric_limits<float>::min();
// #define PIraw = 3.14159265358979323846;
// const double PI= 3.14159265358979323846;
// # define PI 3.14159265358979323846  /* pi */
//const double MAX_NUM = std::numeric_limits<float>::max();

double LLMIN = 1.0e+99;

////////////////////////////////////////////////////////////////////////////////
/// PLINK
////////////////////////////////////////////////////////////////////////////////

#define PACK_DENSITY 4
#define PLINK_NA 3

#define PLINK_PHENO_MISSING -9

// The BED file magic numbers
#define PLINK_OFFSET 3

#define COVAR_ACTION_TRAIN_TEST 0
#define COVAR_ACTION_TRAIN_ONLY 1


/* 3 is 11 in binary, we need a 2 bit mask for each of the 4 positions */
#define MASK0 3	  // 3 << 2 * 0 //
#define MASK1 12  // 3 << 2 * 1 //
#define MASK2 48  // 3 << 2 * 2 //
#define MASK3 192 // 3 << 2 * 3 //

#define BUFSIZE 100

/*
 *                   plink BED           sparsnp
 * minor homozyous:  00 => numeric 0     10 => numeric 2
 * heterozygous:     10 => numeric 2     01 => numeric 1
 * major homozygous: 11 => numeric 3     00 => numeric 0
 * missing:          01 => numeric 1     11 => numeric 3
 *
 *
 * http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml says,
 * The bytes in plink are read backwards HGFEDCBA, not GHEFCDAB, but we read
 * them forwards as a character (a proper byte)
 *
 * By default, plink usage dosage of the *major* allele, since allele A1 is
 * usually the minor allele and the code "1" refers to the second allele A2,
 * so that "11" is A2/A2 or major/major.
 *
 * We always use minor allele dosage, to be consistent with the output from
 * plink --recodeA which used minor allele dosage by default.
 *
 * out: array of genotypes
 * in: array of packed genotypes (bytes)
 * n: number of bytes in input
 *
 */
void decode_plink(unsigned char *out,
      const unsigned char *in, const unsigned int n)
{
   unsigned int i, k;
   unsigned char tmp, geno;
   unsigned int a1, a2;

   for(i = 0 ; i < n ; ++i){
      tmp = in[i];
      k = PACK_DENSITY * i;

      /* geno is interpreted as a char, however a1 and a2 are bits for allele 1 and
       * allele 2. The final genotype is the sum of the alleles, except for 01
       * which denotes missing.
       */
      geno = (tmp & MASK0);
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK1) >> 2;
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK2) >> 4;
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
      k++;

      geno = (tmp & MASK3) >> 6;
      a1 = !(geno & 1);
      a2 = !(geno >> 1);
      out[k] = (geno == 1) ? 3 : a1 + a2;
   }
}

// [[Rcpp::export]]
int BMreadbed(std::string bedfile,
              SEXP G,
              unsigned int N=1011, // number of individuals
              unsigned int p=83794, // number of SNPs
              bool verbose=true){
   ////////////////////////////////////////////////////
   // Opening backing file
   Rcpp::XPtr<BigMatrix> bigMat(G);
   MatrixAccessor<double> macc(*bigMat);
   ////////////////////////////////////////////////////
   // Initializing
   unsigned long len;
   unsigned int np, nsnps, ncovar;
  ////////////////////////////////////////////////////
  // Opening binary file
   std::ifstream in(bedfile, std::ios::in | std::ios::binary);
   if(!in){
      std::cerr << "[read_bed] Error reading file " << bedfile << std::endl;
      throw std::runtime_error("io error");
   }
   in.seekg(0, std::ifstream::end);
   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;
   // std::cout << "The size in bytes of the file is: " << len<< std::endl;
  ////////////////////////////////////////////////////
  // Iterating

   // size of packed data, in bytes, per SNP
   np = (unsigned int)std::ceil((double)N / (double)PACK_DENSITY);
   // std::cout << "Size in bytes of SNP packs is: " << np << std::endl;
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

  // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
    for(unsigned int j = 0 ; j < nsnps ; j++){
    in.read((char*)tmp, sizeof(char) * np);
    decode_plink(tmp2, tmp, np);
      for(unsigned int i = 0 ; i < N ; i++){
        // std::cout << (double)tmp2[i] << std::endl;
        macc[j][i] = (double)tmp2[i];
      }
    }
   in.close();
  return 1;
}


// [[Rcpp::export]]
arma::Mat<double> readbed(std::string bedfile,
                          int N, int p,
                          const arma::uvec & myrows,
                          const arma::uvec & mycols,
                          bool verbose=true){
   ////////////////////////////////////////////////////
   // Opening backing file
   arma::Mat<double> X(N,p);
   ////////////////////////////////////////////////////
   // Initializing
   unsigned long len;
   unsigned int np, nsnps, ncovar;
  ////////////////////////////////////////////////////
  // Opening binary file
   std::ifstream in(bedfile, std::ios::in | std::ios::binary);
   if(!in){
      std::cerr << "[read_bed] Error reading file " << bedfile << std::endl;
      throw std::runtime_error("io error");
   }
   in.seekg(0, std::ifstream::end);
   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;
   // std::cout << "The size in bytes of the file is: " << len<< std::endl;

  ////////////////////////////////////////////////////
  // Iterating
   // size of packed data, in bytes, per SNP
   np = (unsigned int)std::ceil((double)N / (double)PACK_DENSITY);
   // std::cout << "Size in bytes of SNP packs is: " << np << std::endl;
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

  // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
   double val=0;
    for(unsigned int j = 0 ; j < nsnps ; j++){
      arma::vec tmp3(N,  arma::fill::zeros);
      in.read((char*)tmp, sizeof(char) * np);
      decode_plink(tmp2, tmp, np);
        for(unsigned int i = 0 ; i < N ; i++){
          val=(double)tmp2[i];
          if(val != PLINK_NA) tmp3[i] = val; // default zeroes
        }
      X.col(j)=tmp3;
    }
   in.close();

    // Subset matrix
    if(myrows.n_elem == X.n_rows){
      X=X.cols(mycols-1);
    }else if(mycols.n_elem == X.n_rows){
      X=X.rows(myrows-1);
    }else{
      X=X.submat(myrows-1,mycols-1);
    }

  return X;
}


// [[Rcpp::export]]
int printbed(std::string bedfile,
              unsigned int N=1011, // number of individuals
              unsigned int p=83794, // number of SNPs
              bool verbose=true){
   ////////////////////////////////////////////////////
   // Initializing
   unsigned long len;
   unsigned int np, nsnps, ncovar;
  ////////////////////////////////////////////////////
  // Opening binary file
   std::ifstream in(bedfile, std::ios::in | std::ios::binary);
   if(!in){
      std::cerr << "[read_bed] Error reading file " << bedfile << std::endl;
      throw std::runtime_error("io error");
   }
   in.seekg(0, std::ifstream::end);
   // file size in bytes, ignoring first 3 bytes (2byte magic number + 1byte mode)
   len = (unsigned int)in.tellg() - 3;
   std::cout << "The size in bytes of the file is: " << len<< std::endl;
  ////////////////////////////////////////////////////
  // Iterating

   // size of packed data, in bytes, per SNP
   np = (unsigned int)std::ceil((double)N / (double)PACK_DENSITY);
   std::cout << "Size in bytes of SNP packs is: " << np << std::endl;
   nsnps = len / np;
   in.seekg(3, std::ifstream::beg);

   unsigned char* tmp = new unsigned char[np];

 // Allocate more than the sample size since data must take up whole bytes
   unsigned char* tmp2 = new unsigned char[np * PACK_DENSITY];
   // Rcpp::NumericVector tmp3(N); // in original used Eigen

    for(unsigned int j = 0 ; j < nsnps ; j++){ // but only first 10 blocks
    in.read((char*)tmp, sizeof(char) * np);
    // std::cout<< tmp << std::endl; // prints as cat in terminal
    decode_plink(tmp2, tmp, np);
      for(unsigned int i = 0 ; i < N ; i++){ // go over 10 individuals individuals
        // tmp3(i) = (double)tmp2[i];
        std::cout << (double)tmp2[i] << std::endl;
      }
    }
   in.close();
  return 1;
}



////////////////////////////////////////////////////////////////////////////////
/// Utilities with matrices
////////////////////////////////////////////////////////////////////////////////



// [[Rcpp::export]]
bool BMrecode(SEXP A, double & a, const double & b ){
    Rcpp::XPtr<BigMatrix> bigMat(A);
    MatrixAccessor<double> macc(*bigMat);
    int i, j;
    double val=0;
    for (j = 0; j <bigMat->ncol(); j++) {
      for (i = 0; i < bigMat->nrow(); i++) {
        val= (macc[j][i]);
        if(val == a){
        macc[j][i] = b;
       }
      }
    }
    return Rcpp::wrap(true);
}

// [[Rcpp::export]]
bool BMsimulate(SEXP A){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);

      Rcpp::NumericVector maf = Rcpp::runif( bigMat->ncol(), 0,0.49);

      int i, j;
      for (j = 0; j <bigMat->ncol(); j++) {
        for (i = 0; i < bigMat->nrow(); i++) {
          if(Rcpp::runif(1)(0) < maf(j)) macc[j][i] = 2;
        }
      }
      return Rcpp::wrap(true);
}

// [[Rcpp::export]]
bool BMwrite012(SEXP A, std::string outfile){
     // Accessor
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);
      // Open file
      ofstream myfile;
      myfile.open(outfile);

      // header looks like ._A
      std::string head;
      std::string spacer=" ";
      for(int l=0;l< bigMat->ncol(); l++) head.append("._A");
      myfile<<head<<"\n";

      // the matrix
      for (int i = 0; i < bigMat->nrow(); i++) {
        for (int j = 0; j <bigMat->ncol(); j++) {
          myfile<<macc[j][i]<<" ";
        }
        myfile<<"\n";
      }
      return Rcpp::wrap(true);
}


// [[Rcpp::export]]
bool BMwritePED(SEXP A, std::string outfile){
     // Accessor
      Rcpp::XPtr<BigMatrix> bigMat(A);
      MatrixAccessor<double> macc(*bigMat);
      // Open file
      ofstream myfile;
      myfile.open(outfile);

      // the matrix
      double val=0;
      for (int i = 0; i < bigMat->nrow(); i++) {
        for (int j = 0; j <bigMat->ncol(); j++) {
          val= (macc[j][i]);
          if (val==0) myfile<< "C C" <<" ";
          else if (val==1) myfile<< "C A" <<" ";
          else if (val==2) myfile<< "A A" <<" ";
          else myfile<< "N N" <<" ";
        }
        myfile<<"\n";
      }
      return Rcpp::wrap(true);
}

// [[Rcpp::export]]
arma::Mat<double> BMsubset(SEXP A,
                           const arma::uvec & myrows,
                           const arma::uvec & mycols ){
      Rcpp::XPtr<BigMatrix> bigMat(A);
      arma::Mat<double> X0((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
                                      // consider saying true, perhaps is faster
      // Subset matrix
    	if(myrows.n_elem == X0.n_rows){
    		X0=X0.cols(mycols-1); // BUG CHECK - position 0 vs 1 in R
    	}else if(mycols.n_elem == X0.n_rows){
    	  X0=X0.rows(myrows-1);// BUG CHECK - position 0 vs 1 in R
    	}else{
    		X0=X0.submat(myrows-1,mycols-1);// BUG CHECK - position 0 vs 1 in R
    	}
      return(X0);
}

// [[Rcpp::export]]
arma::vec upperTmat(const arma::mat mymat){
  arma::vec output((mymat.n_cols*(mymat.n_cols-1))/2);
  arma::mat::const_iterator it = mymat.begin() + mymat.n_rows; //iterator at rows to skip in every step, starts at second column
  long toSkipInVec = 0;
  for(int i = 1; i < mymat.n_cols; i++) //Starts with 1 to skip the diagonal
  {
    std::copy(it, it + i, output.begin() + toSkipInVec);
    toSkipInVec += i;
    it += mymat.n_rows;
  }
  return output;
}

// [[Rcpp::export]]
arma::mat Xmcenter(arma::mat X){
  arma::mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) ;
  }
  return(newX);
}

// [[Rcpp::export]]
arma::mat Xmvcenter(arma::mat X){
  arma::mat newX(X.n_rows,X.n_cols);
  for(int j=0; j<X.n_cols; j++){
   newX.col(j) = (X.col(j) - arma::mean( X.col(j))) /arma::stddev( X.col(j));
  }
  return(newX);
}

// [[Rcpp::export]]
arma::mat LDrelative(SEXP A, arma::uvec  mycols, bool debug = false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  if(bigMat->matrix_type() !=8) Rcpp::stop("Big matrix is not of type double");

  // Read the genome matrix from address
  arma::Mat<double> X((double*) bigMat->matrix(), bigMat->nrow(), bigMat->ncol(), false, false);
  X=X.cols(mycols);

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}

arma::mat LDrelative(arma::mat X, bool debug = false){

  // mean and var center for LD calculation
  X=Xmvcenter(X);
  if(debug) cout << X << endl;

  // Get the relative LD for the proposals
  arma::mat R2 =  arma::trans(X)*X ;
  if(debug)  cout << R2 << endl;
  R2 = R2/ arma::sum(upperTmat(R2));
  if(debug)  cout << arma::sum(upperTmat(R2)) << endl;
  return(R2);
}

////////////////////////////////////////////////////////////////////////////////
/// GWA
////////////////////////////////////////////////////////////////////////////////

// // [[Rcpp::export()]]
// arma::vec rankC(arma::vec x) {
//   // not a typo, rankC -> order R equivalent rankC(rankC) -> rank R equivalent
//     arma::uvec rankarma  = arma::stable_sort_index(x);
//     arma::vec z = arma::conv_to<arma::vec>::from(rankarma);
//     return z;
// }


// // [[Rcpp::export]]
// List accuracies (const arma::vec & y,
//                          const arma::vec & what){
//   // linear model
//   arma::mat xs(what.n_elem, 2);
//   xs.fill(1);
//   xs.col(1)=y;
//   arma::colvec coef = arma::solve(xs, what);
//   double R2= 1- (sum(pow(y-xs*coef,2))/
//                  sum(pow(y-arma::mean(y),2))
//                  );
//   // Pearson r
//   arma::mat m1;
//   m1.insert_cols(m1.n_cols, y);
//   m1.insert_cols(m1.n_cols, what);
//   arma::mat r_=arma::cor(m1);
//   double r = r_(0,1);
//   // Spearman rho
//   arma::mat m2;
//   m2.insert_cols(m2.n_cols, rankC(rankC(y))); // not a typo, rankC -> order R equivalent rankC(rankC) -> rank R equivalent
//   m2.insert_cols(m2.n_cols, rankC(rankC(what)));
//   arma::mat rho_=arma::cor(m2);
//   double rho = rho_(0,1);
//
//   return List::create(Named("a") = coef(0),
//                       Named("b")=coef(1),
//                       Named("R2")=R2,
//                       Named("rho")=rho,
//                       Named("r")=r
//                       );
// }


// [[Rcpp::export]]
arma::colvec My(const arma::colvec & y, const arma::colvec & h){
  /*
  * Mean trait per genotype
  */
  // Declarations
  arma::colvec hunique = unique(h);
  arma::colvec m(hunique.n_elem);
  // cout << hunique << endl;

  for(int i=0; i< hunique.n_elem; i++ ){
    // Create temporal vector
    arma::colvec ytmp;
    // Fill with all values corresponding to the same genotype
    for(int j=0; j<y.n_elem;j++){
      if(h(j) == hunique(i)) {
        ytmp.resize(ytmp.size()+1);
        ytmp(ytmp.size()-1) = y(j);
      }
    }
    // Compute variance
      if(ytmp.n_elem ==1){
       // v(i)=0;
       	m(i)=ytmp(0);
      }else{
      	m(i)=arma::mean( ytmp );
      }
  }
  return(m);
}


// [[Rcpp::export]]
arma::colvec BMs(const SEXP A,
                 const arma::colvec & y,
                 const arma::uvec & mycols,
                 const arma::uvec & myrows,
                 bool debug=false) {
  Rcpp::XPtr<BigMatrix> bigMat(A);
  MatrixAccessor<double> macc(*bigMat);
  arma::vec coef (mycols.n_elem);

  for (int j = 0; j <mycols.n_elem; j++) {
    int n1=0, n0=0;
    double m1=0, m0=0;
      for (int i = 0; i < myrows.n_elem; i++) {
        if(macc[mycols(j)-1][myrows(i)-1] == 0){
          n0++;
          m0 += y(myrows(i)-1);
        }else{
          n1++;
          m1 += y(myrows(i)-1);
        }
      }
      coef(j) = (m1/n1) - (m0/n0);
      if(isna(coef(j))) coef(j)=0;
  }
  double sds = sqrt(arma::var(coef));
  for (int j = 0; j <mycols.n_elem; j++) if(abs(coef(j))>sds*3) coef(j)=0;
  for (int j = 0; j <mycols.n_elem; j++) if(coef(j) < -1) coef(j)=-0.999;
  return coef;
}


// [[Rcpp::export]]
arma::colvec BMmgwa( arma::mat X0,
                     arma::colvec & yraw,
                    bool debug=false) {
      /////////////////////////
 	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;
  // Ordineary Least Squares
  arma::colvec coef(X.n_cols);
	arma::colvec val;
	// cout << "calculate effects" << endl;
  for(int i=0; i< X.n_cols; i++){
		arma::mat Xsub= X.col(i);
		arma::vec val = solve(Xsub,arma::mat(y));
		if(debug) cout << val << endl;
		coef(i) = val(0);
	}
 return(coef);
}



// [[Rcpp::export]]
arma::vec BMridge(arma::mat X0,const arma::colvec & yraw,
             const double & lambda=1) {
 	// Centering
  arma::mat X=Xmcenter(X0);
  arma::vec y = (yraw/arma::mean(yraw) ) - 1;
  // Precompute some values
	int n = X.n_rows;
	int p = X.n_cols;

	arma::vec coef(p, arma::fill::zeros);
	arma::colvec y_tilde = join_cols(y, arma::zeros<arma::colvec>(p));

		arma::mat X_tilde = join_cols(X, sqrt(lambda) * arma::eye(p, p));
		arma::mat Q, R;
		qr_econ(Q, R, X_tilde);
		coef = solve(R, Q.t() * y_tilde);

  return coef;
	// return List::create(Named("coef") = coef,
	//           					Named("lambda") = lambda);
}

////////////////////////////////////////////////////////////////////////////////
/// Likelihood, Probabilities, Proposals
////////////////////////////////////////////////////////////////////////////////

////// PRIOR
// double uniform(const arma::colvec & s, const double min, const double max){ // uniform probability
//   double L= 0;
//   int N=s.n_elem;
//   for(int i=0;i<N ;i++){
//     L+= R::dunif(s(i),min,max,true);
//   }
//   return L;
// }
// double loggaussian(const arma::colvec & s, double svar){ //
//   arma::vec x = log(1+s);
//   // double L = -.5*n*log(2*PI) -.5*n*log(svar) -(1/(2*svar))*sum(arma::pow((x-mean),2));
//   double L=0;
//   for(int i = 0; i<s.n_elem; i++){
//     L+= R::dnorm(x(i),0,svar,true);
//   }
//   return L;
// }
// double logmixgaussian(const arma::colvec & s, double svar, double ss){
//     arma::vec x = log(1+s);
//     double L=0;
//     for(int i = 0; i<s.n_elem; i++){
//       if(x(i)==0){
//         L += log(ss  + (1-ss) * R::dnorm(x(i),0,svar,false)) ;
//       }else{
//         L += (1-ss) * R::dnorm(x(i),0,svar,true);  // L+= R::dnorm(x(i),0,svar,true);
//       }
//     }
//     return L;
//   }
///// prior

// //[[Rcpp::export]]
// double logmixgaussian(const double & s, double svar, double ss){
//     double x = log(1+s);
//     double L;
//       if(s==0){
//         L= (ss  + (1-ss) * R::dnorm(s,0,svar,false)) ;
//       }else{
//         L= ((1-ss) * R::dnorm(s,0,svar,false));  // L+= R::dnorm(x(i),0,svar,true);
//       }
//     return log(L);
// }

/////  LIKELIHOOD
// [[Rcpp::export]]
double LLGaussMix(const double & y, const double & w, const double & v, const double & p){
  double LL;
  if(y==0){
    LL = p  + (1-p) *  R::pnorm(0,w,v,true,false) ;
  }else{
    LL = (1-p) * R::dnorm(y,w,v,false);
  }
return log(LL);
}

// [[Rcpp::export]]
double LIKELIHOOD(const arma::vec & y, // worked in previous iteration
                  const arma::vec & w,
                  const double & b, const double & a,
                  const double & p, const double & mu,
                  const double & epi,
                  bool verbose=false,
                  bool printall=false){
    double loglowest= -1e-300/y.n_elem;
    int countinf=0;
    double L=0;
    double LL;
      for(int i=0; i< y.n_elem ; i ++){ // check for infinity
        LL= LLGaussMix(y(i)/mu,w(i),w(i)*b+a,p);
        if(!std::isnan(LL)) LL=loglowest;
        if(!std::isinf(LL)) LL=loglowest; countinf++;
        L += LL; // if there is one nan, all sum is nan
      }
  cout << "# inf likelihoods " << countinf << endl;
  if(verbose) cout<< L << endl;
  if(std::isinf(L)) L=-1e-300;
  return(L);
}

// [[Rcpp::export]]
double likelihoodC(const Rcpp::NumericVector &  y,
                   const Rcpp::NumericVector & w,
                   const double & b,
                   const double & a,
                   const double & p){
  double mymin = -1022; // this is the minimum exponent value

  double LL=0;
  for(int i=0; i< y.size() ; i ++){ // check for infinity
    double L=0;
    if(y(i)==0){
      L= p  + (1-p) *  R::pnorm(0,w(i),a+w(i)*b,true,false) ;
    }else{
      L=(1-p) * R::dnorm(y(i),w(i),a+w(i)*b,false) ;
    }
    L=log(L);
    if(std::isnan(L)) L = mymin;
    if(std::isinf(L) & L<0) L = mymin;
    if(std::isinf(L) & L>0) L = -mymin;
    LL +=L;
  }
  if(std::isnan(LL)) LL = mymin;
  if(std::isinf(LL) & LL<0) LL = mymin;
  return(LL);
}

double likelihoodC(const arma::vec & y,
                   const arma::vec & w,
                   const double & b,
                   const double & a,
                   const double & p){
  double mymin = -1022; // this is the minimum exponent value

  double LL=0;
  for(int i=0; i< y.n_elem ; i ++){ // check for infinity
    double L=0;
    if(y(i)==0){
      L= p  + (1-p) *  R::pnorm(0,w(i),a+w(i)*b,true,false) ;
    }else{
      L=(1-p) * R::dnorm(y(i),w(i),a+w(i)*b,false) ;
    }
    L=log(L);
    if(std::isnan(L)) L = mymin;
    if(std::isinf(L) & L<0) L = mymin;
    if(std::isinf(L) & L>0) L = -mymin;
    LL +=L;
  }
  if(std::isnan(LL)) LL = mymin;
  if(std::isinf(LL) & LL<0) LL = mymin;
  return(LL);
}

////////////////////////////////////////////////////////////////////////////////
/// Fitness
////////////////////////////////////////////////////////////////////////////////

   // fitness mode
  double wmode(const double &s , const double &x, const int & FITmode){
    double w_;
    switch(FITmode){
      case 1:  // additive
        w_=  (s * x) ;
        break;
      case 2:  // multiplicative
        w_= 1 + (s * x) ;
        break;
      case 3:  // inverse multiplicative
        w_= pow((1 + s),x);
        break;
    }
    return(w_ );
  }
  // add fitness contribution of one SNP to overall fitness
  void wupdate(double &prevw, const double &s,const double &x, const int & FITmode) {     // void operator +*=(double w, double s,int x,int mode)  // operator not needed unless class
    switch(FITmode){
      case 1:  // additive
        prevw= prevw + wmode(s,x, FITmode);
        break;
      default:  // non-additive
        prevw= prevw * wmode(s,x, FITmode);
        break;
    }
  }
  // remove fitness contribution of one SNP to overall fitness (necessary to avoid repeating computations)
  void unwupdate(double &prevw, const double &s,const double &x, const int & FITmode) {
    switch(FITmode){
      case 1 :  // additive
        prevw= prevw - wmode(s,x, FITmode);
        break;
      default:  // non-additive
        prevw= prevw / wmode(s,x, FITmode);
        break;
    }
  }

// [[Rcpp::export]]
arma::vec wC(
               const arma::Mat<double> & X,
               const arma::vec & s,
               const int & mode,
               double epi=1,
               bool verbose=false){
    arma::vec w(X.n_rows); // w= new arma::vec(X.n_rows)
    w.fill(1); // need to fill otherwise unstable floats
    for (int i = 0; i < X.n_cols; i ++) {
        for(int j=0; j < X.n_rows ; j++){
          wupdate(w(j),s(i),X(j,i), mode );
        }
    }
  return(pow(w,epi));
}

// [[Rcpp::export]]
arma::vec wCBM(SEXP A, // worked in previous iteration
                   const arma::vec & s,
                   const arma::uvec & mycols,
                   const arma::uvec & myrows,
                   const int & mode,
                   double epi=1,
                   bool verbose=false){

  Rcpp::XPtr<BigMatrix> bigMat(A);
  MatrixAccessor<double> macc(*bigMat);
  arma::vec w(myrows.n_elem); // w= new arma::vec(X.n_rows)
    w.fill(1); // need to fill otherwise unstable floats
  int i, j;
  double val=0;
  for (j = 0; j <mycols.n_elem; j++) {
    for (i = 0; i < myrows.n_elem; i++) {
      val= (macc[mycols(j)-1][myrows(i)-1] );
      if(!std::isnan(val)) wupdate(w(i),s(j), val, mode );
    }
  }
   for (i = 0; i < myrows.n_elem; i++)  if(w(i)<0) w(i)=0; // CHECK BUGS
  return(pow(w,epi));
}


vector<double> sarr(){
  static const double arr[] ={
1.29116884119727e-02,2.25336598300678e-02,-6.76372108566958e-03,2.67219742814517e-03,-3.96657295675684e-03,-6.95333342601634e-04,-2.55639358064175e-03,1.41475436087206e-02,-7.28387835257605e-03,-9.54053067180183e-03,-1.29007645350121e-03,-4.18489923643839e-03,-2.41892757968837e-02,1.27798446262553e-02,1.30591639089239e-03,1.13324010984746e-02,-2.36743278014688e-02,7.51538764950821e-03,4.31745307097686e-03,-6.36353492617070e-03,-4.96984425104308e-03,4.26996092556298e-03,2.48653075663485e-03,-1.27175105053071e-02,4.04772297590306e-03,-1.07364946635902e-03,-5.76329962239952e-03,-2.36582456635870e-03,1.81280731188393e-03,6.01806606989030e-03,-1.42997225067236e-02,-1.46239518548508e-03,-6.51926182842322e-03,2.12863236838086e-02,-6.37589719098042e-03,3.10460474129060e-03,2.07416616161114e-02,-2.03627819254149e-03,-1.34435757222833e-02,-9.50185204945575e-03,-1.04815165358346e-02,1.79017860030688e-02,2.55571662425313e-03,-2.96722038214892e-03,6.30415546113028e-03,3.88304920021953e-03,-1.00807488679308e-02,4.86283431975298e-03,-2.03672663675765e-02,2.57447274714306e-03,2.89523371186884e-03,6.13356750555294e-03,-2.23182812015876e-03,4.40253955287284e-03,-3.70081680483225e-03,-1.66141398096411e-02,2.30947065783111e-02,-1.90663910754096e-03,4.16676494675605e-03,7.87671679529112e-03,2.33507865711147e-02,-1.90428341656178e-02,4.03455305816536e-03,-8.80749065708297e-03,-1.60717504745538e-02,3.79474464972440e-03,1.35093764785243e-03,6.09164448996680e-03,-1.48079026299666e-02,-1.02941630776154e-03,-1.47790797272853e-02,-6.09425573660516e-03,-4.97421762313188e-03,1.38309859548091e-02,-6.34442609254471e-03,1.21167802892441e-03,8.37140823879690e-03,-1.34802917207789e-02,8.37701362507492e-03,-1.04563291508208e-03,9.27000065548533e-03,-5.09173542861574e-03,-1.79618213631438e-04,2.79556368854394e-03,1.01445543657142e-03,5.56788378748996e-03,-1.94458627794111e-02,-1.88048442299915e-03,-9.19840342639422e-03,-5.80721937583717e-03,1.57553032473443e-02,2.24726169736920e-02,-5.50465036649284e-03,-4.46487288213537e-03,-4.01965759047285e-03,-1.68602468583190e-05,-8.61857748563810e-03,-3.82398636133430e-04,-6.91980025002448e-03,2.02543353822415e-03,-3.21562566783640e-03,-1.99994627424276e-02,4.12214521320231e-03,-9.19621972666262e-04,7.01462576140988e-03,1.32267451852743e-02,-7.45476469110262e-03,1.98882572524806e-03,-1.88492631304393e-03,-3.48590837763307e-03,-2.15630007181211e-03,2.68776821342076e-03,-6.57023402293788e-03,1.15187037173166e-02,1.45169551931656e-02,-1.68627733937078e-02,3.89789914876992e-03,2.13897646649697e-04,1.68829103295920e-02,-1.38791373204128e-02,8.97204727972367e-03,1.05483481199362e-02,2.29189226662463e-03,-1.27401580667859e-02,-3.17966554271265e-03,-2.57027766722984e-04,-4.88046134453834e-03,2.31247736033313e-03,5.50474126998823e-04,-1.74511638616576e-02,1.93970399240340e-03,-6.41216838840997e-04,5.38439373869615e-03,-1.15076963152687e-02,-1.08713623862325e-02,1.29871563804587e-02,4.10234072859850e-03,1.94567427660619e-03,-1.82378346886316e-02,2.39827643982649e-02,2.18037224340517e-02,-4.97310853943822e-03,2.69624488068820e-03,2.23010462381885e-02,2.61954253169563e-02,-5.12349090202291e-03,-9.11729099594127e-03,-6.67438927697017e-03,2.26955024857545e-03,-4.91203511350735e-03,-8.79674127623054e-03,-3.81777259808902e-03,1.11820472248236e-02,-1.40814603104303e-02,-3.23154377845292e-03,4.53277264786811e-03,1.48591451963822e-02,3.84809770841077e-03,-1.12369128484879e-03,-6.64722375081961e-03,6.49223347346339e-03,-8.13027448026116e-03,-9.96236186249910e-03,-6.20930631403271e-03,-8.69433018106469e-03,7.15815266885600e-03,1.06245777871590e-02,2.82209437059988e-03,1.29061984509704e-03,-4.25776984347426e-06,3.60063335558158e-03,-3.25636791844386e-03,-1.76051575637937e-04,1.10032388340389e-02,1.00836484400082e-02,2.06316168742426e-03,8.86692203087391e-03,-6.89957894381288e-03,-2.61803695017604e-03,1.70481989787556e-02,-2.26742623585574e-04,-1.80223819623616e-02,-5.90366438347889e-03,-5.34131097316459e-03,-1.16252533526453e-02,1.11904554252056e-02,-8.75250620428047e-03,1.73879183961989e-03,9.52362446188415e-03,-2.39797264583697e-03,4.70642582707681e-04,4.15351237971628e-03,-9.59109394266411e-03,-6.65160517161767e-03,-1.06747266860652e-03,1.57153976164204e-03,3.93280939331886e-03,-7.53004347156438e-03,1.57054591252572e-02,-9.33866686350837e-03,-4.79716420288190e-03,-2.60364966148263e-03,-7.31965692525405e-03,-1.27557298792552e-02,1.06600807898074e-02,2.04522338171431e-03,-2.57133716170117e-02,9.93382595696657e-03,1.89871106922026e-02,-7.21415348764176e-03,1.24503839996404e-04,-4.55797228171895e-03,4.07741845397047e-03,6.09625491712440e-03,-1.13624695432801e-03,2.06905387827110e-02,-1.13700883891625e-02,-9.07039756184691e-03,-2.76912951674497e-02,8.13903588302067e-03,1.81377888794014e-02,-1.60707861601223e-02,3.05344561885801e-03,1.78638413394649e-03,-1.83147822188529e-02,-7.72309847116659e-03,-6.09986868090995e-04,-7.49812215221213e-03,1.59227170389631e-04,-9.12836359999547e-03,-6.81560677401172e-03,-1.51997099534300e-02,-4.65576707476223e-03,-7.52089624071417e-03,3.01365297518585e-03,3.21105433638547e-03,5.90139565543701e-03,-1.47279802523792e-03,2.43972859792654e-02,-5.35966361405404e-04,4.54740068152959e-03,3.75903766788688e-03,2.21294133152883e-03,6.56890105905461e-03,-1.49897189464387e-03,1.80973921965688e-02,3.27701553405424e-03,-3.01407114740526e-03,-7.40040901441208e-03,5.32060405252754e-03,1.11351899769994e-02,-2.08262436917506e-03,-4.36980233287165e-03,4.35715379272072e-03,-7.63860656202942e-03,-3.84783301769376e-03,8.33914777583322e-03,8.87440947799467e-04,4.63330031825193e-03,3.81833852338231e-03,2.75085227155158e-03,4.62238901951517e-03,-1.92834116298303e-02,2.85087936863437e-03,-2.21033606522800e-02,-2.81148307596268e-03,-1.02333627414549e-02,1.03256582612203e-02,-2.28165545027570e-03,4.77353832524763e-03,-1.97049753801832e-02,9.29876925390705e-04,7.95102225722433e-04,-1.15010752578789e-02,6.22306335879852e-03,5.36824645830425e-03,-1.61891187907937e-03,-2.72254082555712e-04,-1.16644763258090e-02,2.80508424244297e-03,4.61577038394023e-03,1.50787323373258e-04,5.80964292883634e-03,-1.49366542179163e-03,2.39422983757025e-02,5.81473601600435e-03,-2.69054850751560e-03,-5.42627115336403e-03,-1.06302496685271e-02,-1.19120028031822e-02,4.80481133609190e-03,2.14777498098195e-03,2.07150774327736e-03,8.96947116485025e-04,5.79292785806595e-04,1.11004773088519e-03,-1.23244397112215e-02,2.02977223382250e-03,-1.13878170486824e-02,2.58051729283468e-02,-3.08001731643359e-04,-5.98135780055531e-03,9.82606628247451e-04,-4.28022884163792e-03,-5.80312241404179e-03,-5.02920935152351e-03,4.74093531371444e-03,2.22111081736824e-02,4.71285543499311e-03,-9.40753323117560e-03,-1.00683301491817e-02,5.93266976404383e-03,-1.11935050736940e-02,4.62718338560952e-03,-8.91246422221503e-03,-9.57035787876848e-03,-2.67700377183655e-03,8.52793658407491e-03,-1.33797062492182e-02,-9.70848101490385e-03,-1.94577413258532e-02,6.31713334199002e-03,4.13092948797966e-03,4.66358136934186e-04,3.91861507826907e-03,4.38528687917716e-03,-1.48148195546468e-02,-8.47858224819054e-03,-1.88463041341315e-02,-4.45861802051972e-05,-1.90602516576770e-02,5.13244209167785e-03,4.05991547678708e-03,1.10781599284953e-02,4.70801871050308e-03,9.67345925910523e-04,1.06324670484756e-02,-3.72465008094158e-03,-1.16353379521764e-02,-4.12370371213688e-03,4.26877122502445e-03,1.16139673724898e-02,1.98887278425721e-02,1.19286999213970e-02,1.35281787259680e-02,-2.36020083439293e-02,5.72669393884695e-03,-8.04127864246518e-03,1.53628039762099e-02,1.39819124241851e-02,7.06893017153276e-03,-1.45699062783435e-02,-9.54055108163387e-03,-5.23009614422731e-03,6.32996657286711e-03,-2.84660132249959e-03,5.20733544709606e-03,5.49362165279232e-03,3.30557783101026e-03,-9.67007511095874e-03,1.62274831171112e-02,9.56914883431192e-04,-1.27390770567604e-03,-1.26475772912634e-02,-1.09019122855780e-02,-3.03455928838670e-03,7.90961648952271e-03,-9.28383434635183e-03,-2.92869812993334e-02,-7.62997539570653e-04,-3.07668936416505e-03,6.20022589060998e-03,6.73464984866312e-03,3.75678508269495e-03,-7.26025215147452e-03,-7.21303172924437e-03,4.24841695407618e-03,-4.79303889460303e-03,1.88933180980704e-02,-4.64983277629238e-03,5.83614347983841e-03,1.90488515327660e-02,-9.21385505526717e-03,1.89587336136658e-02,9.59829452981809e-03,1.29774475701236e-02,6.91420526714026e-03,2.08248180988768e-03,-1.40351633016873e-02,6.00140995546972e-03,3.66801801005789e-03,1.10791816518312e-02,4.42746515955772e-04,3.54896609106703e-04,-1.35803655565025e-02,-7.27867284894157e-04,-1.42923731832489e-03,-1.43172434749643e-02,2.02425932774586e-02,2.44906759037100e-05,5.14789037485497e-03,-1.68021091389290e-02,-9.05369259917810e-03,-1.12214324240889e-03,-7.87209827679369e-03,-4.83399958044473e-05,1.51343553599972e-02,1.91400627678529e-03,1.10687258230993e-02,1.62373244087899e-03,-3.40440406796527e-03,4.55437114301915e-03,-8.66325607549878e-03,-2.95323081272014e-03,5.30788722989106e-03,-1.40315723289766e-02,-1.22901129212671e-02,9.17412157103414e-03,-5.27916277024132e-03,-5.32807844835248e-03,-1.35645003598819e-02,-1.22300939449820e-02,8.21556163302462e-03,3.78418606391451e-03,-6.86182115932521e-03,-1.11181721178318e-02,1.47814084815998e-02,1.68735840871086e-02,1.41125365947992e-03,-4.69748441852169e-03,-1.53305601414478e-02,-4.49775484790671e-03,7.20738371452279e-03,-1.44297401828252e-02,3.14857078448805e-03,9.51917109553535e-03,3.63236063177674e-03,-1.28269661226899e-03,-6.67075400267658e-03,8.80477914082012e-03,4.84653684251213e-03,-1.34468975599842e-02,3.45592297663844e-04,2.31050933190358e-03,8.07008143603460e-03,-4.90267405408884e-03,-2.31646625076687e-04,4.94610589370414e-03,-7.28943299710705e-03,3.82865797978549e-03,-7.88500054614150e-03,5.43039490448005e-03,-1.07894416288974e-02,1.18374505245575e-02,-1.44239226512026e-02,2.41146102328371e-03,3.57599920022245e-03,-1.04962027820681e-02,-1.74972635607873e-03,-1.67687492682391e-02,1.13885717588280e-03,1.91604021993936e-02,-1.58253410579565e-02,1.07142069625994e-02,-1.31090966042476e-02,-3.57488498619107e-03,-1.92085273826397e-02,4.11217209357173e-03,1.27897406831261e-02,-4.81009160747548e-03,-6.11983316650178e-03,-1.41697725179278e-02,2.30704429740776e-02,-4.87727507574731e-03,-1.35944271121525e-02,6.35423397563795e-03,-5.61806091260941e-03,1.98200949381522e-03,-7.97334669769478e-03,1.21531793330716e-02,-1.97475586715146e-03,-8.65117395011450e-03,-2.02945354851414e-03,-7.03141055712975e-03,-1.01201532039937e-02,1.56848812843533e-02,-8.47188985954794e-03,9.93346330493106e-03,-7.20737767509372e-03,-1.89145500222754e-04,-5.01376543794685e-03,5.78294070532670e-03,-1.85948087076682e-02,-7.31281098886127e-03,6.04244466641934e-03,1.82163318309843e-02,-1.53630987663109e-02,-9.26364548095693e-03,2.70522190095521e-02,-6.08359847240458e-03,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00,0.00000000000000e+00
  };
  // static const int arr[] = {16,2,77,29};
  vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
  return vec;
}

// [[Rcpp::export]]
arma::vec ssimC(arma::uvec snps,
                double svar){
  // arma::vec news;
  // news.load("../databig/s_svar01.txt");
  arma::vec news(sarr());

  news= exp( log(1+news) * (svar/0.01) )-1 ;

  arma::vec subnews(snps.n_elem,arma::fill::zeros);
  for( int i=0; i<snps.n_elem;i++){
    subnews(i) = news(snps(i)-1);
  }
  // if(svar==0.05){
  //   news.load("databig/s_svar05.txt");
  // }else if(svar==0.25){
  //   news.load("databig/s_svar25.txt");
  // }else if(svar==0.5){
  //   news.load("databig/s_svar50.txt");
  // } else{
  //   cout << "S values not found for accuracy calculations!" << endl;
  //   news=exp(Rcpp::rnorm(nsnp,0,svar))-1;
  // }
  return subnews;
}


////////////////////////////////////////////////////////////////////////////////
/// Spatial Projected Gradient
////////////////////////////////////////////////////////////////////////////////
// Modified from E Birgin software Tango
/* =================================================================
Module: Spectral Projected Gradient. Problem definition.
=================================================================

Last update of any of the component of this module:

March 14, 2008.

Users are encouraged to download periodically updated versions of
this code at the TANGO Project web page:

www.ime.usp.br/~egbirgin/tango/

=================================================================
================================================================= */

#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))

extern double *l,*u;

#define gamma 1.0e-04

#define lmax  1.0e+30
#define lmin  1.0e-30
#define M     100

double ptr(double i){ return 1.0+(-log((1.0-i)/(i))) ;}
double iptr(double i){ return (1.0/(1.0+exp(1.0-i))) ;}
double iptr05(double i){ return (0.5/(1.0+exp(1.0-i))) ;}
double iptr2(double i){ return (2/(1.0+exp(1.0-i))) ;}
double seltr(double i){ return 1.0+ log(1.0+i) ;}
double iseltr(double i){ return exp(i-1.0) -1.0 ;}



/***********************************************************************
 **********************************************************************/

void reperr(int inform) {

   char str[80];

   strcpy(str,"*** There was an error in the user supplied ");

   if ( inform == -90 )
      strcat(str,"EVALF");

   else if ( inform == -91 )
      strcat(str,"EVALG");

   else if ( inform == -92 )
      strcat(str,"PROJ");

   strcat(str," subroutine ***");

   printf("\n\n%s\n\n",str);
}

/***********************************************************************
 **********************************************************************/

void evalf(int n,double *x,
           const arma::Mat<double> &X,
           const arma::vec &y,
           double epi,
           double mod,
           double *f,int *flag) {
  int i;
  double b,a,p,mu;
  arma::vec sel(n-4,arma::fill::zeros);
  arma::vec w;

  *flag = 0;

  // Extract and transform parameters
  b=iptr(x[n-4]);
  a=iptr(x[n-3]);
  p=iptr05(x[n-2]);
  mu=iptr2(x[n-1]);
  for ( i=0; i<sel.n_elem; i++ ){
    sel[i] = iseltr(x[i]);
  }
   w = wC(X, sel,mod,epi);
  *f = 0.0;
  *f = (-likelihoodC(y/mu,w,b,a,p));

}

void sevalf(int n,double *x,
           const arma::Mat<double> &X,
           const arma::vec &y,
           double epi,
           double mod,
           double *f,int *inform){

   int flag;
  // evalf(n,x,f,&flag);
  evalf(n,x,X,y,epi,mod,f,&flag); // use pointer but to fsmall

   /* This is true if f if Inf, - Inf or NaN */ // I have already done checks
   // if ( ! ( *f > - 1.0e+99 ) || ! ( *f < 1.0e+99 ) )
   //    *f = 1.0e+99;

   if ( flag != 0 ) {
      *inform = -90;
      reperr(*inform);
   }
}

/***********************************************************************
 **********************************************************************/
void evalg(int n,double *x,
           const arma::Mat<double> &X,
           const arma::vec &y,
           const double epi,
           const double mod,
           double *f,int *flag,
           double *g,
           const double eps){
   int i;

  *flag = 0;
   double fsmall=0.0;
   const double fref= *f;

   for ( i=0; i<n; i++ ){
      x[i] = x[i] + eps; // add stepsize
      evalf(n,x,X,y,epi,mod,&fsmall,flag); // use pointer but to fsmall
      g[i] = (fsmall- fref ) / eps;
      x[i] = x[i] - eps; // remove stepsize
   }
}

void sevalg (int n,double *x,
             const arma::Mat<double> &X,
             const arma::vec &y,
             const double epi,
             const double mod,
             double *f, int *inform,
             double *g,
             const double eps){
   int flag;

   evalg(n,x,X,y,epi,mod,f,&flag,g,eps); //  update with custom numeric gradient

   if ( flag != 0 ) {
      *inform = -91;
      reperr(*inform);
   }
}

/***********************************************************************
 **********************************************************************/
// void proj(int n,double *x,int *flag) {
//    int i;
//
//   *flag = 0;
//
//    for ( i=0; i<n; i++ ) // I do not want box constrints!
//       x[i] = max(l[i], min(x[i], u[i]));
// }
//
//
// void sproj(int n,double *x,int *inform) {
//
//    int flag;
//
//    proj(n,x,&flag);
//
//    if ( flag != 0 ) {
//       *inform = -92;
//       reperr(*inform);
//    }
// }


/***********************************************************************
 **********************************************************************/

void ls(int n,double *x,
        const arma::Mat<double> &X,
        const arma::vec &y,
        double epi,
        double mod,
        double f,double *g,double *d,double *lastfv,
        int maxfc,int *fcnt,double *fnew,double *xnew,int *lsinfo,
        int *inform) {

/* Nonmonotone line search with safeguarded quadratic interpolation

   lsinfo:

   0: Armijo-like criterion satisfied
   2: Maximum number of functional evaluations reached */

   int i;
   double alpha,atmp,fmax,gtd;

   fmax = lastfv[0];
   for ( i=1; i<M; i++ ) fmax = max (fmax, lastfv[i]);

   gtd = 0.0;
   for ( i=0; i<n; i++ ) gtd+= g[i] * d[i];

   alpha = 1.0;
   for ( i=0; i<n; i++ ) xnew[i] = x[i] + alpha * d[i];

   sevalf(n,xnew,
           X,y,epi,mod,
           fnew,inform); (*fcnt)++;
   if ( *inform != 0 ) return;

///// Main loop /////
   while ( *fnew > fmax + gamma * alpha * gtd && *fcnt < maxfc ) { // I think here is the problem
      // Safeguarded quadratic interpolation
      if ( alpha <= 0.1 ){
        alpha /= 2.0;
      }else {
        atmp = ( - gtd * alpha * alpha ) /
  	         ( 2.0 * ( *fnew - f - alpha * gtd ) );
         if ( atmp < 0.1 || atmp > 0.9 * alpha ){
            atmp = alpha / 2.0;
         }
       alpha = atmp;
      }

      // New trial
      for ( i=0; i<n; i++ ) xnew[i] = x[i] + alpha * d[i];

      sevalf(n,xnew,
             X,y,epi,mod,
             fnew,inform); (*fcnt)++;
      if ( *inform != 0 ) return;
  }
///// End of main loop /////

/* Termination flag */
   if ( *fnew <= fmax + gamma * alpha * gtd )
      *lsinfo = 0;
   else if ( *fcnt >= maxfc )
      *lsinfo = 2;
}



/***********************************************************************
 **********************************************************************/

// void spg(int n,double *x,double epsopt,int maxit,int maxfc,int iprint,
// 	 double *f,double *gpsupn,int *iter,int *fcnt,int *spginfo,
// 	 int *inform) {

void spg(int n,double *x,
         const arma::Mat<double> &X,
         const arma::vec &y,
         double epi,
         double mod,
         double epsopt,int maxit,int maxfc,int iprint,
      	 double *f,double *gpsupn,int *iter,int *fcnt,int *spginfo,
      	 int *inform,
      	 double eps) {

/* Subroutine SPG implements the Spectral Projected Gradient Method
   (Version 2: "Feasible continuous projected path") to find a
   local minimizers of a given function with convex constraints,
   described in

   E.G. Birgin, J.M. Martinez and M. Raydan, "Nonmonotone spectral
   projected gradient methods for convex sets", SIAM Journal on
   Optimization 10, pp. 1196-1211, 2000.

   The user must supply the external subroutines evalf, evalg and
   proj to evaluate the objective function and its gradient and to
   project an arbitrary point onto the feasible region.

   This version 14 MARCH 2008 by E.G.Birgin, J.M.Martinez and M.Raydan.

   Other parameters (i means input, o means output):

   n       (i)   number of variables
   x       (i/o) initial guess on input, solution on output
   epsopt  (i)   tolerance for the convergence criterion
   maxit   (i)   maximum number of iterations
   maxfc   (i)   maximum number of functional evaluations
   iprint  (i)   controls output level (0 = no print)
   f       (o)   functional value at the solution
   gpsupn  (o)   sup-norm of the projected gradient at the solution
   iter    (o)   number of iterations
   fcnt    (o)   number of functional evaluations
   spginfo (o)   indicates the reason for stopping
   inform  (o)   indicates an error in an user supplied subroutine

   spginfo:

   0: Small continuous-projected-gradient norm
   1: Maximum number of iterations reached
   2: Maximum number of functional evaluations reached

   spginfo remains unset if inform is not equal to zero on output

   inform:

     0: ok
   -90: error in the user supplied evalf subroutine
   -91: error in the user supplied evalg subroutine
   -92: error in the user supplied proj  subroutine */

   int i,lsinfo;
   double fbest,fnew,lambda,sts,yty,sty;
   double *d,*g,*gnew,*gp,*lastfv,*s,*xbest,*xnew,*j;
   FILE *tabline;

   char presentation[]=
  "\n============================================================================\n"
     " Non-Additive Polygenic model (NAP) for inference of genome-wide        \n"
     " selection coefficients under additive and non-additive fitness. \n"
     "============================================================================\n\n";

   //// Initialization


// Print problem information

   if ( iprint > 0 ) {
      printf("%s",presentation);
      printf(" Entry to SPG optimization\n");
      printf(" Number of Variables: %d\n",n);
   }

// Get memory
   lastfv = (double *) malloc (M * sizeof(double));
   d      = (double *) malloc (n * sizeof(double));
   g      = (double *) malloc (n * sizeof(double));
   gnew   = (double *) malloc (n * sizeof(double));
   gp     = (double *) malloc (n * sizeof(double));
   s      = (double *) malloc (n * sizeof(double));
   xbest  = (double *) malloc (n * sizeof(double));
   xnew   = (double *) malloc (n * sizeof(double));
   j      = (double *) malloc (n * sizeof(double));

//// Set some initial values:
// error tracker
   *inform = 0;
// for counting number of iterations as well as functional evaluations
   *iter = 0;
   *fcnt = 0;
// for the non-monotone line search //
   for ( i=0; i<M; i++ ) lastfv[i] = -1.0e+99;

// Project initial guess
  // if ( iprint > 0 ) printf(" Project initial guess\n");
  //  sproj(n,x,inform);
  //  if ( *inform  != 0 ) return;

// Compute function and gradient at the initial point
   sevalf(n,x,
          X,y,epi,mod,
          f,inform);  (*fcnt)++ ;
   if ( iprint > 0 ) printf(" Function value at initial point %e\n",*f);
   if ( *inform != 0 ) return;
   sevalg (n,x,
           X,y,epi,mod,
           f, inform,
           g,eps);
   if ( iprint > 0 )printf(" Compute gradient at initial point %e\n",*g);
   if ( *inform != 0 ) return;

// Store functional value for the non-monotone line search
   lastfv[0] = *f;

// Compute continuous-project-gradient and its sup-norm
   for ( i=0; i<n; i++ )
      gp[i] = x[i] - g[i];

   // sproj(n,gp,inform);
   // if ( *inform != 0 ) return;

   *gpsupn = 0.0;
   for ( i=0; i<n; i++ ) {
      gp[i] -= x[i];
      *gpsupn = max( *gpsupn, fabs( gp[i] ) );
   }

// Initial steplength //
   if ( *gpsupn != 0.0 )
      lambda = min( lmax, max( lmin, 1.0 / *gpsupn ) );
   else
      lambda = 0.0;

// Initiate best solution and functional value
   fbest = *f;

   for ( i=0; i<n; i++ )
      xbest[i] = x[i];

// Print initial information //
   if ( iprint >  0 ) {
       printf("\n ITER\t F\t\t GPSUPNORM\n");
      printf(" %d\t %e\t %e\n",*iter,*f,*gpsupn);
   }

/* ==================================================================
   Main loop
   ================================================================== */

   // while( *gpsupn > epsopt && *iter < maxit && *fcnt < maxfc ) {
   while( *gpsupn > epsopt && *iter < maxit && *fcnt < maxfc ) {
    // in R there is also ftol, the change between successive f iterations

      // Iteration
      (*iter)++;

      // Compute search direction
      for ( i=0; i<n; i++ ) d[i] = x[i] - lambda * g[i];

      // sproj(n,d,inform);
      // if ( *inform != 0 ) return;

      for ( i=0; i<n; i++ ) d[i]-= x[i];

      // Perform safeguarded quadratic interpolation along the
      //spectral continuous projected gradient

      ls(n,x,X,y,epi,mod,*f,g,d,lastfv,maxfc,fcnt,&fnew,xnew,&lsinfo,inform);
      if ( *inform != 0 ) return;

      // Set new functional value and save it for the non-monotone line search
      *f = fnew;
      lastfv[(*iter) % M] = *f;

      // Gradient at the new iterate
      sevalg (n,xnew,
           X,y,epi,mod,
           f, inform,
           gnew,eps);
      if ( *inform != 0 ) return;

      // // bug check
      // for( i=0; i<5; i++ ) std::cout << "gnew" << gnew[i] << std::endl;
      // // bug check

      // Compute s = xnew - x and y = gnew - g, <s,s>, <s,y>, the
      //   continuous-projected-gradient and its sup-norm
      sts = 0.0;
      sty = 0.0;
      yty = 0.0;
      for ( i=0; i<n; i++ ) {
        s[i]  = xnew[i] - x[i];
        j[i]  = gnew[i] - g[i];
        sts  += s[i] * s[i];
        yty  += j[i] * j[i];
        sty  += s[i] * j[i];
        x[i]  = xnew[i];
        g[i]  = gnew[i];
        gp[i] = x[i] - g[i];
      }

      // sproj(n,gp,inform);
      // if ( *inform != 0 ) return;

      *gpsupn = 0.0;
      for ( i=0; i<n; i++ ) {
        gp[i] -= x[i];
        *gpsupn = max( *gpsupn, fabs( gp[i] ) );
      }

      // Spectral steplength // Different in R version
        /*  if (method == 1)
              lambda <- if (sts == 0 | sty < 0)
                  lmax
              else min(lmax, max(lmin, sts/sty))
          else if (method == 2)
              lambda <- if (sty < 0 | yty == 0)
                  lmax
              else min(lmax, max(lmin, sty/yty))
          else if (method == 3) # DEFAULT
              lambda <- if (sts == 0 | yty == 0)
                  lmax
              else min(lmax, max(lmin, sqrt(sts/yty))) */
      // Original
      /*if ( sty <= 0 ) {
        lambda = lmax;
      }else{
        lambda = max( lmin, min( sts / sty, lmax ) );
      }*/
      // Trying to implement mehtod 3
      if (sts == 0 || yty == 0){
        lambda = lmax;
      }else{
        lambda = min(lmax, max(lmin, sqrt(sts/yty)));
      }

      // Best solution and functional value
      if ( *f < fbest ) {
        fbest = *f;
        for ( i=0; i<n; i++ )xbest[i] = x[i];
      }

      // Print iteration information
      if ( iprint >  0 ) {
        if ( (*iter) % 10 == 0 ){
          printf(" %d\t %e\t %e\n",*iter,*f,*gpsupn);
        }
      }
   }
/* ==================================================================
   End of main loop
   ================================================================== */

/* Finish returning the best point */

   *f = fbest;

   for ( i=0; i<n; i++ )
      x[i] = xbest[i];

/* Write statistics */

   if ( iprint > 0 ) {
      printf("\n");
      printf (" Number of iterations               : %d\n",*iter);
      printf (" Number of functional evaluations   : %d\n",*fcnt);
      printf (" Objective function value           : %e\n",*f);
      printf (" Sup-norm of the projected gradient : %e\n",*gpsupn);
   }

/* Free memory */
   free(lastfv);
   free(d);
   free(g);
   free(gnew);
   free(gp);
   free(s);
   free(xbest);
   free(xnew);
   free(j);

/* Termination flag */

   if ( *gpsupn <= epsopt ) {
      *spginfo = 0;
      if ( iprint > 0 )
	 printf("\n Flag of SPG: Solution was found.\n");
   }

   else if ( *iter >= maxit ) {
      *spginfo = 1;
      if ( iprint > 0 )
	 printf("\n Flag of SPG: Maximum of iterations reached.\n");
   }

   else {
      *spginfo = 2;
      if ( iprint > 0 )
	 printf("\n Flag of SPG: Maximum of functional evaluations reached.\n");
   }
}


////////////////////////////////////////////////////////////////////////////////
/* Main program */
// [[Rcpp::export]]
Rcpp::List napspgC(/* input data */
            std::string bedfile,
            int N, int p,
            const arma::uvec & myrows,
            const arma::uvec & mycols,
            arma::vec y,
            arma::vec par,
            double epi,
            double mod,
            int maxit=20,
            bool verbose=1
            ){
   ///////////////////////////////////////////////
   // std::cout << "Initializing objects" <<std::endl;
   int fcnt,iprint,i,inform,iter,maxfc,n,spginfo;
   double epsopt,f,gpsupn,*x;
   double eps;
   n = par.n_elem ;
   iprint = verbose;
   maxfc  = 10000;
   epsopt = 1.0e-05;
   eps= 1.0e-07;

   // Get memory
   x = (double *) malloc ( n * sizeof(double));
   // l = (double *) malloc ( n * sizeof(double));
   // u = (double *) malloc ( n * sizeof(double));
   // Define bounds
   // for ( i=0; i<n; i++ ) {
   //    l[i] = -Inf;
   //    u[i] =  Inf;
   // }

  // Define initial Guess //
   for ( i=0; i<n; i++ )
      x[i] = par[i];
   ///////////////////////////////////////////////
   // Get genome matrix
    if(verbose) std::cout << " Reading .bed file (missing data as reference)" <<std::endl;
   //* load matrix */
    arma::Mat<double> X = readbed(bedfile,N,p,myrows, mycols);
    std::cout << "Analysing " <<  X.n_rows << " individuals and " <<
              X.n_cols << " variants" <<std::endl;
   ///////////////////////////////////////////////
   // Call SPG //
   spg(n,x,X,y,epi,mod,epsopt,maxit,maxfc,iprint,&f,&gpsupn,&iter,&fcnt,&spginfo,&inform,eps);
    // std::cout << "Finished" <<std::endl;

  ///////////////////////////////////////////////
  // get solutions
   std::cout << "Calculating solutions" <<std::endl;
   arma::uvec allrows(N,arma::fill::zeros);
   X = readbed(bedfile,N,p,allrows, mycols); // dummy variable so no subset

  double b,a,pi,mu;
  arma::vec sel(n-4,arma::fill::zeros);
  arma::vec w;

  b=iptr(x[n-4]);
  a=iptr(x[n-3]);
  pi=iptr05(x[n-2]);
  mu=iptr2(x[n-1]);
  for ( i=0; i<sel.n_elem; i++ ) sel[i] = iseltr(x[i]);
  w = mu * wC(X, sel,mod,epi);

  ///////////////////////////////////////////////
   // Free memory //
  // get vector of fitness
	return Rcpp::List::create(
	                    Rcpp::Named("s") = sel,
	                    Rcpp::Named("a") = a,
	                    Rcpp::Named("b") = b,
	                    Rcpp::Named("p") = pi,
	                    Rcpp::Named("mu") = mu,
	                    Rcpp::Named("w") = w,
	          					Rcpp::Named("f") = f);
}


// [[Rcpp::export]]
arma::vec polyscore(std::string bedfile,
                    int N, int p,
                    const arma::vec betas,
                    const arma::uvec & mycols,
                    double mu=1,
                    double mod=1,
                    double epi=1
                    ){
  arma::Mat<double> X(N,p);
  arma::uvec allrows(N) ;
  arma::vec w;
  X=readbed(bedfile,N,p,allrows,mycols);
  w = mu * wC(X, betas,mod,epi);
  return w;
}
