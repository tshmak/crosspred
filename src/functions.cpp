/**
   lassosum
   functions.cpp
   Purpose: functions to perform p-thresholding

   @author Timothy Mak

   @version 0.1

 */
// [[Rcpp::interfaces(r, cpp)]]

#include <stdio.h>
#include <string>
#include <bitset>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
#include <lassosum.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace lassosum;

/**
 	Opens a Plink binary files

	@s file name
	@BIT ifstream
	@return is plink file in major mode

*/

bool openPlinkBinaryFile(const std::string s, std::ifstream &BIT) {
  BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  if (!BIT.is_open()) {
    throw "Cannot open the bed file";
  }

  // 2) else check for 0.99 SNP/Ind coding
  // 3) else print warning that file is too old
  char ch[1];
  BIT.read(ch, 1);
  std::bitset<8> b;
  b = ch[0];
  bool bfile_SNP_major = false;
  bool v1_bfile = true;
  // If v1.00 file format
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  // std::cerr << "check magic number" << std::endl;
  if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
    // Next number
    BIT.read(ch, 1);
    b = ch[0];
    if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
      // Read SNP/Ind major coding
      BIT.read(ch, 1);
      b = ch[0];
      if (b[0])
        bfile_SNP_major = true;
      else
        bfile_SNP_major = false;

      // if (bfile_SNP_major) std::cerr << "Detected that binary PED file is
      // v1.00 SNP-major mode" << std::endl;
      // else std::cerr << "Detected that binary PED file is v1.00
      // individual-major mode" << std::endl;

    } else
      v1_bfile = false;

  } else
    v1_bfile = false;
  // Reset file if < v1
  if (!v1_bfile) {
    Rcerr << "Warning, old BED file <v1.00 : will try to recover..."
              << std::endl;
    Rcerr << "  but you should --make-bed from PED )" << std::endl;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
    BIT.read(ch, 1);
    b = ch[0];
  }
  // If 0.99 file format
  if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
    Rcerr << std::endl
              << " *** Possible problem: guessing that BED is < v0.99      *** "
              << std::endl;
    Rcerr << " *** High chance of data corruption, spurious results    *** "
              << std::endl;
    Rcerr
        << " *** Unless you are _sure_ this really is an old BED file *** "
        << std::endl;
    Rcerr << " *** you should recreate PED -> BED                      *** "
              << std::endl
              << std::endl;
    bfile_SNP_major = false;
    BIT.close();
    BIT.clear();
    BIT.open(s.c_str(), std::ios::in | std::ios::binary);
  } else if (!v1_bfile) {
    if (b[0])
      bfile_SNP_major = true;
    else
      bfile_SNP_major = false;
    Rcerr << "Binary PED file is v0.99" << std::endl;
    if (bfile_SNP_major)
      Rcerr << "Detected that binary PED file is in SNP-major mode"
                << std::endl;
    else
      Rcerr << "Detected that binary PED file is in individual-major mode"
                << std::endl;
  }
  return bfile_SNP_major;
}

//' Mutiply a bed file with a vector with p-value thresholds bins
// [[Rcpp::export]]
arma::mat multiBed4(const std::string fileName, int N, int P,
                    const arma::vec weights, arma::Col<int> pbin, int nbin,
                    const arma::Col<int> col_skip_pos, const arma::Col<int> col_skip,
                    const arma::Col<int> keepbytes, const arma::Col<int> keepoffset, 
                    const int trace) {

  // Similar to multiBed3 but only for p-value thresholding

  std::ifstream bedFile;
  bool snpMajor = openPlinkBinaryFile(fileName, bedFile);
  if (!snpMajor)
    throw std::runtime_error("We currently have no plans of implementing the "
                               "individual-major mode. Please use the snp-major "
                               "format");

  int i = 0;
  int ii = 0;
  int iii = 0;
  const bool colskip = (col_skip_pos.n_elem > 0);
  unsigned long long int Nbytes = ceil(N / 4.0);
  const bool selectrow = (keepbytes.n_elem > 0);
  int n;
  if (selectrow)
    n = keepbytes.n_elem;
  else
    n = N;
  int jj;

  arma::mat result = arma::mat(n, nbin, arma::fill::zeros);
  // pbin.transform( [](int x) {return (x - 1); }); // Change from one-based to zero-based vector

  std::bitset<8> b; // Initiate the bit array
  char ch[Nbytes];

  int chunk;
  double step;
  double Step = 0; 
  if(trace > 0) {
    chunk = weights.n_elem / pow(10,trace); 
    step = 100 / pow(10,trace); 
    // Rcout << chunk << " " << step << " " << trace << " " << (10^trace) << "\n"; 
    // Rcout << "Started C++ program \n"; 
  }
  
  while (i < P) {
    Rcpp::checkUserInterrupt();
    if (colskip) {
      if (ii < col_skip.n_elem) {
        if (i == col_skip_pos[ii]) {
          bedFile.seekg(col_skip[ii] * Nbytes, bedFile.cur);
          i = i + col_skip[ii];
          ii++;
          continue;
        }
      }
    }
    
    if(trace > 0) {
      if (iii % chunk == 0) {
        Rcout << Step << "% done\n";
        Step = Step + step; 
      }
    }
    
    bedFile.read(ch, Nbytes); // Read the information
    if (!bedFile)
      throw std::runtime_error(
          "Problem with the BED file...has the FAM/BIM file been changed?");

    int j = 0;
    if (!selectrow) {
      for (jj = 0; jj < Nbytes; jj++) {
        b = ch[jj];

        int c = 0;
        while (c < 7 && j < N) { // from the original PLINK: 7 because of 8 bits
          int first = b[c++];
          int second = b[c++];
          if (first == 0) {
            if (weights[iii] != 0.0) {
              result(j, pbin[iii]) += (2 - second) * weights[iii];
            }
          }
          j++;
        }
      }
    } else {
      for (jj = 0; jj < keepbytes.n_elem; jj++) {
        b = ch[keepbytes[jj]];

        int c = keepoffset[jj];
        int first = b[c++];
        int second = b[c];
        if (first == 0) {
          if (weights[iii] != 0.0) {
            result(j, pbin[iii]) += (2 - second) * weights[iii];
          }
        }
        j++;
      }
    }

    i++;
    iii++;
  }

  arma::mat result2 = cumsum(result, 1);
  return result2;
}

//' An overall beta for cross-prediction 
//' 
//' @param fileName location of bam file
//' @param N number of subjects 
//' @param P number of positions 
//' @param col_skip_pos which variants should we skip
//' @param col_skip which variants should we skip
//' @param keepbytes which bytes to keep
//' @param keepoffset what is the offset
//' @param pred The expected prediction
//' @param meanbeta mean beta in the cross-prediction
//' @return an armadillo vector
//' @keywords internal
//' 
// [[Rcpp::export]]
arma::vec overallbeta(const std::string fileName, int N, int P,
                      arma::Col<int> col_skip_pos, arma::Col<int> col_skip, 
                      arma::Col<int> keepbytes, arma::Col<int> keepoffset, // const int fillmissing
                      arma::vec pred, arma::vec meanbeta, 
                      const std::string save = "", const std::string load = "") {
  
  arma::mat X = genotypeMatrix(fileName, N, P, col_skip_pos, col_skip,
                               keepbytes, keepoffset, 1);

// arma::mat X = arma::mat(N, P, arma::fill::zeros);
  
  arma::mat inv; 
  if(load == "") {
    inv = arma::pinv(X);
    if(save != "") inv.save(save); 
  } 
  else {
    inv.load(load);
  }
  
  arma::vec Xmeanbeta = X * meanbeta; 
  arma::vec right = inv * (Xmeanbeta -  pred); 
  arma::vec result = meanbeta - right; 
  
  return result;
}
