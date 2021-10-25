
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// factorial
int fact(int n){
  int prod = 1;
  for(int i = n; i > 0; i--){
    prod *= i;
  }
  return prod;
}

int choose(int n, int k){
  return fact(n) / (fact(k) * fact(n - k));
}

IntegerMatrix combnIndices2(int n){
  IntegerMatrix mat(choose(n, 2), 2);
  int counter = 0;
  for(int i = 0; i < (n - 1); ++i){
    for(int j = i + 1; j < n; ++j){
      mat(counter, 0) = i;
      mat(counter, 1) = j;
      counter++;
    }
  }
  return mat;
}


NumericVector getAbc(const NumericMatrix &count_mat, const NumericMatrix &het_mat, const IntegerVector &n_i, const IntegerVector &cols){

  NumericVector n_i2(cols.length());
  for(int i = 0; i < cols.length(); i++){
    n_i2[i] = n_i[cols[i]];
  }
  
  double r = (double) n_i2.length(); // the number of populations
  double n_bar = sum(n_i2) / r;
  double n_c = ((r * n_bar) - (sum(pow(n_i2,2)) / (r * n_bar))) / (r - 1);

  NumericVector p_bar(count_mat.nrow());
  NumericVector s2(count_mat.nrow());
  for(int j = 0; j < count_mat.nrow(); ++j){
    NumericVector p_i(cols.length());
    NumericVector p_bar_numerator(cols.length());
    for(int k = 0; k < cols.length(); ++k){
      p_i[k] = count_mat(j,cols[k]) / (2 * n_i2[k]);  
      p_bar_numerator[k] = p_i[k] * n_i2[k];
    }
    p_bar[j] = sum(p_bar_numerator) / (r * n_bar);
    NumericVector s2_numerator(p_i.length());
    for(int k = 0; k < p_i.length(); ++k){
      s2_numerator[k] = pow(p_i[k] - p_bar[j], 2) * n_i2[k];  
    }
    s2[j] = sum(s2_numerator) / ((r - 1) * n_bar);
  }
  
  if(r == 1){
    s2 = NumericVector(count_mat.nrow()); // if there's only one population set 's2' to be all 0s
  }
  NumericVector h_bar(het_mat.nrow());
  for(int j = 0; j < het_mat.nrow(); ++j){
    double sum = 0;
    for(int k = 0; k < cols.length(); ++k){
      sum += het_mat(j,cols[k]);
    }
    h_bar[j] = sum / (r * n_bar);
  }
  // Rcout << "r: " << r << "\n";
  // Rcout << "n_bar: " << n_bar << "\n";
  // Rcout << "n_c: " << n_c << "\n";
  // Rcout << "p_bar: " << p_bar << "\n";
  // Rcout << "s2: " << s2 << "\n";
  // Rcout << "h_bar: " << h_bar << "\n";
  NumericVector a = (n_bar / n_c) * (s2 - (1 / (n_bar - 1)) * ((p_bar * (1 - p_bar)) - (((r - 1) / r) * s2) - (h_bar / 4)));
  NumericVector b = (n_bar / (n_bar - 1)) * ((p_bar * (1 - p_bar)) - (((r - 1) / r) * s2)  - ((2 * n_bar - 1) / (4 * n_bar)) * h_bar);
  NumericVector c = h_bar / 2;
  double a_sum = sum(a);
  double b_sum = sum(b);
  double c_sum = sum(c);  
  
  NumericVector abc = NumericVector::create(a_sum, b_sum, c_sum);
  return abc;
}

NumericVector getFstatsFromAbc(const NumericMatrix &abc_mat){
  double fit = 1 - sum(abc_mat(_,2)) / sum(abc_mat);
  double fst = sum(abc_mat(_,0)) / sum(abc_mat);
  double fis = 1 - sum(abc_mat(_,2)) / (sum(abc_mat(_,1)) + sum(abc_mat(_,2)));
  NumericVector fstats = NumericVector::create(fit, fst, fis);
  fstats.names() = CharacterVector({"Fit", "Fst", "Fis"});
  return fstats;
}

// [[Rcpp::export]]
List wcCpp(const IntegerMatrix als, const IntegerVector pop) {
  // get the size of each population
  IntegerVector n_i = table(pop);
  
  IntegerVector pop_ind = seq(0,n_i.length()-1);
  pop_ind.names() = n_i.names();

  NumericMatrix abc_mat(als.ncol()/2, 3);
  IntegerMatrix pairs = combnIndices2(n_i.length());
  std::vector<NumericMatrix> pw_abc_mats(pairs.nrow());
  for(int i = 0; i < pairs.nrow(); ++i){
    pw_abc_mats[i] = NumericMatrix(als.ncol()/2, 3);
  }
  for(int i = 0; i < (als.ncol() - 1); i += 2){
    IntegerVector al_unq1 = unique(als(_, i));
    IntegerVector al_unq2 = unique(als(_, i+1));
    IntegerVector al_unq = union_(al_unq1, al_unq2);
    IntegerVector al_ind = seq(0,al_unq.length()-1);
    
    al_ind.names() = al_unq;
    int al_n = al_unq.length();
    NumericMatrix count_mat(al_n, n_i.length());
    NumericMatrix het_mat(al_n, n_i.length());
    
    for(int rw = 0; rw < als.nrow(); rw++) {
      std::string pop_i = std::to_string(pop[rw]);
      std::string al1 = std::to_string(als(rw, i));
      std::string al2 = std::to_string(als(rw, i + 1));
      if(al1 != al2){
        het_mat(al_ind[al1], pop_ind[pop_i])++;
        het_mat(al_ind[al2], pop_ind[pop_i])++;
      }
      for(int cn = i; cn <= i+1; cn++){
        std::string al_i = std::to_string(als(rw, cn));
        count_mat(al_ind[al_i], pop_ind[pop_i])++;
      }
    }
    
    // Rcout << "count_mat\n" << count_mat << "\n\n";
    // Rcout << "het_mat\n" << het_mat << "\n\n";
    abc_mat(i/2,_) = getAbc(count_mat, het_mat, n_i, seq(0,count_mat.ncol()-1));
    for(int j = 0; j < pairs.nrow(); ++j){
      IntegerVector cols = IntegerVector::create(pairs(j,0), pairs(j,1));
      pw_abc_mats[j](i/2,_) = getAbc(count_mat, het_mat, n_i, cols);
    }
  }
  
  NumericMatrix f_loc(als.ncol()/2, 3);
  colnames(f_loc) = Rcpp::CharacterVector({"Fit","Fst","Fis"}); //name the columns
  f_loc(_, 0) = 1 - (abc_mat(_,2) / (abc_mat(_,0) + abc_mat(_,1) + abc_mat(_,2)));
  f_loc(_, 1) = abc_mat(_,0) / (abc_mat(_,0) + abc_mat(_,1) + abc_mat(_,2));
  f_loc(_, 2) = 1 - (abc_mat(_,2) / (abc_mat(_,1) + abc_mat(_,2)));
  
  // Rcout << "abc_mat\n" << abc_mat << "\n\n";
  NumericVector fstats = getFstatsFromAbc(abc_mat);
  NumericMatrix pw_fst(n_i.length(), n_i.length());
  rownames(pw_fst) = CharacterVector(n_i.names());
  colnames(pw_fst) = CharacterVector(n_i.names());
  for(size_t i = 0; i < pw_abc_mats.size(); ++i){
    double fst = getFstatsFromAbc(pw_abc_mats[i])["Fst"];
    pw_fst(pairs(i,0), pairs(i,1)) = fst;
    pw_fst(pairs(i,1), pairs(i,0)) = fst;
  }
  
  return List::create(Named("fstats") = fstats,
                      Named("fstats_loci") = f_loc,
                      Named("pw_fst") = pw_fst);
}



// 
// 
// // [[Rcpp::export]]
// List wcCpp(const IntegerMatrix als, const IntegerVector pop) {
//   // get the size of each population
//   //Rcout << "check1\n";
//   auto n_i = table(pop);
//   
//   CharacterVector n_i_nms = n_i.names();
//   //Rcout << n_i_nms;
//   //Rcout << "n_i\n" << n_i_nms << "\n" << n_i << "\n\n";
//   //Rcout << "check2\n";
//   double r = (double) n_i.length(); // the number of populations
//   //Rcout << "r\n" << r << "\n";
//   double n_bar = sum(n_i) / r;
//   //Rcout << "n_bar\n" << n_bar << "\n";
//   double n_c = ((r * n_bar) - (sum(pow(n_i,2)) / (r * n_bar))) / (r - 1);
//   //Rcout << "n_c\n" << n_c << "\n";
//   //Rcout << "check3\n";
//   IntegerVector pop_ind = seq(0,n_i.length()-1);
//   
//   //Rcout << "check4\n";
//   
//   pop_ind.names() = n_i.names();
//   CharacterVector pop_ind_nms = pop_ind.names();
//   //Rcout << "pop_ind\n" << pop_ind_nms << "\n" << pop_ind << "\n\n";
//   
//   //Rcout << "check5\n";
//   NumericMatrix f_loc(als.ncol()/2, 3);
//   colnames(f_loc) = Rcpp::CharacterVector({"Fit","Fst","Fis"}); //name the columns
//   NumericMatrix abc_mat(als.ncol()/2, 3);
//   //Rcout << "check6\n";
//   
//   for(int i = 0; i < (als.ncol() - 1); i += 2){
//     //Rcout << "   check7\n";
//     //IntegerVector al_unq = unique(als( Range(0, als.nrow()) , Range(i, i + 1) ));
//     IntegerVector al_unq1 = unique(als(_, i));
//     IntegerVector al_unq2 = unique(als(_, i+1));
//     //Rcout << "   check8\n";
//     IntegerVector al_unq = union_(al_unq1, al_unq2);
//     //Rcout << "   check9\n";
//     IntegerVector al_ind = seq(0,al_unq.length()-1);
//     
//     //Rcout << "   check10\n";
//     al_ind.names() = al_unq;
//     CharacterVector al_ind_nms = al_ind.names();
//     //Rcout << "al_ind\n" << al_ind_nms << "\n" << al_ind << "\n\n";
//     //Rcout << "   check11\n";
//     int al_n = al_unq.length();
//     NumericMatrix count_mat(al_n, r);
//     NumericMatrix het_mat(al_n, r);
//     //Rcout << "   check12\n";
//     for(int rw = 0; rw < als.nrow(); rw++) {
//       std::string pop_i = std::to_string(pop[rw]);
//       std::string al1 = std::to_string(als(rw, i));
//       std::string al2 = std::to_string(als(rw, i + 1));
//       if(al1 != al2){
//         het_mat(al_ind[al1], pop_ind[pop_i])++;
//         het_mat(al_ind[al2], pop_ind[pop_i])++;
//       }
//       for(int cn = i; cn <= i+1; cn++){
//         std::string al_i = std::to_string(als(rw, cn));
//         count_mat(al_ind[al_i], pop_ind[pop_i])++;
//       }
//     }
//     //Rcout << "count_mat\n" << count_mat << "\n\n";
//     //Rcout << "het_mat\n" << het_mat << "\n\n";
//     
//     //Rcout << "   check13\n";
//     NumericMatrix p_i(count_mat.nrow(), count_mat.ncol());
//     for(int j = 0; j < count_mat.nrow(); ++j){
//       for(int k = 0; k < count_mat.ncol(); ++k){
//         p_i(j,k) = count_mat(j,k) / (2 * n_i[k]);  
//       }
//     }
//     //Rcout << "p_i\n" << p_i << "\n\n";
//     //Rcout << "   check14\n";
//     NumericVector p_bar(p_i.nrow());
//     for(int j = 0; j < p_i.nrow(); ++j){
//       double sum = 0;
//       for(int k = 0; k < p_i.ncol(); ++k){
//         sum += n_i[k] * p_i(j,k);
//       }
//       p_bar[j] = sum / (r * n_bar);
//     }
//     //Rcout << "p_bar\n" << p_bar << "\n\n";
//     
//     //Rcout << "   check15\n";
//     NumericMatrix s2_0(p_i.nrow(), p_i.ncol());
//     for(int j = 0; j < p_i.nrow(); ++j){
//       for(int k = 0; k < p_i.ncol(); ++k){
//         s2_0(j,k) = pow(p_i(j,k) - p_bar[j], 2);  
//       }
//     }
//     //Rcout << "s2_0\n" << s2_0 << "\n\n";
//     //Rcout << "   check16\n";
//     NumericVector s2(s2_0.nrow());
//     for(int j = 0; j < s2_0.nrow(); ++j){
//       double sum = 0;
//       for(int k = 0; k < s2_0.ncol(); ++k){
//         sum += s2_0(j,k) * n_i[k];
//       }
//       s2[j] = sum / ((r - 1) * n_bar);
//     }
//     //Rcout << "s2\n" << s2 << "\n\n";
//     //Rcout << "   check17\n";
//     NumericVector h_bar(het_mat.nrow());
//     for(int j = 0; j < het_mat.nrow(); ++j){
//       double sum = 0;
//       for(int k = 0; k < het_mat.ncol(); ++k){
//         sum += het_mat(j,k);
//       }
//       h_bar[j] = sum / (r * n_bar);
//     }
//     //Rcout << "h_bar\n" << h_bar << "\n\n";
//     //Rcout << "   check18\n";
//     NumericVector a = (n_bar / n_c) * (s2 - (1 / (n_bar - 1)) * ((p_bar * (1 - p_bar)) - (((r - 1) / r) * s2) - (h_bar / 4)));
//     NumericVector b = (n_bar / (n_bar - 1)) * ((p_bar * (1 - p_bar)) - (((r - 1) / r) * s2)  - ((2 * n_bar - 1) / (4 * n_bar)) * h_bar);
//     NumericVector c = h_bar / 2;
//     //Rcout << "   check19\n";
//     double a_sum = sum(a);
//     double b_sum = sum(b);
//     double c_sum = sum(c);
//     //Rcout << "   check20\n";
//     f_loc(i/2, 0) = 1 - (c_sum / (a_sum + b_sum + c_sum));
//     f_loc(i/2, 1) = a_sum / (a_sum + b_sum + c_sum);
//     f_loc(i/2, 2) = 1 - (c_sum / (b_sum + c_sum));
//     abc_mat(i/2,_) = NumericVector::create(a_sum, b_sum, c_sum);
//     //Rcout << "   check21\n";
//   }
//   //Rcout << "check22\n";
//   //Rcout << "sum(abc_mat(_,1))\n" << sum(abc_mat(_,1)) << "\n\n";
//   //Rcout << "sum(abc_mat)\n" << sum(abc_mat) << "\n\n";
//   double fit = 1 - sum(abc_mat(_,2)) / sum(abc_mat);
//   double fst = sum(abc_mat(_,0)) / sum(abc_mat);
//   double fis = 1 - sum(abc_mat(_,2)) / (sum(abc_mat(_,1)) + sum(abc_mat(_,2)));
//   NumericVector fstats = NumericVector::create(fit, fst, fis);
//   fstats.names() = CharacterVector({"Fit", "Fst", "Fis"});
//   //Rcout << "check23\n";
//   return List::create(Named("fstats") = fstats, Named("fstats_loci") = f_loc);
// }
// 
