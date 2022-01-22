
#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// factorial
long int fact(int n){
  long int prod = 1;
  for(int i = n; i > 0; i--){
    prod *= i;
  }
  return prod;
}

// int chooseCust(int n, int k){
//   Rcout << "fact(n) = " << fact(n) << "\n";
//   Rcout << "fact(k) = " << fact(k) << "\n";
//   Rcout << "fact(n-k) = " << fact(n-k) << "\n";
//   auto val = fact(n) / (fact(k) * fact(n - k));
//   Rcout << "fact(n) / (fact(k) * fact(n - k) = " << val << "\n";
//   return fact(n) / (fact(k) * fact(n - k));
// }

// function choose(n, k)
//   if k == 0 return 1
//   return (n * choose(n - 1, k - 1)) / k

unsigned choose( unsigned n, unsigned k )
{
  if (k > n) return 0;
  if (k * 2 > n) k = n-k;
  if (k == 0) return 1;
  
  int result = n;
  for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
  }
  return result;
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


NumericVector getAbc(const NumericMatrix &p_mat, const NumericMatrix &het_mat, const NumericVector &n_i, const IntegerVector &cols){
  
  LogicalVector non_zero_cols(cols.length());
  for(int i = 0; i < cols.length(); i++){
    non_zero_cols[i] = n_i[cols[i]] != 0;
  }
  IntegerVector cols2 = cols[non_zero_cols];
  if(cols2.length() < 2){
    return NumericVector::create(NA_REAL, NA_REAL, NA_REAL);
  }
  
  // I think I need to add a check here - I could check each value of 'n_i' - 
  // if a value is 0 then we exclude the corresponding column (remove the column
  // from 'cols')
  NumericVector n_i2(cols2.length());
  for(int i = 0; i < cols2.length(); i++){
    n_i2[i] = n_i[cols2[i]];
  }
  
  double r = (double) n_i2.length(); // the number of populations
  double n_bar = sum(n_i2) / r;
  double n_c = ((r * n_bar) - (sum(pow(n_i2,2)) / (r * n_bar))) / (r - 1);

  NumericVector p_bar(p_mat.nrow());
  NumericVector s2(p_mat.nrow());
  for(int j = 0; j < p_mat.nrow(); ++j){
    NumericVector p_bar_numerator(cols2.length());
    for(int k = 0; k < cols2.length(); ++k){
      p_bar_numerator[k] = p_mat(j, cols2[k]) * n_i2[k];
    }
    p_bar[j] = sum(p_bar_numerator) / (r * n_bar);
    // p_bar[j] = sum(na_omit(p_bar_numerator)) / (r * n_bar);
    
    NumericVector s2_numerator(cols2.length());
    for(int k = 0; k < cols2.length(); ++k){
      s2_numerator[k] = pow(p_mat(j, cols2[k]) - p_bar[j], 2) * n_i2[k];  
    }
    s2[j] = sum(s2_numerator) / ((r - 1) * n_bar);
    // s2[j] = sum(na_omit(s2_numerator)) / ((r - 1) * n_bar);
  }
  
  if(r == 1){
    s2 = NumericVector(p_mat.nrow()); // if there's only one population set 's2' to be all 0s
  }
  NumericVector h_bar(het_mat.nrow());
  for(int j = 0; j < het_mat.nrow(); ++j){
    double sum = 0;
    for(int k = 0; k < cols2.length(); ++k){
      sum += het_mat(j,cols2[k]);
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
  // double a_sum = sum(na_omit(a));
  // double b_sum = sum(na_omit(b));
  // double c_sum = sum(na_omit(c));  

  NumericVector abc = NumericVector::create(a_sum, b_sum, c_sum);

  return abc;
}

NumericVector getFstatsFromAbc(const NumericMatrix &abc_mat){
  double sum0 = sum(na_omit(abc_mat(_,0)));
  double sum1 = sum(na_omit(abc_mat(_,1)));
  double sum2 = sum(na_omit(abc_mat(_,2)));
  double sum_all = sum0 + sum1 + sum2;
  
  double fit = 1 - sum2 / sum_all;
  double fst = sum0 / sum_all;
  double fis = 1 - sum2 / (sum1 + sum2);
  NumericVector fstats = NumericVector::create(fit, fst, fis);
  fstats.names() = CharacterVector({"Fit", "Fst", "Fis"});
  return fstats;
  
  // double fit = 1 - sum(abc_mat(_,2)) / sum(abc_mat);
  // double fst = sum(abc_mat(_,0)) / sum(abc_mat);
  // double fis = 1 - sum(abc_mat(_,2)) / (sum(abc_mat(_,1)) + sum(abc_mat(_,2)));
  // NumericVector fstats = NumericVector::create(fit, fst, fis);
  // fstats.names() = CharacterVector({"Fit", "Fst", "Fis"});
  // return fstats;
}

// [[Rcpp::export]]
List wcCpp(const IntegerMatrix als, const IntegerVector pop) {
  // get the size of each population
  //IntegerVector n_i = table(pop);
  NumericVector n_i = as<NumericVector>(table(pop));
  //NumericVector n_i_num = as<NumericVector>(n_i);
  // Rcout << "n_i: " << n_i << "\n";
  IntegerVector pop_ind = seq(0,n_i.length()-1);
  pop_ind.names() = n_i.names();
  
  NumericMatrix pop_stats(n_i.length(), 2);
  colnames(pop_stats) = Rcpp::CharacterVector({"Ho", "Hs"}); //name the columns
  rownames(pop_stats) = Rcpp::CharacterVector(n_i.names());
  NumericMatrix pop_ho(als.ncol()/2, n_i.length());
  NumericMatrix pop_hs(als.ncol()/2, n_i.length());
  // NumericMatrix n_rows_all(als.ncol()/2, n_i.length()); // # of non-NA rows for each population and locus
  
  NumericMatrix abc_mat(als.ncol()/2, 3);
  NumericMatrix loci_stats(als.ncol()/2, 6);
  colnames(loci_stats) = Rcpp::CharacterVector({"N_al", "Ho", "Hs", "Fit","Fst","Fis"}); //name the columns
  // NumericVector ho_vec(als.ncol()/2);
  // NumericVector he_vec(als.ncol()/2);
  
  IntegerMatrix pairs = combnIndices2(n_i.length());
  std::vector<NumericMatrix> pw_abc_mats(pairs.nrow());
  for(int i = 0; i < pairs.nrow(); ++i){
    pw_abc_mats[i] = NumericMatrix(als.ncol()/2, 3);
  }
  
  for(int i = 0; i < (als.ncol() - 1); i += 2){
    IntegerVector al_unq1 = unique(als(_, i));
    IntegerVector al_unq2 = unique(als(_, i+1));
    IntegerVector al_unq = union_(al_unq1, al_unq2);

    al_unq = al_unq[!is_na(al_unq)];
    IntegerVector al_ind = seq(0,al_unq.length()-1);
    al_ind.names() = al_unq;
    int al_n = al_unq.length();
    
    NumericMatrix count_mat(al_n, n_i.length());
    NumericMatrix het_mat(al_n, n_i.length());
    NumericVector het_n_pop(n_i.length());
    NumericVector al_counts(al_n);
    NumericVector n_rows_pop(n_i.length());
    
    int het_n = 0;
    int n_rows = 0; // number of non-NA rows
    
    
    for(int rw = 0; rw < als.nrow(); rw++) {
      if(!IntegerVector::is_na(als(rw, i)) && !IntegerVector::is_na(als(rw, i + 1))){
        
        std::string pop_i = std::to_string(pop[rw]);
        int pi = pop_ind[pop_i];
        n_rows_pop[pi]++;
        n_rows++;
        
        std::string al1 = std::to_string(als(rw, i));
        std::string al2 = std::to_string(als(rw, i + 1));
        
        if(al1 != al2){
          int a1i = al_ind[al1];
          int a2i = al_ind[al2];

          het_mat(a1i, pi)++;
          het_mat(a2i, pi)++;
          het_n++;
          het_n_pop(pi)++;
        }
        
        for(int cn = i; cn <= i+1; cn++){
          
          std::string al_i = std::to_string(als(rw, cn));
          //p_i[k] = p_mat(j,cols[k]) / (2 * n_i2[k]);
          int ai = al_ind[al_i];
          
          count_mat(ai, pi)++;
          al_counts[ai]++;
          // p_mat(al_ind[al_i], pop_ind[pop_i]) += 1/(2*n_pop_i);
        }
      }
    }

    NumericMatrix p_mat(count_mat.nrow(), count_mat.ncol());
    for(int j = 0; j < count_mat.nrow(); ++j){
      for(int k = 0; k < count_mat.ncol(); ++k){
        // p_mat(j,k) = count_mat(j,k) / (2 * n_i[k]);
        p_mat(j,k) = count_mat(j,k) / (2 * n_rows_pop[k]);
      }
    }
    
    // for(int j = 0; j < count_mat.ncol(); ++j){
    //   if(n_rows_pop[j] > 0){
    //     for(int k = 0; k < count_mat.nrow(); ++k){
    //       // p_mat(j,k) = count_mat(j,k) / (2 * n_i[k]);
    //       p_mat(k,j) = count_mat(k,j) / (2 * n_rows_pop[j]);
    //     }
    //   }
    // }
    // 
    // Rcout << p_mat << "\n";
    // Rcout << het_mat << "\n";
    // Rcout << n_rows_pop << "\n";
    abc_mat(i/2,_) = getAbc(p_mat, het_mat, n_rows_pop, seq(0,p_mat.ncol()-1));
    
    for(int j = 0; j < pairs.nrow(); ++j){
      IntegerVector cols = IntegerVector::create(pairs(j,0), pairs(j,1));
      // pw_abc_mats[j](i/2,_) = getAbc(p_mat, het_mat, n_i, cols);
      pw_abc_mats[j](i/2,_) = getAbc(p_mat, het_mat, n_rows_pop, cols);
    }
    
    // ho_vec[i/2] = (double) het_n / als.nrow();
    // he_vec[i/2] = 1 - sum(pow(al_counts / (2 * als.nrow()), 2));
    loci_stats(i/2, 0) = al_n;
    // loci_stats(i/2, 1) = (double) het_n / als.nrow();
    loci_stats(i/2, 1) = (double) het_n / n_rows;
    // loci_stats(i/2, 2) = 1 - sum(pow(al_counts / (2 * als.nrow()), 2));
    loci_stats(i/2, 2) = 1 - sum(pow(al_counts / (2 * n_rows), 2));
    
    // pop_ho(i/2, _) = het_n_pop / n_i;
    pop_ho(i/2, _) = het_n_pop / n_rows_pop;
    for(int j = 0; j < p_mat.ncol(); ++j){
      pop_hs(i/2, j) = 1 - sum(pow(p_mat(_, j), 2));
    }
  }
  
  for(int i = 0; i < pop_ho.ncol(); ++i){
    pop_stats(i,0) = mean(na_omit(pop_ho(_,i)));
    pop_stats(i,1) = mean(na_omit(pop_hs(_,i)));
    // pop_stats(i,0) = sum(pop_ho(_,i) * n_rows_all(_,i)) / sum(n_rows_all(_,i));
    // pop_stats(i,1) = sum(pop_hs(_,i) * n_rows_all(_,i)) / sum(n_rows_all(_,i));
  }
  // NumericMatrix f_loc(als.ncol()/2, 3);
  // colnames(f_loc) = Rcpp::CharacterVector({"Fit","Fst","Fis"}); //name the columns
  Rcout << abc_mat << "\n";
  loci_stats(_, 3) = 1 - (abc_mat(_,2) / (abc_mat(_,0) + abc_mat(_,1) + abc_mat(_,2)));
  loci_stats(_, 4) = abc_mat(_,0) / (abc_mat(_,0) + abc_mat(_,1) + abc_mat(_,2));
  loci_stats(_, 5) = 1 - (abc_mat(_,2) / (abc_mat(_,1) + abc_mat(_,2)));
  
  NumericVector fstats = getFstatsFromAbc(abc_mat);
  NumericVector stats(5);
  stats.names() = Rcpp::CharacterVector({"Ho", "Hs", "Fit", "Fst", "Fis"}); //name the columns
  stats[0] = mean(loci_stats(_,1));
  stats[1] = mean(loci_stats(_,2));
  stats[2] = fstats[0];
  stats[3] = fstats[1];
  stats[4] = fstats[2];
  
  NumericMatrix pw_fst(n_i.length(), n_i.length());
  rownames(pw_fst) = CharacterVector(n_i.names());
  colnames(pw_fst) = CharacterVector(n_i.names());
  for(size_t i = 0; i < pw_abc_mats.size(); ++i){
    double fst = getFstatsFromAbc(pw_abc_mats[i])["Fst"];
    pw_fst(pairs(i,0), pairs(i,1)) = fst;
    pw_fst(pairs(i,1), pairs(i,0)) = fst;
  }
  return List::create(Named("stats") = stats,
                      Named("loci_stats") = loci_stats,
                      Named("pw_fst") = pw_fst,
                      Named("pop_stats") = pop_stats);
}


