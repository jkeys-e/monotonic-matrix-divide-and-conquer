// MatrixSearch-Skeleton.cpp
// Joe Song
// Created: Sept 11, 2016 
/*******************/
// MatrixSearch-jkeys.cpp
// Jeremy Keys (based on skeleton provided by Dr. Joe Song, see above)
// Last modified: Oct 13, 2016 (Sept 18, 2017 add comment header)

#include <Rcpp.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include <climits>
#include <cassert>
#include <cmath>
#include <iomanip> //to get setw, so cout can have uniform width (better print fcn)
using namespace std;


template <typename T>
class matrix {
private:
  vector<T> m_vec;
  size_t m_nrow;
  size_t m_ncol;
  
public:
  matrix(const vector<T> &v, size_t nrow, size_t ncol):
    m_vec(v), m_nrow(nrow), m_ncol(ncol) {
  } 
  
  size_t nrow() const { return m_nrow; }
  size_t ncol() const { return m_ncol; }
  
  T & operator () (size_t i, size_t j) {
    return m_vec[i + j * m_nrow];
  }
  
  const T & operator () (const size_t i, const size_t j) const {
    return m_vec[i + j * m_nrow];
  }
  
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const matrix<T> & m)
{
  for(size_t i=0; i<m.nrow(); ++i) {
    for(size_t j=0; j<m.ncol(); ++j) {
      out << setw(6) << m(i, j);
    }
    out << endl << endl;
  }
  return out;
} 

template<typename T> 
vector<T> find_row_maxima_itr(const matrix<T> & m) 
{
  // Your code here:
  // Iteratively find the maximum of each row row by row:
  vector<T> result; //0th index corresponds to 1st row, 1st index-2nd row, ... n-1th index to nth row
 
  int max_col = 0;
  
  for (unsigned int i = 0; i !=  m.nrow(); i++) {
	  for (unsigned int j = 0; j != m.ncol(); j++) {
		  if(m(i,j) > m(i, max_col))
			  max_col = j;
		  
	  }
	  
	  result.push_back(m(i, max_col));
	  max_col = 0;
  }

  
  return result;
}


template<typename T> 
vector<T> find_row_maxima(const matrix<T> & m) 
{
  // Your code here:
  // divide-and-conquer on monotonic matrix
  
  return(find_row_maxima_helper(m, 0, m.nrow() - 1, 0, m.ncol() - 1, 0));

 
  
}

void print_matrix(const vector<double> & v, size_t nrow, size_t ncol)
{
  matrix<double> mat(v, nrow, ncol);

  cout << "Input matrix is:" << endl << endl;
  cout << mat << endl;
}

vector<double> row_maxima_itr(const vector<double> & v, size_t nrow, size_t ncol)
{
  matrix<double> mat(v, nrow, ncol);
  return find_row_maxima_itr(mat);
}

template<typename T> 
vector<T> find_row_maxima_helper(const matrix<T> & m, int row_start, int row_end, int col_start, int col_end, int depth) {
	vector<T> result; //empty vector

	if(row_start > row_end) { //n == 0  (subproblem which doesn't exist)
		return result;
	}
	else if(row_start == row_end) {	
		int col_max = col_start;
		for (int i = col_start; i <= col_end; i++) {
			if ( m(row_start, i) > m(row_start, col_max)) {
				col_max = i;
			}
		}
		
		result.push_back(m(row_start, col_max));
		return result;
	}
	
	//find the middle row
	int row_mid = floor((row_end + row_start)/2); 	
	int col_mid = col_start;	
	
	for(int i = col_start; i <= col_end; i++) {
		if( m(row_mid, i) > m(row_mid, col_mid) ) {
			col_mid = i;
		}
	}
	
	int sub1_row_start = row_start;
	int sub1_row_end   = row_mid - 1;
	int sub1_col_start = col_start;
	int sub1_col_end   = col_mid;
	
	int sub2_row_start = row_mid + 1;
	int sub2_row_end   = row_end;
	int sub2_col_start = col_mid;
	int sub2_col_end   = col_end;
	
	T mid_max = m(row_mid, col_mid);
	vector<T> sub1 = find_row_maxima_helper(m, sub1_row_start, sub1_row_end, sub1_col_start, sub1_col_end, depth+4);
	vector<T> sub2 = find_row_maxima_helper(m, sub2_row_start, sub2_row_end, sub2_col_start, sub2_col_end, depth+4);

	for(unsigned int i = 0; i != sub1.size(); i++)
		result.push_back(sub1[i]);	
	
	result.push_back(mid_max);	
	
	for(unsigned int i = 0; i != sub2.size(); i++)
		result.push_back(sub2[i]);
	
	/*
	cerr << setw(depth) << "end of row_maxima_helper" << endl;
	for(unsigned int i = 0; i != result.size(); i++)
		std::cerr << setw(5) << result[i];
	
	std::cerr << endl;
	*/
	
	return result;
}

vector<double> row_maxima(const vector<double> & v, size_t nrow, size_t ncol)
{
  matrix<double> mat(v, nrow, ncol);
    
  return find_row_maxima(mat);
}

bool test_row_maxima(vector<double> (*rmfun) (const vector<double> & v, size_t nrow, size_t ncol)) 
{
  bool passed = true;
  /* Monotonic matrix example 1:
   0,  4,  -1,  2.5,   -4, 
  -3,  8, -10,    2,    7,
  -4, -3,  -1, -100, -5.5,
   0,  2, 0.3,   -3,  2.5,
   1,  0,   1,    2,    3,
  -8,  9,   2,    5,   10};
  */
  // x is column major vectorization of the matrix
  double x[] = {   0,  -3,   -4,   0, 1, -8,  //-8 == m(0, 5)
                   4,   8,   -3,   2, 0,  9,  //9  == m(1, 5)
                  -1, -10,   -1, 0.3, 1,  2,  //2  == m(2, 5)
                 2.5,   2, -100,  -3, 2,  5,  //5  == m(3, 5)	
                  -4,   7, -5.5, 2.5, 3, 10}; //10 == m(4, 5)	
  
  vector<double> v(x, x+30);
	
  double rmax_truth[] = {4, 8, -1, 2.5, 3, 10};
  
  print_matrix(v, 6, 5);
  
  if(rmfun(v, 6, 5) != vector<double>(rmax_truth, rmax_truth+6)) {
    cout << "ERROR: failed test 1!" << endl;
    passed = false;
  }
  
  return passed;
}

// [[Rcpp::export]]
bool testall() 
{
  bool passed = true;
  if(!test_row_maxima(row_maxima_itr)) {
    cout << "ERROR: row_maxima_itr() failed some test!" << endl;
    passed = false;
  } 

  if(!test_row_maxima(row_maxima)) {
    cout << "ERROR: row_maxima() failed some test!" << endl;
    passed = false;
  } 
  
  if(passed) {
    cout << "All tests passed. Congratulations!" << endl;
  }
  return passed;
}


/*
 
// [Rcpp::export]
void WrapperItr(NumericVector nmv, int nrow, int ncol) {
  //run and throw away result, this is just for measuring runtime
  vector<double> v;
  
  for (int i = 0; i != nmv.size(); i++)
    v.push_back(nmv[i]);
  
  row_maxima_itr(v, nrow, ncol);
}

// [Rcpp::export]
void WrapperDiv(NumericVector nmv, int nrow, int ncol) {
  //run and throw away result, this is just for measuring runtime
  vector<double> v;
  
  for (int i = 0; i != nmv.size(); i++)
    v.push_back((double) nmv[i]);
  
  row_maxima(v, nrow, ncol); 
}
*/

 

/*
// [Rcpp::export]
void WrapperItr(const matrix<double> & m) {
  //run and throw away result, this is just for measuring runtime
  find_row_maxima_itr(m);
}

// [Rcpp::export]
void WrapperDiv(const matrix<double> & m) {
  //run and throw away result, this is just for measuring runtime
  find_row_maxima(m); 
}
*/


//[[Rcpp::export]]
void WrapperItr(NumericMatrix& m, int nrow, int ncol) {
  //run and throw away result, this is just for measuring runtime
  vector<double> v;
  
  for (int i = 0; i != nrow*ncol; i++)
    v.push_back(m[i]);
  
  
    row_maxima_itr(v, nrow, ncol);
    return;
}

//[[Rcpp::export]]
void WrapperDiv(NumericMatrix& m, int nrow, int ncol) {
  //run and throw away result, this is just for measuring runtime
  vector<double> v;
  
  for (int i = 0; i != nrow*ncol; i++)
    v.push_back(m[i]);
  
  
  row_maxima(v, nrow, ncol);
  return;
}


int main()
{
	
  testall();
  return 0;
}
  
/*** R

if(!testall()) stop()

random.monotone.matrix <- function(nrow, ncol) 
{
  m <- matrix(rnorm(nrow*ncol), nrow=nrow, ncol=ncol)
  row.maxima.indices <- apply(m, 1, which.max)
  o <- order(row.maxima.indices)
  m <- m[o, ]
}

# Your R code for run time evaluation and visualization

num_samples	<- 25
num_steps 	<- 10
n_step 		<- 500
ns 			<- seq(1999, 6499, n_step)

runtime_div 	<- vector(length=num_samples)
runtime_itr 	<- vector(length=num_samples)
runtime_div_avg <- vector(length=num_steps)
runtime_itr_avg <- vector(length=num_steps)

vn1 <- vector(length=num_steps)
vn2 <- vector(length=num_steps)

c1 <- .000000001
c2 <- .000000001

for (i in 1:num_steps) {
  m <- random.monotone.matrix(ns[i],ns[i])
  n1 <- c1 * (ns[i] * ns[i])       #O(n*m)
  n2 <- c2 * (log2(ns[i])*ns[i])   #O(nlogm)
  vn1[i] <- n1
  vn2[i] <- n2
  for (j in 1:num_samples) {
    runtime_itr[j] 	 <- system.time(WrapperItr(m, as.integer(ns[i]), as.integer(ns[i])))[["user.self"]]
    runtime_div[j] 	 <- system.time(WrapperDiv(m, as.integer(ns[i]), as.integer(ns[i])))[["user.self"]]
  }
  
  runtime_itr_avg[i] <- sum(runtime_itr) / num_samples
  runtime_div_avg[i] <- sum(runtime_div) / num_samples
}

plot(ns, runtime_itr_avg, xlab="NxN", ylab = "seconds", col = "red", main = "Iter (Red) vs. Div & Conquer (Blue)")
points(ns, runtime_div_avg, col = "blue")
points(ns, vn1, col = "yellow")
points(ns, vn2, col ="green")
#curve(  , 0, ns[num_steps], add = TRUE )
#curve( (c2 * (log2(ns[i])*ns[i])) , 0, ns[num_steps], add = TRUE )

fcn1 <- function() {
  n <- (c1 * (ns[i] * ns[i]))
}

grid(col="blue")
*/
