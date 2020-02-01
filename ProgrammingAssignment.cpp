/*
Author: Linus Chen
AP Physics Assignment 0
C Programming Assignment (Bonus included)
Language: C++
Function: Accepts keyed entry of data and solves for a general polynomial without using built-in functions
Method:
Sets up Vandermonde matrix
Computes transpose
Multiplies the two matricies together
Compute inverse
Solves for and prints out the general polynomial
*/
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;
//2-dimensional vector used to store matrix
typedef vector<vector<double> > vvd;

//Prints given matrix
void printmat(vvd mat){
	cout.precision(3);
	for(int i=0; i<mat.size(); i++){
		for(int j=0; j<mat[0].size(); j++)
			cout << setw(10) << mat[i][j];
		cout << endl;
	} cout << endl;
}

//Returns the Vandermonde matrix for any degree
vvd vandermonde(int degree, vector<double> x){
	vvd v(x.size());
	for(int i=0; i<x.size(); i++)
		for(int k=degree; k>=0; k--) //Supports kth degree
			v[i].push_back(pow(x[i],k));
	return v;
}

//Returns the transpose of a given matrix
vvd transpose(vvd matrix){
	vvd t(matrix[0].size());
	//Change rows into columns and columns into rows
	for(int i=0; i<t.size(); i++)
		for(int j=0; j<matrix.size(); j++)
			t[i].push_back(matrix[j][i]);
	return t;
}

//Returns the product of matrix a and matrix b
vvd multiply(vvd a, vvd b){
	vvd product(a.size());
	for(int i=0; i<a.size(); i++)
		for(int j=0; j<b[0].size(); j++)
			product[i].push_back(0);
	if(a[0].size()!=b.size())
		return product; //Error: invalid operation
	double cur;
	for(int i=0; i<a.size(); i++)
		for(int j=0; j<b[0].size(); j++)
			for(int k=0; k<a[0].size(); k++)
				product[i][j] += a[i][k] * b[k][j];
	return product;
}

//Finds determinant of matrix by using Gaussian elimination
double det(vvd matrix){
	if(matrix.size()!=matrix[0].size())
		return 0; //Error, not a square matrix
	if(matrix.size()==0) //Check if 0x0 matrix
		return 1; //The determinant is 1
	if(matrix.size()==1) //Check if 1x1 matrix
		return matrix[0][0]; //Determinant is single value

	double cur, cur2, determinant=1, temp, divlater=1;

	//Gaussian elimination
	for(int i=0; i<matrix.size(); i++){
    		cur = matrix[i][i];
		if(cur==0) //Swap rows
      		for(int m=i+1; m<matrix.size(); m++){
        		if(matrix[m][i]!=0){
          			for(int l=0; l<matrix.size(); l++){
        				temp = matrix[i][l];
        				matrix[i][l] = matrix[m][l];
        				matrix[m][l] = temp;
          			}
          			cur = matrix[i][i];
          			divlater *= -1;
          			break;
       			}		
      		}		
    		for(int j=i+1; j<matrix.size(); j++){ //Row reduction
      			divlater *= matrix[i][i];
     			cur2 = matrix[j][i];
      			for(int k=0; k<matrix.size(); k++){
          			matrix[j][k] *= cur;
          			matrix[j][k] -= cur2 * matrix[i][k];
      			}		
    		}
  	}//Row-echelon form reached
  	for(int i=0; i<matrix.size(); i++)
  		determinant *= matrix[i][i];
  	determinant /= divlater;
  	return determinant;
}

//Returns the inverse of a given matrix
//Gauss-Jordan elimination with augmented Matrix
vvd inv(vvd matrix){
	//set up identity matrix for augmented matrix
	vvd inverse(matrix.size());
	for(int i=0; i<matrix.size(); i++)
		for(int j=0; j<matrix.size(); j++)
			inverse[i].push_back(0);
	for(int i=0; i<matrix.size(); i++) inverse[i][i] = 1;

	if(matrix.size()!=matrix[0].size()||det(matrix)==0)
		return matrix; //Error

	//Gauss-Jordan elimination
	//Perform same operations on
	//original and identity matrices
	double cur, cur2;
	for(int i=0; i<matrix.size(); i++){
    		cur = matrix[i][i];
		if(cur==0) //Swap rows
      		for(int m=i+1; m<matrix.size(); m++){
        		if(matrix[m][i]!=0){
          			for(int l=0; l<matrix.size(); l++){
            				swap(matrix[i][l],matrix[m][l]);
            				swap(inverse[i][l],inverse[m][l]);
          			}
          			cur = matrix[i][i];
          			break;
        		}
      		}
    		for(int j=i+1; j<matrix.size(); j++){ //Row reduction
      			cur2 = matrix[j][i];
      			for(int k=0; k<matrix.size(); k++){
        			matrix[j][k] *= cur;
        			inverse[j][k] *= cur;
        			matrix[j][k] -= cur2 * matrix[i][k];
        			inverse[j][k] -= cur2 * inverse[i][k];
      			}
    		}
	}//At this point, matrix is now in row-echelon form

 	//Reduce matrix into identity metrix
  	for(int i=matrix.size()-1; i>=0; i--){
    		cur = matrix[i][i];
    		for(int j=0; j<matrix.size(); j++) {
    			matrix[i][j] /= cur;
    			inverse[i][j] /= cur;
    		}
    		for(int j=i-1; j>=0; j--){
      			cur2 = matrix[j][i];
      			for(int k=0; k<matrix.size(); k++) {
      				matrix[j][k] -= cur2 * matrix[i][k];
      				inverse[j][k] -= cur2 * inverse[i][k];
      			}
    		}
  	}
  	//At this point, original matrix is identity matrix
  	//and the identity matrix is now the inverse
  	return inverse;
}

//Finds the coefficient vector and prints results
void solve(vvd vt, vvd vtv, vector<double> y){
	vvd inverse = inv(vtv);
	vvd y_v(y.size());
	for(int i=0; i<y_v.size(); i++) y_v[i].push_back(y[i]);

	//Uses the equation: c = (VtV)^(-1)*Vt*y
	vvd coeffs = multiply(inverse,vt);
	coeffs = multiply(coeffs, y_v);

	printf("Coefficients of polynomial in order of a_k, a_(k-1), ... , a_0:\n");
	printmat(coeffs);

	//Display the polynomial
	cout << "The polynomial is" << endl;
	for(int i=0; i<coeffs.size(); i++){
		cout << "(" << coeffs[i][0] << ")";
		if(i<coeffs.size()-1)
			cout << "x";
		if(i<coeffs.size()-2)
			cout << "^" << (coeffs.size()-1-i);
		if(i<coeffs.size()-1)
			cout << " + ";
	}
	cout << endl;
}

int main(){
	//Allows for fast input
	ios_base::sync_with_stdio(false);
	cin.tie(NULL);
	int N; //holds number of data points
	printf("Enter number of data points: ");
	cin >> N;
	//Input data
	vector<double> x(N), y(N);
	printf("Enter data:\n");
	printf("X: ");
	for(int i=0; i<N; i++) cin >> x[i];
	printf("Y: ");
	for(int i=0; i<N; i++) cin >> y[i];
	cout << endl << endl;

	//Allows solving for general polynomial
	printf("Enter degree of polynomial: ");
	int deg; cin >> deg;
	vvd v = vandermonde(deg,x);
	vvd vt = transpose(v);
	vvd vtv = multiply(vt,v);

	//Print results
	printf("\nInput Data:\n");
	cout << "X: ";
	for(int i=0; i<N; i++) cout << x[i] << " "; cout << endl;
	cout << "Y: ";
	for(int i=0; i<N; i++) cout << y[i] << " "; cout << endl;
	cout << endl;
	printf("V:\n");
	printmat(v);
	printf("Vt:\n");
	printmat(vt);
	printf("VtV:\n");
	printmat(vtv);
	solve(vt, vtv, y);
	cout << "\nEnd of program reached\n";
}
