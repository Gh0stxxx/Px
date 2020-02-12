#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <fstream>
#include <unsupported/Eigen/MatrixFunctions>
#include <utility>

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>mat;
typedef Eigen::SparseMatrix<double>smat;





mat puissance(const mat & M,int n){
	if(n==0){return mat::Identity(M.rows(),M.cols());}
	else if(n==1){return M;}
	else{
		mat N(M.rows(),M.cols());
		N=puissance(M,n-1);
		return N*M;
		}
	}
	
	
mat puissance2(const mat & M,int n){
	if(n==0){return mat::Identity(M.rows(),M.cols());}
	else if(n%2==0){mat K(M.rows(),M.cols());
		K=puissance2(M,n/2);
		return K*K;}
	else{
		mat K(M.rows(),M.cols());
		K=puissance2(M,(n-1)/2);
		return M*K*K;
		}
	}
	
smat puissance_sparse(const smat & M,int n){
	
	if(n==0){smat I(M.rows(),M.cols());
		I.setIdentity();
		return I;}
	else if(n%2==0){smat K(M.rows(),M.cols());
		K=puissance_sparse(M,n/2);
		return K*K;}
	else{
		smat K(M.rows(),M.cols());
		K=puissance_sparse(M,(n-1)/2);
		return M*K*K;}
	}
		
		
mat crochet_lie(const mat & A, const mat & B){return A*B-B*A;}

std::pair<double,double> f(const mat & X, const mat & Y){
	mat P=X.exp()*Y.exp();
	mat S=(X+Y).exp();
	mat A=P-S;
	mat B=P-S*crochet_lie(X,Y);
	
	
	return std::make_pair(A.norm(),B.norm());
	
}
		
			
int main(){

	/*
	for(int i=0;i<30;i++){
		for(int j=0;j<30;j++){
	matdat>>dat(i,j);
	}
}
	mat A(3,3);
	A<<0.4,0.6,0,0.75,0.25,0,0,0,1;
	
	auto t1 = std::chrono::system_clock::now();
	mat B=puissance(dat,1000);
	auto t2 = std::chrono::system_clock::now();
	std::chrono::duration<double> diff=t2-t1;
	std::cout<<B<<std::endl<<"il s'est écoulé "<<diff.count()<<"s."<<"pour le calcul 1"<<std::endl;
	
	auto t3 = std::chrono::system_clock::now();
	mat C=puissance2(dat,1000);
	auto t4 = std::chrono::system_clock::now();
	std::chrono::duration<double> diff2=t4-t3;
	std::cout<<C<<std::endl<<"il s'est écoulé "<<diff2.count()<<"s."<<"pour le calcul 2"<<std::endl;
	*/ 
	

mat dat(30,30);
	std::ifstream matdat("matrice.dat");
 
 Eigen::SparseMatrix<double> dmat(30,30);
	double x;
	for(int i=0;i<30;i++){
		for(int j=0;j<30;j++){
			matdat>>x;
			if(x>1.e-10){dmat.coeffRef(i,j)=x;}
	}
	}
	
	std::cout<<puissance_sparse(dmat,1000)<<std::endl;
	}





