#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <unsupported/Eigen/MatrixFunctions>
#include <random>
#include <Eigen/Eigenvalues>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixDouble;
typedef Eigen::SparseMatrix<double> SparseDouble;

//Question 2.1
MatrixDouble puissance(const MatrixDouble & M, int n){
if(n==0){
    return MatrixDouble::Identity(M.rows(),M.cols());
}
else if(n==1){
    return M;
}
else return M*puissance(M, n-1);
}

//Question 2.4
MatrixDouble puissance2(const MatrixDouble & M, int n){
if(n==0){
	return MatrixDouble::Identity(M.rows(),M.cols());
}
else if(n==1){
    return M;
}
else if(n%2 == 0){
    MatrixDouble N(M.rows(),M.cols());
    N = puissance2(M, n/2);
    return N*N;
}
else { //case n%2==1
    MatrixDouble N(M.rows(),M.cols());
    N = puissance2(M, (n-1)/2);
    return M*N*N;
}
	
}

//Question 2.6
SparseDouble puissance_sparse(const SparseDouble & M, int n){
if(n==0){
	SparseDouble Id(M.rows(),M.cols());
	Id.setIdentity();	//cf. https://eigen.tuxfamily.org/dox/classEigen_1_1SparseMatrix.html pour la documentation
	return Id;
}
else if(n==1){
    return M;
}
else if(n%2 == 0){
    SparseDouble N(M.rows(),M.cols());
    N = puissance_sparse(M, n/2);
    return N*N;
}
else { //cas n%2==1
    SparseDouble N(M.rows(),M.cols());
    N = puissance_sparse(M, (n-1)/2);
    return M*N*N;
	/*	ATTENTION !!!
	 * Surtout ne pas écrire: return puissance_sparse(M, (n-1)/2)*puissance_sparse(M, (n-1)/2)*M;
	 * car sinon le calcul de puissance_sparse(M, (n-1)/2) serait fait deux fois et nous n'aurions aucune amélioration 
	 * de performances. La solution proposée le calcule une fois et ensuite on utilise deux fois le résultat. 
	 * */
}
}

//Question 2.8
MatrixDouble crochet_lie(const MatrixDouble & A, const MatrixDouble & B){
return A*B-B*A;
}

//Question 2.9
std::pair<double,double> f(const MatrixDouble & X, const MatrixDouble & Y){
    MatrixDouble M1(X.cols(),X.cols());
    MatrixDouble M2(X.cols(),X.cols());
    MatrixDouble S(X.cols(),X.cols());
    MatrixDouble P(X.cols(),X.cols());
    S=X+Y;
    P=X.exp()*Y.exp();
    M1 = P-S.exp();
    M2 = P-(S+0.5*crochet_lie(X,Y)).exp();
    /* Remarque:
     * il est utile de précalculer S et P car sinon X.exp() est calculé plusieurs fois inutilement. 
     * Le seul inconvénient est un léger encombrement mémoire supplémentaire mais en général ce n'est pas le plus important.
     * */
    std::pair<double,double> res;
    res.first=M1.norm();
    res.second=M2.norm();
    return res;
}

int main(){
    //auto start = std::chrono::system_clock::now();
    //Question 2.2
    Eigen::Matrix<double, 3,3> A;
    A << 0.4, 0.6, 0.,
        0.75, 0.25, 0.,
        0., 0., 1;
    std::cout << " A = " << std::endl << A << std::endl;
    std::cout << " A^100 = " << std::endl << puissance(A,100) << std::endl;
    std::cout << " Exp(A) = " << std::endl << A.exp() << std::endl;

    //Question 2.3
    //L'esperluette '&' sert à passer l'argument en référence.
    //Sans esperluette, tous les coefficients de la matrice et de ses puissances 1 à n sont copiés :
    //cela représente, dans le cas de A^100, 900 copies.
    //

    //Questions 2.5 et 2.6
    int N = 30;
    std::ifstream o("matrice.dat");
    MatrixDouble B(N,N);
    for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
                o >> B(i,j) ;
            }
    }

    MatrixDouble PuissB(N,N);
    auto t1 = std::chrono::system_clock::now();
    PuissB=puissance(B,1000);
    auto t2 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff1 = t2 - t1;
    std::cout << "Temps ecoule pour calculer B^1000 avec la premiere methode = " << diff1.count() << " s" << std::endl;
    /* 
     * Sur un ordinateur de bureau classique, un peu vieux mais de bonne qualité: 
     * nous obtenons autour de 0.013s. 
     * 
     * */
	MatrixDouble PuissB2(N,N);
    auto t3 = std::chrono::system_clock::now();
    PuissB2 = puissance2(B,1000);
    auto t4 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff2 = t4 - t3;
    std::cout << "Temps ecoule pour calculer B^1000 avec la deuxieme methode = " << diff2.count() << " s" << std::endl;
    /* 
     * Sur un ordinateur de bureau classique, un peu vieux mais de bonne qualité: 
     * nous obtenons autour de 0.000153824 s
     * cela est beaucoup mieux !
     * */
   
    
    SparseDouble C(N,N);// Conversion de B en matrice sparse
    for(int i=0; i<N; i++){
            for(int j=0; j<N; j++){
				if( B(i,j) >= 1.e-10 ) C.coeffRef(i,j) = B(i,j) ;
				/* ATTENTION:
				 * le test est nécessaire ! En effet, l'avantage des matrices sparse est de ne stocker que 
				 * certains coefficients (ceux qui sont non nuls). La responsabilité est laissée à
				 * l'utilisateur. Chaque appel de C.coeffRef ajoute un coefficient (quel que soit sa valeur).
				 * 
				 * Ici le test n'est pas ==0. pour se prémunir d'éventuelles erreurs d'arrondi dans la matrice B.
				 * */
				 
				/* Remarque: 
				 * La documentation de Eigen dit .coeff(i,j) est un accesseur, alors que
				 * le mutateur se nomme .coeffRef(i,j).
				 * */
            }
    }
	SparseDouble PuissC(N,N);
    auto t5 = std::chrono::system_clock::now();
    PuissC = puissance_sparse(C,1000);
    auto t6 = std::chrono::system_clock::now();
    std::chrono::duration<double> diff3 = t6-t5;
    std::cout << "Temps ecoule pour calculer B^1000 en sparse = " << diff3.count() << " s" << std::endl;
    /* 
     * Sur un ordinateur de bureau classique, un peu vieux mais de bonne qualité: 
     * nous obtenons autour de 5.4e-5 s :
     * cela est très largement mieux !
     * */


    //Question 2.10
    std::mt19937 G(time(NULL));
    std::uniform_real_distribution<double> Loi(-100,100);
    double somme1 = 0.;
    double somme2 = 0.;
    const int nb_samples=10000;
    std::pair<double,double> norms;
    MatrixDouble X(3,3);
    MatrixDouble Y(3,3);
    // Initialisation
    for(int i=0; i<3; i++){
		for (int j = 0; j < 3; j++)
		{
			X(i,j)=0.;
			Y(i,j)=0.;
		}
	}
    for(int i=0; i<nb_samples; i++){
		// Mise à jour des coefficients qui changent:
        X(0,1) = Loi(G);
        X(0,2) = Loi(G);
        X(1,2) = Loi(G);
        Y(0,1) = Loi(G);
        Y(0,2) = Loi(G);
        Y(1,2) = Loi(G);
        // Calcul des normes:
        norms = f(X,Y);
        somme1+=norms.first;
        somme2+=norms.second;
    }
    std::cout << "M_A = " << somme1/double(nb_samples)<< std::endl;
    std::cout << "M_B = " << somme2/double(nb_samples) << std::endl;



    //Question 2.11
    int N2 = 150;
    std::normal_distribution<double> Loi_diag(0,1);
    std::normal_distribution<double> Loi_hors_diag(0,2);
    
    
    int nb_boxes = 20;
    double a=-3.;
    double b=3.;
    int simul = 50;//On peut remplacer simul par un nombre plus grand, comme 40 ou 50, mais pas trop non plus sinon le temps de calcul devient trés long
    std::vector<double> hist(nb_boxes,0);
    MatrixDouble GOE(N2,N2);
    
	//Génération des matrices et calcul des valeurs propres
    for(int i=0; i<simul; i++){

        for(int j=0; j<N2; j++){
			GOE(j,j) = Loi_diag(G);
            for(int k=j+1; k<N2; k++){
                GOE(j,k) = Loi_hors_diag(G);
                GOE(k,j) = GOE(j,k);
            }
        }

        Eigen::EigenSolver<MatrixDouble> Solver(GOE);
        
        double lambda_norm;
        for(int i=0; i<N2; i++){
			lambda_norm= Solver.eigenvalues()[i].real()/(2.*sqrt(N2));
            int indice = floor( (lambda_norm-a)/(b-a)*nb_boxes); // calcul de l'indice par la partie entière
            if( (indice >= 0) && (indice < nb_boxes) ){
                hist[indice]+= 1./double(simul*N2);
            }
        }
    }

    std::ofstream o2("eigenvalues.dat");
    for(int i = 0; i< nb_boxes; i++){
        o2 << a+(2*i+1)*(b-a)/(2.*nb_boxes) << "\t" << hist[i] << std::endl;
    }

return 0;
}
