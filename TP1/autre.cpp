#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>

struct Fiche {
		std::string prenom;
		std::string ville;
		int age;
		double temps;
	};

int main(){
	std::ofstream fichier("toulousains.txt");
	int N=2500;
	
	std::ifstream F ("smalldata.txt");

	std::vector<Fiche> vdata(N);
	for(int i=0;i<2500;i++){
		F>>vdata[i].prenom >>vdata[i].ville>>vdata[i].age>>vdata[i].temps;
	}



	int nb_lyonnais=std::count_if(vdata.begin(),vdata.end(),
		[](Fiche vdata){return vdata.ville=="Lyon";});
		
		std::cout<<"Le nb de lyonnais est : "<<nb_lyonnais<<std::endl;
	
	
	int nb_jeunes_lyonnais=std::count_if(vdata.begin(),vdata.end(),
		[](Fiche vdata){return vdata.ville=="Lyon" && vdata.age<30;});
	
	std::cout<<"Le nb de jeunes lyonnais est : "<<nb_jeunes_lyonnais<<std::endl;


	bool toulousain70=std::any_of(vdata.begin(),vdata.end(),
	[](Fiche vdata){return vdata.ville=="Toulouse" && vdata.age>70;});
	if(toulousain70){std::cout<<"il existe un toulousain de plus de 70 ans"<<std::endl;};
	if(!toulousain70){std::cout<<"il n'existe pas de toulousain de plus de 70 ans"<<std::endl;};


	auto couple=std::minmax_element(vdata.begin(),vdata.end(),[](const Fiche & f1,const Fiche & f2){return(f1.age<f2.age);});
	std::cout<<"L'age minimum est : "<<couple.first->age<<std::endl;
	std::cout<<"L'age maximum est : "<<couple.second->age<<std::endl;
	std::cout<<"Le plus jeune s'apelle : "<<couple.first->prenom<<std::endl;
	std::cout<<"Le plus agé s'apelle : "<<couple.second->prenom<<std::endl;
	
	
	int sum=std::accumulate(vdata.begin(),vdata.end(),0,[](int i,int vdata.age){return i+f;});
	std::cout<<"La moyenne d'age est de : "<<double(sum)/double(N)<<std::endl;
	
}
