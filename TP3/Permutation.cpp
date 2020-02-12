#include "Permutation.hpp"

Permutation::Permutation(int n0): n(n0), images(n0){
	for(int i=0;i<images.size();i++){
		images[i]=i;
	}
}

Permutation::Permutation(const std::vector<int> & v):n(v.size()), images(v){}

Permutation Permutation::extend(int m) const {
	if(m<=n){return *this;}
	else{
		Permutation newperm(m);
		for(int i=0;i<n;i++){
		newperm[i]=images[i];
		}
		for(int i=n;i<m;i++){
			newperm[i]=i;
			}
		}
	}
