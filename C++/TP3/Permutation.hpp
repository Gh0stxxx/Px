#include <vector>

class Permutation{
	private:
	int n;
	std::vector<int> images;
	
	public:
	Permutation(int n0);
	Permutation(const std::vector<int> & v);
	
	
	int size() const {return n;};
	int operator[](int i) const {return images[i];};
	int & operator[](int i){return images[i] ;};
	Permutation extend(int k) const;
	};
