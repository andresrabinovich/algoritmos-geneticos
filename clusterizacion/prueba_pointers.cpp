#include <iostream>
#include <vector>
using namespace std;

class punto{
	public:
	float x;
	float y;
};

void cargar(vector<punto *> &numeros){
	numeros.push_back(new punto());
}
int main(){

	vector<punto *> numeros;
	for(int i = 0; i < 100000000.0L; i++){
		cargar(numeros);
	}
	for(vector<punto*>::const_iterator it = numeros.begin(); it != numeros.end(); it++)
	{
	    delete *it;
	} 
	numeros.clear();
	return 0;
}
