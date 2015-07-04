#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

int main()
{ 
	stringstream buffer;
	string dataset = "s1";
	buffer << "datasets/" << dataset << ".txt";
	
	ifstream fin(buffer.str().c_str());

    float var1;
    float var2;

    while (fin >> var1 >> var2)
    {
      
        cout << var1 << "\t" << var2 << "\n";
    }
    
}
