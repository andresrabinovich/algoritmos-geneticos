#include <iostream>
#include <random>
#include <functional> //Para std::bind
#include <vector>
#include <cmath> //Para infinity
#include <string>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <thread>

using namespace std;

/**************
CONFIGURACIONES
**************/

/*-------------------------
Configuraciones del dataset
-------------------------*/
int centros_por_dimension = 4; //los centros son los lugares alrededor de los cuales se van a armar los clusters
int puntos_por_cluster = 50;
float parametro_de_red = 1;
float ancho_del_cluster = 0.1; //lo que mide el cluster en x
float alto_del_cluster = 0.1; //lo que mide el cluster en y

/*-------------------
Configuraciones de AG
-------------------*/
int poblacion = 40;
float pm = 0.1; //probabilidad de mutacion
float pc = 0.3; //probabilidad de single-point crossover
int pp = 3; //Cromosomas a seleccionar aleatoreamente en cada busqueda de padres
int k_max = 20; //Maxima cantidad de clusters a buscar
float alfa = 0.1; //Amplitud de mutacion de valores de activacion
float epsilon = 0.3; //Amplitud de mutacion de centroides
int soluciones_de_elite = 4; //Las mejores soluciones pasan sin alteraciones a la proxima generacion
int generaciones = 5000;

/*--------------------------
Configuraciones del programa
--------------------------*/
bool multithreding = true;

/*--------------------------------
Calcula algunas cantidades previas
--------------------------------*/
int cantidad_de_centros = centros_por_dimension*centros_por_dimension;
int cantidad_de_puntos  = cantidad_de_centros*puntos_por_cluster;

/*-------------------------------------------
Inicializa el generador de numeros aleatorios
-------------------------------------------*/
//default_random_engine generador_aleatorio(123456);

/*******************
COMIENZAN LAS CLASES
*******************/

/*---------
Clase punto
---------*/
class c_punto{
	public:
	float x;
	float y;
	
	c_punto();
	c_punto(float x, float y);
};
//Constructores de punto
c_punto::c_punto(){
	x = 0;
	y = 0;
}	

c_punto::c_punto(float x, float y){
	c_punto::x = x;
	c_punto::y = y;
}	

/*--------------
Clase rectangulo
--------------*/
class c_rectangulo{
	public:
	c_punto inferior_izquierdo;
	c_punto superior_derecho;
	
	c_rectangulo(float x1, float y1, float x2, float y2);
};
//Constructor de rectangulo
c_rectangulo::c_rectangulo(float x1, float y1, float x2, float y2){
	inferior_izquierdo.x 	= x1;
	inferior_izquierdo.y 	= y1;
	superior_derecho.x 	= x2;
	superior_derecho.y 	= y2;
}
/*-----------------------
Clase generador aleatorio
-----------------------*/
class c_generador_aleatorio{
	public:
	default_random_engine generador_aleatorio;
	float pm;
	float alfa;
	float epsilon_x;
	float epsilon_y;
	c_rectangulo* valores_limite;
	c_generador_aleatorio();
	c_generador_aleatorio(int seed);
	void cargar_configuraciones_de_cromosomas(float pm, float alfa, float epsilon, int k_max, c_rectangulo* valores_limite);

	float runif(); //Llama a runif(0, 1);
	float runif(float limite_inferior, float limite_superior);	
	
	//Todas estas funciones generan reales aleatorios entre los limites establecidos al cargar la configuracion de los cromosomas
	float runif_x();
	float runif_y();
	float runif_alfa();
	float runif_epsilon_x();
	float runif_epsilon_y();
};

c_generador_aleatorio::c_generador_aleatorio(){
	generador_aleatorio.seed(std::random_device{}());
}

c_generador_aleatorio::c_generador_aleatorio(int seed){
	generador_aleatorio.seed(seed);
}

void c_generador_aleatorio::cargar_configuraciones_de_cromosomas(float pm, float alfa, float epsilon, int k_max, c_rectangulo* valores_limite){
	c_generador_aleatorio::generador_aleatorio = generador_aleatorio;
	c_generador_aleatorio::valores_limite = valores_limite;
	c_generador_aleatorio::pm = pm;
	c_generador_aleatorio::alfa = alfa;
	c_generador_aleatorio::epsilon_x = epsilon*(valores_limite->superior_derecho.x-valores_limite->inferior_izquierdo.x);
	c_generador_aleatorio::epsilon_y = epsilon*(valores_limite->superior_derecho.y-valores_limite->inferior_izquierdo.y);	
	return;
}

float c_generador_aleatorio::runif(){
	return c_generador_aleatorio::runif(0, 1);
}
float c_generador_aleatorio::runif(float limite_inferior, float limite_superior){
	uniform_real_distribution<double> urd(limite_inferior, limite_superior);
	//auto runif 	= bind(urd, ref(generador_aleatorio));
	return urd(generador_aleatorio);
}	

//Todas estas funciones generan reales aleatorios entre los limites establecidos al cargar la configuracion de los cromosomas
float c_generador_aleatorio::runif_x(){
	uniform_real_distribution<double> urd_x(valores_limite->inferior_izquierdo.x, valores_limite->superior_derecho.x);
	//auto runif_x 	= bind(urd_x, ref(generador_aleatorio));
	return urd_x(generador_aleatorio);
}
float c_generador_aleatorio::runif_y(){
	uniform_real_distribution<double> urd_y(valores_limite->inferior_izquierdo.y, valores_limite->superior_derecho.y);
	//auto runif_y 	= bind(urd_y, ref(generador_aleatorio));		
	return urd_y(generador_aleatorio);
}
float c_generador_aleatorio::runif_alfa(){
	uniform_real_distribution<double> urd_alfa(-alfa, alfa);
	//auto runif_alfa 	= bind(urd_alfa, ref(generador_aleatorio));			
	return urd_alfa(generador_aleatorio);
}
float c_generador_aleatorio::runif_epsilon_x(){
	uniform_real_distribution<double> urd_epsilon_x(-epsilon_x, epsilon_x);
	//auto runif_epsilon_x 	= bind(urd_epsilon_x, ref(generador_aleatorio));			
	return urd_epsilon_x(generador_aleatorio);
}
float c_generador_aleatorio::runif_epsilon_y(){
	uniform_real_distribution<double> urd_epsilon_y(-epsilon_y, epsilon_y);		
	//auto runif_epsilon_y 	= bind(urd_epsilon_y, ref(generador_aleatorio));					
	return urd_epsilon_y(generador_aleatorio);
}

/*-------------
Clase cromosoma
-------------*/
class c_cromosoma {
	public:
	int k_max;
	int clusters;
	vector<float> valores_de_activacion;
	vector<c_punto> centroides;
	float promedio_al_centro_del_dataset;
	vector<c_punto>* puntos;
	double fitness;
	c_generador_aleatorio* generador_aleatorio;
	c_rectangulo* valores_limite;
	vector<int> solucion_formato_cluster; //Representacion de la solucion en formato cluster(a que cluster va cada punto)
		
	c_cromosoma();
	c_cromosoma(int k_max, c_rectangulo* valores_limite, float promedio_al_centro_del_dataset, vector<c_punto>* puntos, c_generador_aleatorio* generador_aleatorio);
	c_cromosoma(int k_max, c_rectangulo* valores_limite, float promedio_al_centro_del_dataset, vector<c_punto>* puntos, c_generador_aleatorio* generador_aleatorio, const vector<float> &valores_de_activacion, const vector<c_punto> &centroides);
	double calcular_fitness();
	string to_string();
	void mutar();
	string to_cluster();
};

//El constructor por defecto viene vacio. Sirve para definir cromosomas temporales
c_cromosoma::c_cromosoma(){}

//Este constructor genera las soluciones de forma random. Es el constructor para generar la tirada inicial de cromosomas
c_cromosoma::c_cromosoma(int k_max, c_rectangulo* valores_limite, float promedio_al_centro_del_dataset, vector<c_punto>* puntos, c_generador_aleatorio* generador_aleatorio){
	
	c_cromosoma::k_max = k_max;
	c_cromosoma::promedio_al_centro_del_dataset = promedio_al_centro_del_dataset;
	c_cromosoma::puntos = puntos;
	c_cromosoma::valores_limite = valores_limite;
	c_cromosoma::generador_aleatorio = generador_aleatorio;

	//Genera el cromosoma al azar de la corrida. Primero los valores de activacion.
	//Los dos primeros siempre estan activos para que haya al menos dos activos
	valores_de_activacion.push_back(0.5);
	valores_de_activacion.push_back(0.5);

	for(int valor = 2; valor < k_max; valor++){
		valores_de_activacion.push_back(generador_aleatorio->runif());
	}
	
	//Ahora genera los centroides. Los genera dentro de los rangos maximos y minimos de los valores limite
	for(int cluster=0; cluster < k_max; cluster++){
		centroides.push_back(c_punto(generador_aleatorio->runif_x(), generador_aleatorio->runif_y()));
	}
	
	//Calcula el fitness inicial
	calcular_fitness();
}

//Override del constructor por defecto para cargar valores_de_activacion y centroides
c_cromosoma::c_cromosoma(int k_max, c_rectangulo* valores_limite, float promedio_al_centro_del_dataset, vector<c_punto>* puntos, c_generador_aleatorio* generador_aleatorio, const vector<float> &valores_de_activacion, const vector<c_punto> &centroides){
	
	c_cromosoma::k_max = k_max;
	c_cromosoma::promedio_al_centro_del_dataset = promedio_al_centro_del_dataset;
	c_cromosoma::puntos = puntos;
	c_cromosoma::valores_limite = valores_limite;

	//Genera el cromosoma al azar de la corrida. Primero los valores de activacion.
	c_cromosoma::valores_de_activacion = valores_de_activacion;
	//Ahora genera los centroides. Los genera dentro de los rangos maximos y minimos de los valores limite
	c_cromosoma::centroides = centroides;	
	
	//Calcula el fitness inicial
	calcular_fitness();
}

//Devuelve la solucion codificada en una tira de puntos pertenecientes a clusters
string c_cromosoma::to_cluster(){
	
	//Incializa la tira que representa la solucion en el formato de clusters
	solucion_formato_cluster.clear();

	//Asigna cada punto a un cluster por proximidad
	for(int punto = 0; punto < puntos->size(); punto++){
		//Inicializa la distancia al cluster al que pertenece (infinity = no pertenece a ninguno)		
		double distancia_a_k = INFINITY;
		int k;
		for(int i = 0; i < k_max; i++){
			//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
			if(valores_de_activacion[i] >= 0.5){
				//cout << (*puntos)[punto].x << "|" << (*puntos)[punto].y << "|" << centroides[i].x << "|" << centroides[i].y << "|" << ((*puntos)[punto].x - centroides[i].x) << "|" << ((*puntos)[punto].y - centroides[i].y) << "\n";
				//Si la distancia a este centroide es la mas chica, se lo asigno al centroide
				double distancia_a_i = (pow(((*puntos)[punto].x - centroides[i].x), 2) + pow(((*puntos)[punto].y - centroides[i].y), 2));
				if(distancia_a_i < distancia_a_k) {
					distancia_a_k = distancia_a_i;
					k = i;
				}
			}		
		}
		solucion_formato_cluster.push_back(k);
	}
	
	/*
	//Ordena la lista para que no queden agujeros entre los clusters (el 111133335555 tiene que ser 111122223333)
	int k_actual = 1;
	for(int i = 0; i < k_max; i++){
		if(valores_de_activacion[i] >= 0.5){
			bool cambio_k = false;
			for(int solucion = 0; solucion < solucion_formato_cluster.size(); solucion++){
				if(solucion_formato_cluster[solucion] == i){
					solucion_formato_cluster[solucion] = k_actual;
					cambio_k = true;
				}
			}
			if(cambio_k) k_actual++;
		}
	}
	*/
	stringstream buffer;
	for(int solucion = 0; solucion < (solucion_formato_cluster.size()-1); solucion++){
		buffer << solucion_formato_cluster[solucion] << "\t";
	}
	buffer << solucion_formato_cluster[solucion_formato_cluster.size()-1];		
	return buffer.str();
}
/*-------------------------------------------------------------------------------------------------
Funcion to_string simple que devuelve como tira para poder imprimir en pantalla o exportar a un txt
-------------------------------------------------------------------------------------------------*/
string c_cromosoma::to_string(){
	stringstream buffer;
	for(int k = 0; k < k_max; k++){
		buffer << valores_de_activacion[k] << ",";
	}
	for(int k = 0; k < k_max; k++){
		buffer << centroides[k].x << ", " << centroides[k].y << ",";
	}

	return buffer.str();
}

/*-----------------
Funcion de mutacion
-----------------*/
void c_cromosoma::mutar(){

	//Muta los valores de activacion y la ubicacion de los centroides
	for(int k = 0; k < k_max; k++){
		//Muta valor_de_activacion kesimo
		if(generador_aleatorio->runif() < pm){ //Muta con probabilidad pm
	
			//muta
			float nuevo_va = valores_de_activacion[k] + generador_aleatorio->runif_alfa();
			if(nuevo_va > 1){			
				valores_de_activacion[k] = 1;
			}else if(nuevo_va < 0){
				valores_de_activacion[k] = 0;
			}else{
				valores_de_activacion[k] = nuevo_va;
			}
		}
		//Muta centroide kesimo
		if(generador_aleatorio->runif() < pm){ //Muta con probabilidad pm
	
			float nuevo_x = centroides[k].x + generador_aleatorio->runif_epsilon_x();
			float nuevo_y = centroides[k].y + generador_aleatorio->runif_epsilon_y();
			if(nuevo_x > valores_limite->superior_derecho.x){			
				nuevo_x = valores_limite->superior_derecho.x;
			}else if(nuevo_x < valores_limite->inferior_izquierdo.x){
				nuevo_x = valores_limite->inferior_izquierdo.x;
			}

			if(nuevo_y > valores_limite->superior_derecho.y){			
				nuevo_y = valores_limite->superior_derecho.y;
			}else if(nuevo_y < valores_limite->inferior_izquierdo.y){
				nuevo_y = valores_limite->inferior_izquierdo.y;
			}
			centroides[k].x = nuevo_x;
			centroides[k].y = nuevo_y;			

		}		
		
	}
	
	//Actualiza el fitness de la solucion
	calcular_fitness();
	return;
}

/*----------------
Funcion de fitness
----------------*/
double c_cromosoma::calcular_fitness(){

	//Va a usar como funcion de fitness el indice Calinski-Harabasz 

	//Cuantos clusters hay habilitados
	clusters = 0;
	for(int k = 0; k < k_max; k++){
		//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
		if(valores_de_activacion[k] >= 0.5) clusters++;
	}
	//si solo hay un cluster, devuelve 0
	if(clusters < 2) return 0;
	
	//Va a sumar todas las distancias intra-clusters
	double distancia_intra_cluster = 0;

	//Asigna primero cada punto a un cluster por proximidad
	for(int punto = 0; punto < puntos->size(); punto++){
		//Inicializa la distancia al cluster al que pertenece (infinity = no pertenece a ninguno)		
		double distancia_a_k = INFINITY;
		for(int i=0; i < k_max; i++){
			//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
			if(valores_de_activacion[i] >= 0.5){
				//cout << (*puntos)[punto].x << "|" << (*puntos)[punto].y << "|" << centroides[i].x << "|" << centroides[i].y << "|" << ((*puntos)[punto].x - centroides[i].x) << "|" << ((*puntos)[punto].y - centroides[i].y) << "\n";
				//Si la distancia a este centroide es la mas chica, se lo asigno al centroide
				double distancia_a_i = (pow(((*puntos)[punto].x - centroides[i].x), 2) + pow(((*puntos)[punto].y - centroides[i].y), 2));
				if(distancia_a_i < distancia_a_k) {
					distancia_a_k = distancia_a_i;
				}
			}		
		}
		distancia_intra_cluster	+= distancia_a_k;//La distancia cuadriatica del punto al cluster que pertenece	
	}
		
	//Teniendo todo solo resta calcular el indice CH
	fitness = ((promedio_al_centro_del_dataset - distancia_intra_cluster)/distancia_intra_cluster)*((puntos->size()-clusters)/(clusters-1));
	return fitness;
}

/**********************
COMIENZAN LAS FUNCIONES
**********************/

/*-------------------------------------
Funcion para generar el dataset inicial
-------------------------------------*/
void generar_dataset(vector<c_punto> &centros, vector<c_punto> &puntos, c_rectangulo &valores_limite, c_generador_aleatorio &generador_aleatorio){

	//genera la red equiespaciada
	int centro = 0;	
	for(int x = 1; x <= centros_por_dimension; x++){
		for(int y = 1; y <= centros_por_dimension; y++){
			centros[centro].x = x;
			centros[centro].y = y;
			centro++;
		}
	}	

	//genera los puntos de datos alrededor de los centros y tambien busca los limites
	int punto = 0;
	for(centro = 0; centro < cantidad_de_centros; centro++){
		for(int punto_en_cluster = 0; punto_en_cluster < puntos_por_cluster; punto_en_cluster++){
			puntos[punto].x = centros[centro].x + generador_aleatorio.runif(-ancho_del_cluster, ancho_del_cluster);
			puntos[punto].y = centros[centro].y + generador_aleatorio.runif(-alto_del_cluster, alto_del_cluster);
			//Busca maximo y minimo en x
			if(puntos[punto].x > valores_limite.superior_derecho.x){
				valores_limite.superior_derecho.x = puntos[punto].x;
			}else if (puntos[punto].x < valores_limite.inferior_izquierdo.x){
				valores_limite.inferior_izquierdo.x = puntos[punto].x;			
			}
			//Busca maximo y minimo en y
			if(puntos[punto].y > valores_limite.superior_derecho.y){
				valores_limite.superior_derecho.y = puntos[punto].y;
			}else if (puntos[punto].y < valores_limite.inferior_izquierdo.y){
				valores_limite.inferior_izquierdo.y = puntos[punto].y;			
			}
			punto++;
		}		
	}	
	return;
}
/*-------------------------------------------------------
Funcion que encuentra el centroide para todos los puntos.
-------------------------------------------------------*/
c_punto encontrar_centroide(const vector<c_punto> &puntos){
	
	c_punto o_punto;
	
	for(int punto = 0; punto < puntos.size(); punto++){
		o_punto.x += puntos[punto].x;
		o_punto.y += puntos[punto].y;		
	}
	o_punto.x = o_punto.x/puntos.size();
	o_punto.y = o_punto.y/puntos.size();	
	return o_punto;
}

/*----------------------------------------------------
Funcion que encuentra el promedio de todos los puntos.
----------------------------------------------------*/
float promedio_dataset(const vector<c_punto> &puntos, c_punto centroide){
	float promedio = 0;
	for(int punto = 0; punto < puntos.size(); punto++){	
		promedio += ( pow((puntos[punto].x-centroide.x),2) + pow((puntos[punto].y-centroide.y),2) );
	}
	return promedio;
}

/*---------------------------------------------------
Funcion que elige una pareja en funcion de su fitness
---------------------------------------------------*/
void elegir_pareja(const vector<c_cromosoma> &cromosomas, vector<int> &pareja, c_generador_aleatorio &generador_aleatorio){
	
	//Inicializa un vector con los indices de cada cromosoma y lo shufflea
	vector<int> padres(poblacion);
	for(int padre = 0; padre < poblacion; padre++){
		padres[padre] = padre;
	}

	shuffle(padres.begin(), padres.end(), generador_aleatorio.generador_aleatorio);
	pareja[0]			= 0;
	pareja[1]			= 1;
	float mejor_fitness 		= -INFINITY;
	float segundo_mejor_fitness 	= -INFINITY;

	//Toma los primeros pp elementos del vector de indices de cromosomas y se queda con los cromosomas
	//que tengan mejor fitness y los guarda en pareja
	for(int padre = 0; padre < pp; padre++){
		if(cromosomas[padres[padre]].fitness > mejor_fitness){
			segundo_mejor_fitness = mejor_fitness;
			pareja[1] = pareja[0];
			mejor_fitness = cromosomas[padres[padre]].fitness;
			pareja[0] = padres[padre];
		}else if(cromosomas[padres[padre]].fitness > segundo_mejor_fitness){
			segundo_mejor_fitness = cromosomas[padres[padre]].fitness;
			pareja[1] = padres[padre];			
		}	
	}
	return;
	
}

/*-------------------
Funcion de crossover.
-------------------*/
void cruzar(const vector<c_cromosoma> &cromosomas_padres, vector<c_cromosoma> &cromosomas_hijos, c_generador_aleatorio &generador_aleatorio){

	//Vector auxiliar que va a contener la pareja seleccionada para cruza
	vector<int> pareja(2);

	//Elige la pareja a cruzar
	elegir_pareja(cromosomas_padres, pareja, generador_aleatorio);	

	//Genera los hijos como copias de los padres que pueden ser crossoveriados despues
	cromosomas_hijos.push_back(cromosomas_padres[pareja[0]]);
	cromosomas_hijos.push_back(cromosomas_padres[pareja[1]]); 	 
	//cout << cromosomas_hijos.size() << "\n" << cromosomas_padres[pareja[0]].to_string() << "\n" << cromosomas_padres[pareja[1]].to_string() << "\n" << cromosomas_hijos[cromosomas_hijos.size()-2].to_string() << "\n" << cromosomas_hijos[cromosomas_hijos.size()-1].to_string() << "\n";

	//Hace crossover entre dos cromosomas con probabilidad pc para cada elemento
	for(int k = 0; k < k_max; k++){
		if(generador_aleatorio.runif() < pc){ //intercambiamos los locus con probabilidad pc
			cromosomas_hijos[cromosomas_hijos.size()-2].centroides[k] = cromosomas_padres[pareja[1]].centroides[k];
			cromosomas_hijos[cromosomas_hijos.size()-1].centroides[k] = cromosomas_padres[pareja[0]].centroides[k];
		}
		
	}

	return;
}

/*------------------------------------------------
Calcula el promedio del fitness de las soluciones.
------------------------------------------------*/
float promedio_fitness(const vector<c_cromosoma> &cromosomas){
	float promedio = 0;
	for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
		promedio += cromosomas[cromosoma].fitness;	
	}
	return promedio/cromosomas.size();	
}

/*---------------------------------------------------
Comparador para hacer sorting manteniendo los indices
---------------------------------------------------*/
typedef std::pair<float,int> mypair;
bool comparator ( const mypair& l, const mypair& r){ return l.first > r.first; }

/*-----------------------------
Devuelve las mejores soluciones
-----------------------------*/
void reservar_soluciones_de_elite(const vector<c_cromosoma> &cromosomas, vector<c_cromosoma> &cromosomas_de_elite){

	//Armamos un vector de pares, fitness, indice, ordenamos por fitnes y devolvemos el 
	//indice de las soluciones_de_elite mejores
	vector<pair<float,int>> fitness(cromosomas.size());
	for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
		fitness[cromosoma].first  = cromosomas[cromosoma].fitness;
		fitness[cromosoma].second = cromosoma;
	}
	
	sort(fitness.begin(), fitness.end(), comparator);
	
	for(int cromosoma = 0; cromosoma < cromosomas_de_elite.size(); cromosoma++){
		cromosomas_de_elite[cromosoma] = cromosomas[fitness[cromosoma].second];	
	}
	
	return;	
}

/*-------------------------------------
Devuelve el indice de la mejor solucion
-------------------------------------*/
int mejor_solucion(const vector<c_cromosoma> &cromosomas){

	//Armamos un vector de pares, fitness, indice, ordenamos por fitnes y devolvemos el 
	//indice de las soluciones_de_elite mejores
	int idx_mejor_solucion 	= 0;
	float mejor_fitness 	= -INFINITY;
	for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
		if(cromosomas[cromosoma].fitness > mejor_fitness){
			idx_mejor_solucion = cromosoma;
			mejor_fitness = cromosomas[cromosoma].fitness;
		}
	}

	return idx_mejor_solucion;
}

/*--------------------------------------------------------------------
Pisa soluciones_de_elite cromosomas random con los cromosomas de elite
--------------------------------------------------------------------*/
void aplicar_soluciones_de_elite(vector<c_cromosoma> &cromosomas, vector<c_cromosoma> &cromosomas_de_elite, c_generador_aleatorio &generador_aleatorio){

	//Inicializa un vector con los indices de cada cromosoma y lo shufflea
	vector<int> indices_de_cromosomas(poblacion);
	for(int i = 0; i < cromosomas.size(); i++){
		indices_de_cromosomas[i] = i;
	}

	shuffle(indices_de_cromosomas.begin(), indices_de_cromosomas.end(), generador_aleatorio.generador_aleatorio);
	
	for(int cromosoma = 0; cromosoma < cromosomas_de_elite.size(); cromosoma++){
		cromosomas[indices_de_cromosomas[cromosoma]] = cromosomas_de_elite[cromosoma];	
	}
	
	return;	
}

/*---------------------------------
Funcion para mutar usando threading
---------------------------------*/
void mutar(vector<c_cromosoma> &cromosomas, int desde, int hasta){
	for(int cromosoma = desde; cromosoma < hasta; cromosoma++){
		cromosomas[cromosoma].mutar();
	}
}

/*******************
COMIENZA EL PROGRAMA
*******************/
int main(){

	//Archivo de salida
	ofstream archivo_soluciones;
	ofstream archivo_puntos;
	archivo_soluciones.open ("soluciones.txt");
	archivo_puntos.open ("puntos.txt");

	//Habilita el generador de numeros aleatorios para elegir los valores dentro de los cuales se van a mover los centroides
	c_generador_aleatorio generador_aleatorio;

	//Definimos los vectores que contienen el dataset
	vector<c_punto> centros(cantidad_de_centros);
	vector<c_punto> puntos(cantidad_de_puntos);
	c_rectangulo valores_limite(INFINITY, INFINITY, 0, 0);//Va a contener los valores limite que tomaron los puntos, en forma de un rectangulo (o sea, dos puntos, el angulo inferior izquierdo y el superior derecho)

	//Va a contener la nueva generacion de cromosomas despues de cada cruza
	vector<c_cromosoma> cromosomas_nuevos;
	
	//Vector de cromosomas de elite que van a pasar si o si a la proxima generacion
	vector<c_cromosoma> cromosomas_de_elite(soluciones_de_elite);
	
	//Genera el dataset
	generar_dataset(centros, puntos, valores_limite, generador_aleatorio);
	
	//Carga en el generador aleatorio las configuraciones limite para generar los distintos numeros
	generador_aleatorio.cargar_configuraciones_de_cromosomas(pm, alfa, epsilon, k_max, &valores_limite); 

	//Para el fitness de los cromosomas vamos a necesitar el promedio 
	//de la distancia de los puntos al centroide. Calculamos esto una sola vez.
	c_punto centroide = encontrar_centroide(puntos);
	float promedio = promedio_dataset(puntos, centroide);
	
	//Genera la poblacion inicial de cromosomas
	//Cada cromosoma tiene k_max valores de activacion y k_max centroides
	vector<c_cromosoma> cromosomas;
	for(int cromosoma = 0; cromosoma < poblacion; cromosoma++){
		cromosomas.push_back(c_cromosoma(k_max, &valores_limite, promedio, &puntos, &generador_aleatorio));
	}

	//Construye el cromosoma objetivo para ver el fitness al que apuntamos (ground truth)
	vector<float> vda(k_max, 0);
	for(int k = 0; k < cantidad_de_centros; k++){
		vda[k] = 0.5;
	}
	vector<c_punto> centroides_objetivo(k_max, c_punto());
	for(int k = 0; k < cantidad_de_centros; k++){
		centroides_objetivo[k] = centros[k];
	}
	c_cromosoma cromosoma_objetivo(k_max, &valores_limite, promedio, &puntos, &generador_aleatorio, vda, centroides_objetivo); 	

	//Comienza la evolucion
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	for(int generacion = 0; generacion < generaciones; generacion++){
	
		//Se queda con los mejores soluciones_de_elite
		reservar_soluciones_de_elite(cromosomas, cromosomas_de_elite);
		
		//Cruza los cromosomas de acuerdo a su fitness. Cuanto mas fitness mas probabilidad de cruza. 
		//Elige poblacion/2 parejas. Vacia previamente el vector de cromosomas nuevos
		cromosomas_nuevos.clear();

		//Genera la poblacion nueva
		for(int cruza = 0; cruza < (poblacion/2); cruza++){
			cruzar(cromosomas, cromosomas_nuevos, generador_aleatorio);
		}
			
		//Muta la poblacion (eso calcula ademas el nuevo fitness)
		if(multithreding){
			//Muta la poblacion nueva usando dos threads para paralelismo. Los argumentos que 
			//queremos que se modifiquen tienen que ir entre std::ref()
			std::thread primera(mutar, ref(cromosomas_nuevos), 0, cromosomas_nuevos.size()/2);
			std::thread segunda(mutar, ref(cromosomas_nuevos), cromosomas_nuevos.size()/2, cromosomas_nuevos.size());  
		
			//Pausa hasta que terminan los dos threads
			primera.join();
			segunda.join();
		}else{
			for(int cromosoma = 0; cromosoma < poblacion; cromosoma++){
				cromosomas_nuevos[cromosoma].mutar();	
			}
		}				
		//Actualiza la poblacion con la poblacion nueva
		cromosomas = cromosomas_nuevos;

		//Pisa soluciones al azar con los cromosomas de elite	
		aplicar_soluciones_de_elite(cromosomas, cromosomas_de_elite, generador_aleatorio);
	  	
	  	//Guarda en archivo el promedio del fitness
	  	if(generacion % 100 == 0){
		  	if(generacion == 100){
		  		std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
		  		auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
				double duracion_corrida = (duration/(100*1000000.0L))*generaciones;
				cout << "Duracion estimada de corrida: " << duracion_corrida << "\n";
		  	}
		  	int idx_mejor_solucion = mejor_solucion(cromosomas);
			cout << "mejor solucion|mejor fitness|clusters: " << idx_mejor_solucion << "|" << cromosomas[idx_mejor_solucion].fitness << "|" << cromosomas[idx_mejor_solucion].clusters << "\n";
			cout << "Promedio/Optima (paso): " << promedio_fitness(cromosomas) << "/" << cromosoma_objetivo.fitness << " (" << generacion << "/" << generaciones << ")\n";
		}

	}

  	int idx_mejor_solucion = mejor_solucion(cromosomas);
	cout << "mejor solucion|mejor fitness|clusters: " << idx_mejor_solucion << "|" << cromosomas[idx_mejor_solucion].fitness << "|" << cromosomas[idx_mejor_solucion].clusters << "\n";
	cout << "Promedio/Optima  (paso): " << promedio_fitness(cromosomas) << "/" << cromosoma_objetivo.fitness << " (" << generaciones << "/" << generaciones << ")\n";	
	cout << cromosomas[idx_mejor_solucion].to_string() << "\n";
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Duracion real: " << duration/1000000.0L << "\n";

 	for(int c = 0; c < cromosomas.size(); c++){
		archivo_soluciones << cromosomas[c].to_cluster() << "\n";
	}	
 	archivo_soluciones.close();
 	
 	for(int p = 0; p < puntos.size(); p++){
		archivo_puntos << puntos[p].x << "\t" << puntos[p].y << "\n";
	}
 	archivo_puntos.close();
	
	return 0;
}


