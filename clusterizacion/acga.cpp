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
#include <stdlib.h>
#include <unistd.h>

using namespace std;

/**************
CONFIGURACIONES
**************/

/*-------------------------
Configuraciones del dataset
-------------------------*/
int centros_por_dimension = 2; //los centros son los lugares alrededor de los cuales se van a armar los clusters
int puntos_por_cluster = 10;
double parametro_de_red = 1;
double ancho_del_cluster = 0.1; //lo que mide el cluster en x
double alto_del_cluster = 0.1; //lo que mide el cluster en y

/*-------------------
Configuraciones de AG
-------------------*/
int especies = 1; //Una especie esta compuesta por toda una poblacion y las especies se van evolucionando en paralelo
int poblacion = 40;
double pm = 0.1; //probabilidad de mutacion
double pc = 0.3; //probabilidad de single-point crossover
int pp = 3; //Cromosomas a seleccionar aleatoreamente en cada busqueda de padres
int k_max = 20; //Maxima cantidad de clusters a buscar
double alfa = 0.1; //Amplitud de mutacion de valores de activacion
double epsilon = 0.3; //Amplitud de mutacion de centroides
int soluciones_de_elite = 1; //Las mejores soluciones pasan sin alteraciones a la proxima generacion
int generaciones = 3000;
int cada_cuantas_generaciones_intercambiar_especies = 500;
int soluciones_de_elite_a_intercambiar = 4;
int cada_cuantas_generaciones_splitear_clusters = 250;

/*--------------------------
Configuraciones del programa
--------------------------*/
bool multithreding 				= true;
string dataset 					= "s3"; //Dataset a clusterizar. Si es vacio genera una red al azar
int cutoff					= 500; //Cuantas iteraciones tienen que pasar con el mismo max fitness para que corte
int cantidad_de_clusters_inicial 		= 3; //Fuerza la cantidad de clusters habilitados que tiene que tener una solucion en la primer generacion como minimo
int cada_cuantas_generaciones_mostrar_datos 	= 100;
unsigned int cantidad_de_threads = std::thread::hardware_concurrency(); //std::thread::hardware_concurrency() autodetecta el maximo de procesadores

/*--------------------------------
Calcula algunas cantidades previas
--------------------------------*/
int cantidad_de_centros = centros_por_dimension*centros_por_dimension;
int cantidad_de_puntos  = cantidad_de_centros*puntos_por_cluster;

/*******************
COMIENZAN LAS CLASES
*******************/

/*---------
Clase punto
---------*/
class c_punto{
	public:
	double x;
	double y;
	
	c_punto();
	c_punto(double x, double y);
};
//Constructores de punto
c_punto::c_punto(){
	x = 0;
	y = 0;
}	

c_punto::c_punto(double x, double y){
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
	
	c_rectangulo(double x1, double y1, double x2, double y2);
};
//Constructor de rectangulo
c_rectangulo::c_rectangulo(double x1, double y1, double x2, double y2){
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
	double pm;
	double alfa;
	double epsilon_x;
	double epsilon_y;
	c_rectangulo* valores_limite;
	c_generador_aleatorio();
	c_generador_aleatorio(int seed);
	void cargar_configuraciones_de_cromosomas(double pm, double alfa, double epsilon, int k_max, c_rectangulo* valores_limite);
	
	uniform_real_distribution<double> urd;
	uniform_real_distribution<double> urd_x;
	uniform_real_distribution<double> urd_y;
	uniform_real_distribution<double> urd_alfa;
	uniform_real_distribution<double> urd_epsilon_x;
	uniform_real_distribution<double> urd_epsilon_y;
	
	double runif(); //Llama a runif(0, 1);
	double runif(double limite_inferior, double limite_superior);	
	
	//Todas estas funciones generan reales aleatorios entre los limites establecidos al cargar la configuracion de los cromosomas
	double runif_x();
	double runif_y();
	double runif_alfa();
	double runif_epsilon_x();
	double runif_epsilon_y();
};

c_generador_aleatorio::c_generador_aleatorio(){
	generador_aleatorio.seed(std::random_device{}());
}

c_generador_aleatorio::c_generador_aleatorio(int seed){
	generador_aleatorio.seed(seed);
}

void c_generador_aleatorio::cargar_configuraciones_de_cromosomas(double pm, double alfa, double epsilon, int k_max, c_rectangulo* valores_limite){
	c_generador_aleatorio::generador_aleatorio = generador_aleatorio;
	c_generador_aleatorio::valores_limite = valores_limite;
	c_generador_aleatorio::pm = pm;
	c_generador_aleatorio::alfa = alfa;
	c_generador_aleatorio::epsilon_x = epsilon*(valores_limite->superior_derecho.x-valores_limite->inferior_izquierdo.x);
	c_generador_aleatorio::epsilon_y = epsilon*(valores_limite->superior_derecho.y-valores_limite->inferior_izquierdo.y);	
	urd 		= uniform_real_distribution<double>(0, 1);
	urd_x 		= uniform_real_distribution<double>(valores_limite->inferior_izquierdo.x, valores_limite->superior_derecho.x);
	urd_y		= uniform_real_distribution<double>(valores_limite->inferior_izquierdo.y, valores_limite->superior_derecho.y);
	urd_alfa 	= uniform_real_distribution<double>(-alfa, alfa);
	urd_epsilon_x 	= uniform_real_distribution<double>(-epsilon_x, epsilon_x);
	urd_epsilon_y 	= uniform_real_distribution<double>(-epsilon_y, epsilon_y);	
	return;
}

double c_generador_aleatorio::runif(){
	return urd(generador_aleatorio);
}
double c_generador_aleatorio::runif(double limite_inferior, double limite_superior){
	uniform_real_distribution<double> urd_especial(limite_inferior, limite_superior);
	return urd_especial(generador_aleatorio);
}	

//Todas estas funciones generan reales aleatorios entre los limites establecidos al cargar la configuracion de los cromosomas
double c_generador_aleatorio::runif_x(){
	return urd_x(generador_aleatorio);
}
double c_generador_aleatorio::runif_y(){
	return urd_y(generador_aleatorio);
}
double c_generador_aleatorio::runif_alfa(){
	return urd_alfa(generador_aleatorio);
}
double c_generador_aleatorio::runif_epsilon_x(){
	return urd_epsilon_x(generador_aleatorio);
}
double c_generador_aleatorio::runif_epsilon_y(){
	return urd_epsilon_y(generador_aleatorio);
}

/*-------------
Clase cromosoma
-------------*/
class c_cromosoma {
	public:
	int k_max;
	int clusters;
	vector<double> valores_de_activacion;
	vector<c_punto> centroides;
	double promedio_al_centro_del_dataset;
	vector<c_punto*>* puntos;
	vector<int> puntos_clusterizados;
	double fitness;
	c_generador_aleatorio* generador_aleatorio;
	c_rectangulo* valores_limite;
	vector<int> solucion_formato_cluster; //Representacion de la solucion en formato cluster(a que cluster va cada punto)
	bool primera_corrida;
		
	c_cromosoma();
	c_cromosoma(int k_max, c_rectangulo* valores_limite, double promedio_al_centro_del_dataset, vector<c_punto*> *puntos, c_generador_aleatorio* generador_aleatorio);
	c_cromosoma(int k_max, c_rectangulo* valores_limite, double promedio_al_centro_del_dataset, vector<c_punto*> *puntos, c_generador_aleatorio* generador_aleatorio, const vector<double> &valores_de_activacion, const vector<c_punto> &centroides);
	double calcular_fitness();
	void mutar();
	string to_string();
	string to_cluster();
	void splitear_el_cluster_mas_grande();
};

//El constructor por defecto viene vacio. Sirve para definir cromosomas temporales
c_cromosoma::c_cromosoma(){}

//Este constructor genera las soluciones de forma random. Es el constructor para generar la tirada inicial de cromosomas
c_cromosoma::c_cromosoma(int k_max, c_rectangulo* valores_limite, double promedio_al_centro_del_dataset, vector<c_punto*> *puntos, c_generador_aleatorio* generador_aleatorio){
	
	//Inicializa las propiedades del cluster
	c_cromosoma::k_max = k_max;
	c_cromosoma::promedio_al_centro_del_dataset = promedio_al_centro_del_dataset;
	c_cromosoma::puntos = puntos;
	c_cromosoma::valores_limite = valores_limite;
	c_cromosoma::generador_aleatorio = generador_aleatorio;
	puntos_clusterizados = vector<int>(puntos->size());
	
	//Genera el cromosoma al azar de la corrida. Primero los valores de activacion.
	//Fuerza que haya cantidad_de_clusters_inicial habilitados
	for(int cluster_inicial = 0; cluster_inicial < cantidad_de_clusters_inicial; cluster_inicial++){
		valores_de_activacion.push_back(0.5);
	}

	for(int valor = cantidad_de_clusters_inicial; valor < k_max; valor++){
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
c_cromosoma::c_cromosoma(int k_max, c_rectangulo* valores_limite, double promedio_al_centro_del_dataset, vector<c_punto*> *puntos, c_generador_aleatorio* generador_aleatorio, const vector<double> &valores_de_activacion, const vector<c_punto> &centroides){
	
	c_cromosoma::k_max = k_max;
	c_cromosoma::promedio_al_centro_del_dataset = promedio_al_centro_del_dataset;
	c_cromosoma::puntos = puntos;
	c_cromosoma::valores_limite = valores_limite;
	puntos_clusterizados = vector<int>(puntos->size());
	
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

	//Cuantos clusters hay habilitados
	vector<int> clusters_habilitados;
	for(int k = 0; k < k_max; k++){
		//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
		if(valores_de_activacion[k] >= 0.5) clusters_habilitados.push_back(k);
	}

	//Asigna cada punto a un cluster por proximidad
	for(int punto = 0; punto < puntos->size(); punto++){
		//Inicializa la distancia al cluster al que pertenece (infinity = no pertenece a ninguno)		
		double distancia_a_k = INFINITY;
		int k;
		for(int i = 0; i < clusters_habilitados.size(); i++){
			//Si la distancia a este centroide es la mas chica, se lo asigno al centroide
			double distancia_a_i = (pow(((*puntos)[punto]->x - centroides[clusters_habilitados[i]].x), 2) + pow(((*puntos)[punto]->y - centroides[clusters_habilitados[i]].y), 2));
			if(distancia_a_i < distancia_a_k) {
				distancia_a_k = distancia_a_i;
				k = i;
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
	for(int k = 0; k < (k_max-1); k++){
		buffer << centroides[k].x << "," << centroides[k].y << ",";
	}
	buffer << centroides[(k_max-1)].x << "," << centroides[(k_max-1)].y;	

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
			double valor_viejo = valores_de_activacion[k];
			double nuevo_va = valores_de_activacion[k] + generador_aleatorio->runif_alfa();
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
	
			double nuevo_x = centroides[k].x + generador_aleatorio->runif_epsilon_x();
			double nuevo_y = centroides[k].y + generador_aleatorio->runif_epsilon_y();
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
	vector<int> clusters_habilitados;
	for(int k = 0; k < k_max; k++){
		//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
		if(valores_de_activacion[k] >= 0.5) clusters_habilitados.push_back(k);
	}

	clusters = clusters_habilitados.size();
	//si solo hay un cluster, devuelve 0.
	if(clusters < 2) return 0;

	//Va a sumar todas las distancias intra-clusters
	double distancia_intra_cluster = 0;
		
	//Asigna primero cada punto a un cluster por proximidad
	for(int punto = 0; punto < puntos->size(); punto++){
		//Inicializa la distancia al cluster al que pertenece (infinity = no pertenece a ninguno)		
		double distancia_a_k = INFINITY;
		for(int i=0; i < clusters_habilitados.size(); i++){
			//Si la distancia a este centroide es la mas chica, se lo asigno al centroide
			double distancia_a_i = (pow(((*puntos)[punto]->x - centroides[clusters_habilitados[i]].x), 2) + pow(((*puntos)[punto]->y - centroides[clusters_habilitados[i]].y), 2));
			if(distancia_a_i < distancia_a_k) {
				distancia_a_k = distancia_a_i;
				puntos_clusterizados[punto] = clusters_habilitados[i];
			}
		}
		distancia_intra_cluster	+= distancia_a_k;//La distancia cuadriatica del punto al cluster que pertenece	
	}

	//Teniendo todo solo resta calcular el indice CH
	fitness = ((promedio_al_centro_del_dataset - distancia_intra_cluster)/distancia_intra_cluster)*((puntos->size()-clusters)/(clusters-1));
	return fitness;
}

/*----------------
Funcion de spliting
----------------*/
void c_cromosoma::splitear_el_cluster_mas_grande(){

	//Cuantos clusters hay habilitados
	vector<int> clusters_habilitados;
	for(int k = 0; k < k_max; k++){
		//Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
		if(valores_de_activacion[k] >= 0.5) clusters_habilitados.push_back(k);
	}

	clusters = clusters_habilitados.size();
	
	//Si la cantidad de clusters es la maxima, sale
	if(clusters == centroides.size()) return;

	//Va a sumar todas las distancias intra-clusters
	vector<double> distancias_intra_clusters(centroides.size());
	vector<int> puntos_en_clusters(puntos->size());
	double distancia_intra_cluster = 0;
		
	//Asigna primero cada punto a un cluster por proximidad
	for(int punto = 0; punto < puntos->size(); punto++){
		//Inicializa la distancia al cluster al que pertenece (infinity = no pertenece a ninguno)		
		double distancia_a_k = INFINITY;
		for(int i=0; i < clusters_habilitados.size(); i++){
			//Si la distancia a este centroide es la mas chica, se lo asigno al centroide
			double distancia_a_i = (pow(((*puntos)[punto]->x - centroides[clusters_habilitados[i]].x), 2) + pow(((*puntos)[punto]->y - centroides[clusters_habilitados[i]].y), 2));
			if(distancia_a_i < distancia_a_k) {
				distancia_a_k = distancia_a_i;
				puntos_en_clusters[punto] = clusters_habilitados[i];
			}
		}
		distancias_intra_clusters[puntos_en_clusters[punto]] += distancia_a_k;
		distancia_intra_cluster	+= distancia_a_k;//La distancia cuadriatica del punto al cluster que pertenece	
	}

	//buscamos el cluster que mas aporta a las distancias y lo spliteamos por algun lado
	int idx_mayor_aporte 	= 0;
	double mayor_aporte 	= -INFINITY;
	for(int i=0; i < clusters_habilitados.size(); i++){
		if(distancias_intra_clusters[clusters_habilitados[i]] > mayor_aporte){
			idx_mayor_aporte = clusters_habilitados[i];
			mayor_aporte = distancias_intra_clusters[clusters_habilitados[i]];
		}
	}	
	
	//buscamos el ancho y el alto del cluster ordenando los puntos
	c_rectangulo rectangulo(INFINITY, INFINITY, -INFINITY, -INFINITY);
	for(int punto = 0; punto < puntos->size(); punto++){
		if(puntos_en_clusters[punto] == idx_mayor_aporte){
			if((*puntos)[punto]->x < rectangulo.inferior_izquierdo.x){
				rectangulo.inferior_izquierdo.x = (*puntos)[punto]->x;
			}else if((*puntos)[punto]->x > rectangulo.superior_derecho.x){
				rectangulo.superior_derecho.x = (*puntos)[punto]->x;
			}

			if((*puntos)[punto]->y < rectangulo.inferior_izquierdo.y){
				rectangulo.inferior_izquierdo.y = (*puntos)[punto]->y;
			}else if((*puntos)[punto]->y > rectangulo.superior_derecho.y){
				rectangulo.superior_derecho.y = (*puntos)[punto]->y;
			}			
		
		}
	}

	double ancho = rectangulo.superior_derecho.x - rectangulo.inferior_izquierdo.x;
	double alto  = rectangulo.superior_derecho.y - rectangulo.inferior_izquierdo.y;
	
	//generamos dos centroides mas a partir del centroide actual y de la longitud caracteristica
	c_punto centroide_nuevo_1;
	c_punto centroide_nuevo_2;
	c_punto centroide_viejo_1;
	c_punto centroide_viejo_2;	
	if(ancho > alto){
		centroide_nuevo_1.x = centroides[idx_mayor_aporte].x - (ancho/3);
		centroide_nuevo_1.y = centroides[idx_mayor_aporte].y;
		centroide_nuevo_2.x = centroides[idx_mayor_aporte].x + (ancho/3);
		centroide_nuevo_2.y = centroides[idx_mayor_aporte].y;
	}else{
		centroide_nuevo_1.x = centroides[idx_mayor_aporte].x ;
		centroide_nuevo_1.y = centroides[idx_mayor_aporte].y - (alto/3);
		centroide_nuevo_2.x = centroides[idx_mayor_aporte].x;
		centroide_nuevo_2.y = centroides[idx_mayor_aporte].y + (alto/3);
	}

	centroide_viejo_1 = centroides[idx_mayor_aporte];

	//Reemplazamos el centroide actual y alguno apagado por los nuevos
	int idx_viejo;
	double valor_de_activacion_viejo;
	for(int k = 0; k < k_max; k++){
		//Si el valor de activacion es menor a 0.5 (i.e. esta apagado) lo usamos
		if(valores_de_activacion[k] < 0.5){
			valor_de_activacion_viejo = valores_de_activacion[k];
			valores_de_activacion[k] = 0.5;
			centroide_viejo_2 = centroides[k];
			centroides[k] = centroide_nuevo_1;
			idx_viejo = k;
			break;
		}
	}
	centroides[idx_mayor_aporte] = centroide_nuevo_2;
	
	//Calculamos el delta fitness
	double delta_fitness = fitness - calcular_fitness();

	//Montecarleamos con el fitness nuevo (descartamos los cambios con probabilidad 1-exp{-Bdf})
	double r = generador_aleatorio->runif();
	if(r > exp(-delta_fitness)){
		
		valores_de_activacion[idx_viejo] = valor_de_activacion_viejo;
		centroides[idx_viejo] 		 = centroide_viejo_2;
		centroides[idx_mayor_aporte]	 = centroide_viejo_1;
	}
		
}

/*---------------------
Clase fitness landscape
---------------------*/
class c_fitness_landscape{
	public:
	vector<vector<double>> disimilaridad;
	c_fitness_landscape(int cromosomas);
	void calcular_disimilaridad(const vector<c_cromosoma> &cromosomas, const vector<c_punto> &puntos);
	
};

c_fitness_landscape::c_fitness_landscape(int cromosomas){
	disimilaridad = vector<vector<double>> (cromosomas, vector<double>(cromosomas));
}

void c_fitness_landscape::calcular_disimilaridad(const vector<c_cromosoma> &cromosomas, const vector<c_punto> &puntos){
	//Usamos jaccard y calculamos la disimilaridad de forma pairwise
	int M01 = 0; //Numero de atributos en A pero no en B
	int M10 = 0; //Numero de atributos en B pero no en A
	int M11 = 0; //Numero de atributos en ambos a la vez
	
	for(int cromosoma_1 = 0; cromosoma_1 < cromosomas.size(); cromosoma_1++){
		M01 = M10 = M11 = 0;
		for(int cromosoma_2 = (cromosoma_1 + 1); cromosoma_2 < cromosomas.size(); cromosoma_2++){
			for(int punto = 0; punto < puntos.size(); punto++){
				//if(cromosomas[cromosoma_1]
			}
		}
		
	}
	
	
}

/**********************
COMIENZAN LAS FUNCIONES
**********************/

/*-------------------------------------
Funcion para generar el dataset inicial
-------------------------------------*/
void generar_dataset(vector<c_punto> &centros, vector<c_punto*> &puntos, c_rectangulo &valores_limite, c_generador_aleatorio &generador_aleatorio){

	//genera la red equiespaciada
	int centro = 0;	
	for(int x = 1; x <= centros_por_dimension; x++){
		for(int y = 1; y <= centros_por_dimension; y++){
			centros.push_back(c_punto());
			centros[centro].x = x;
			centros[centro].y = y;
			centro++;
		}
	}	

	//genera los puntos de datos alrededor de los centros y tambien busca los limites
	int punto = 0;
	for(centro = 0; centro < cantidad_de_centros; centro++){
		for(int punto_en_cluster = 0; punto_en_cluster < puntos_por_cluster; punto_en_cluster++){
			puntos.push_back(new c_punto());
			puntos[punto]->x = centros[centro].x + generador_aleatorio.runif(-ancho_del_cluster, ancho_del_cluster);
			puntos[punto]->y = centros[centro].y + generador_aleatorio.runif(-alto_del_cluster, alto_del_cluster);
			//Busca maximo y minimo en x
			if(puntos[punto]->x > valores_limite.superior_derecho.x){
				valores_limite.superior_derecho.x = puntos[punto]->x;
			}else if (puntos[punto]->x < valores_limite.inferior_izquierdo.x){
				valores_limite.inferior_izquierdo.x = puntos[punto]->x;			
			}
			//Busca maximo y minimo en y
			if(puntos[punto]->y > valores_limite.superior_derecho.y){
				valores_limite.superior_derecho.y = puntos[punto]->y;
			}else if (puntos[punto]->y < valores_limite.inferior_izquierdo.y){
				valores_limite.inferior_izquierdo.y = puntos[punto]->y;			
			}
			punto++;
		}		
	}	
	return;
}

/*------------------------------------
Funcion para cargar el dataset inicial
------------------------------------*/
void cargar_dataset(vector<c_punto> &centros, vector<c_punto*> &puntos, c_rectangulo &valores_limite, string dataset){

	stringstream buffer;
	buffer << "datasets/" << dataset << "-cb.txt";
	ifstream archivo;
	
	//Lee el ground truth y carga la red
	archivo.open(buffer.str().c_str());
	
	double x;
	double y;

	cantidad_de_centros = 0;
	while (archivo >> x >> y)
	{
		centros.push_back(c_punto());
		centros[cantidad_de_centros].x = x;
		centros[cantidad_de_centros].y = y;
		cantidad_de_centros++;
	}	

	archivo.close();
	//Lee el dataset y carga los puntos. Además busca los limites
	buffer.str("");
	buffer << "datasets/" << dataset << ".txt";
	
	archivo.open(buffer.str().c_str());	
	cantidad_de_puntos = 0;
	while (archivo >> x >> y){
		puntos.push_back(new c_punto());
		puntos[cantidad_de_puntos]->x = x;
		puntos[cantidad_de_puntos]->y = y;
		//Busca maximo y minimo en x
		if(puntos[cantidad_de_puntos]->x > valores_limite.superior_derecho.x){
			valores_limite.superior_derecho.x = puntos[cantidad_de_puntos]->x;
		}else if (puntos[cantidad_de_puntos]->x < valores_limite.inferior_izquierdo.x){
			valores_limite.inferior_izquierdo.x = puntos[cantidad_de_puntos]->x;			
		}
		//Busca maximo y minimo en y
		if(puntos[cantidad_de_puntos]->y > valores_limite.superior_derecho.y){
			valores_limite.superior_derecho.y = puntos[cantidad_de_puntos]->y;
		}else if (puntos[cantidad_de_puntos]->y < valores_limite.inferior_izquierdo.y){
			valores_limite.inferior_izquierdo.y = puntos[cantidad_de_puntos]->y;			
		}
		cantidad_de_puntos++;
	}	
	archivo.close();	
	return;
}

/*-------------------------------------------------------
Funcion que encuentra el centroide para todos los puntos.
-------------------------------------------------------*/
c_punto encontrar_centroide(const vector<c_punto*> &puntos){
	
	c_punto o_punto;
	
	for(int punto = 0; punto < puntos.size(); punto++){
		o_punto.x += puntos[punto]->x;
		o_punto.y += puntos[punto]->y;		
	}
	o_punto.x = o_punto.x/puntos.size();
	o_punto.y = o_punto.y/puntos.size();	
	return o_punto;
}

/*----------------------------------------------------
Funcion que encuentra el promedio de todos los puntos.
----------------------------------------------------*/
double promedio_dataset(const vector<c_punto*> &puntos, c_punto centroide){
	double promedio = 0;
	for(int punto = 0; punto < puntos.size(); punto++){	
		promedio += ( pow((puntos[punto]->x-centroide.x),2) + pow((puntos[punto]->y-centroide.y),2) );
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
	double mejor_fitness 		= -INFINITY;
	double segundo_mejor_fitness 	= -INFINITY;

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
double promedio_fitness(const vector<c_cromosoma> &cromosomas){
	double promedio = 0;
	for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
		promedio += cromosomas[cromosoma].fitness;	
	}
	return promedio/cromosomas.size();	
}

/*---------------------------------------------------
Comparador para hacer sorting manteniendo los indices
---------------------------------------------------*/
typedef std::pair<double,int> mypair;
bool comparator ( const mypair& l, const mypair& r){ return l.first > r.first; }

/*-----------------------------
Devuelve las mejores soluciones
-----------------------------*/
void reservar_soluciones_de_elite(const vector<c_cromosoma> &cromosomas, vector<c_cromosoma> &cromosomas_de_elite, int soluciones_de_elite){

	//Armamos un vector de pares, fitness, indice, ordenamos por fitnes y devolvemos el 
	//indice de las soluciones_de_elite mejores
	vector<pair<double,int>> fitness(cromosomas.size());
	for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
		fitness[cromosoma].first  = cromosomas[cromosoma].fitness;
		fitness[cromosoma].second = cromosoma;
	}
	
	sort(fitness.begin(), fitness.end(), comparator);
	cromosomas_de_elite.clear();
	for(int cromosoma = 0; cromosoma < soluciones_de_elite; cromosoma++){
		cromosomas_de_elite.push_back(cromosomas[fitness[cromosoma].second]);	
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
	double mejor_fitness 	= -INFINITY;
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
	ofstream archivo_mejores_soluciones;
	//ofstream archivo_mejores_soluciones_formato_cluster;
	ofstream archivo_puntos;
	ofstream archivo_soluciones;
	
	archivo_soluciones.open ("soluciones.txt");
	archivo_puntos.open ("puntos.txt");
	archivo_mejores_soluciones.open ("mejores_soluciones.txt");
	//archivo_mejores_soluciones_formato_cluster.open ("mejores_soluciones_formato_cluster.txt");
	
	//Buffer donde vamos a guardar las soluciones en cada paso para no grabarlas a disco que tarda mucho
	stringstream buffer_de_soluciones;
	
	//Habilita el generador de numeros aleatorios para elegir los valores dentro de los cuales se van a mover los centroides
	c_generador_aleatorio generador_aleatorio;

	//Definimos los vectores que contienen el dataset
	vector<c_punto> centros;
	vector<c_punto*> puntos;
	c_rectangulo valores_limite(INFINITY, INFINITY, 0, 0);//Va a contener los valores limite que tomaron los puntos, en forma de un rectangulo (o sea, dos puntos, el angulo inferior izquierdo y el superior derecho)

	vector<vector<c_cromosoma>> cromosomas_nuevos(especies);
	//Va a contener la nueva generacion de cromosomas despues de cada cruza
	//vector<c_cromosoma> cromosomas_nuevos;
	
	//Vector de cromosomas de elite que van a pasar si o si a la proxima generacion
	vector<vector<c_cromosoma>> cromosomas_de_elite(especies);
	//vector<c_cromosoma> cromosomas_de_elite(soluciones_de_elite);

	//Genera el dataset
	if(dataset.size()){
		cargar_dataset(centros, puntos, valores_limite, dataset);
	}else{
		generar_dataset(centros, puntos, valores_limite, generador_aleatorio);
	}
	
	//Carga en el generador aleatorio las configuraciones limite para generar los distintos numeros
	generador_aleatorio.cargar_configuraciones_de_cromosomas(pm, alfa, epsilon, k_max, &valores_limite); 

	//Para el fitness de los cromosomas vamos a necesitar el promedio 
	//de la distancia de los puntos al centroide. Calculamos esto una sola vez.
	c_punto centroide = encontrar_centroide(puntos);
	double promedio = promedio_dataset(puntos, centroide);
	
	//Genera la poblacion inicial de cromosomas
	//Cada cromosoma tiene k_max valores de activacion y k_max centroides
	vector<vector<c_cromosoma>> cromosomas;
	for(int especie = 0; especie < especies; especie++){
		vector<c_cromosoma> cromosomas_auxiliar;
		for(int cromosoma = 0; cromosoma < poblacion; cromosoma++){
			cromosomas_auxiliar.push_back(c_cromosoma(k_max, &valores_limite, promedio, &puntos, &generador_aleatorio));
		}
		cromosomas.push_back(cromosomas_auxiliar);
		cromosomas_auxiliar.clear();
	}
	//vector<c_cromosoma> cromosomas;
	//for(int cromosoma = 0; cromosoma < poblacion; cromosoma++){
	//	cromosomas.push_back(c_cromosoma(k_max, &valores_limite, promedio, &puntos, &generador_aleatorio));
	//}
	
	//Construye el cromosoma objetivo para ver el fitness al que apuntamos (ground truth)
	vector<double> vda(k_max, 0);
	for(int k = 0; k < cantidad_de_centros; k++){
		vda[k] = 0.5;
	}
	vector<c_punto> centroides_objetivo(k_max, c_punto());
	for(int k = 0; k < cantidad_de_centros; k++){
		centroides_objetivo[k] = centros[k];
	}
	c_cromosoma cromosoma_objetivo(k_max, &valores_limite, promedio, &puntos, &generador_aleatorio, vda, centroides_objetivo); 	


	//Calcula cuantos threads va a usar (usa la maxima cantidad de procesadores disponibles)
	if(multithreding && cantidad_de_threads == 0){
		cantidad_de_threads = sysconf(_SC_NPROCESSORS_ONLN);
	}
	std::thread threads[cantidad_de_threads];
	cout << "Usando " << cantidad_de_threads << " threads\n";
	//Comienza la evolucion
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	int iteraciones_con_max_fitness_igual[especies];
	double max_fitness[especies];
	int generacion;
	for(generacion = 0; generacion < generaciones; generacion++){
		for(int especie = 0; especie < especies; especie++){
			
			//Se queda con los mejores soluciones_de_elite
			reservar_soluciones_de_elite(cromosomas[especie], cromosomas_de_elite[especie], soluciones_de_elite);
			
			//Cruza los cromosomas de acuerdo a su fitness. Cuanto mas fitness mas probabilidad de cruza. 
			//Elige poblacion/2 parejas. Vacia previamente el vector de cromosomas nuevos
			cromosomas_nuevos[especie].clear();

			//Genera la poblacion nueva
			for(int cruza = 0; cruza < (poblacion/2); cruza++){
				cruzar(cromosomas[especie], cromosomas_nuevos[especie], generador_aleatorio);
			}
			cromosomas_nuevos[especie] = cromosomas[especie];
				
			//Muta la poblacion (eso calcula ademas el nuevo fitness)
			if(multithreding){
				//Muta la poblacion nueva usando dos threads para paralelismo. Los argumentos que 
				//queremos que se modifiquen tienen que ir entre std::ref()
				for(int n_thread = 0; n_thread < cantidad_de_threads; n_thread++){
					
					//Crea el thread
					threads[n_thread] = std::thread(mutar, ref(cromosomas_nuevos[especie]), (cromosomas_nuevos[especie].size()/cantidad_de_threads)*n_thread, (cromosomas_nuevos[especie].size()/cantidad_de_threads)*(n_thread+1));
				}
				for(int n_thread = 0; n_thread < cantidad_de_threads; n_thread++){
					//Pausa hasta que terminan los threads
					threads[n_thread].join();
				}
			}else{
				mutar(cromosomas_nuevos[especie], 0, cromosomas_nuevos[especie].size());
			}	
				
			//Actualiza la poblacion con la poblacion nueva
			cromosomas[especie] = cromosomas_nuevos[especie];

			//Pisa soluciones al azar con los cromosomas de elite	
			aplicar_soluciones_de_elite(cromosomas[especie], cromosomas_de_elite[especie], generador_aleatorio);
			
			//Guarda todas las soluciones de la generacion
			for(int cromosoma = 0; cromosoma < cromosomas.size(); cromosoma++){
				buffer_de_soluciones << generacion << "\t" << especie << "\t" <<  cromosomas[especie][cromosoma].to_string() << "\n";
			}
			
			//Guarda en archivo el promedio del fitness
			if(generacion % cada_cuantas_generaciones_mostrar_datos == 0){
				if(generacion == cada_cuantas_generaciones_mostrar_datos && especie == 0){
					std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
					auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
					double duracion_corrida = (duration/(cada_cuantas_generaciones_mostrar_datos*1000000.0L))*generaciones;
					cout << "Duracion estimada de corrida: " << duracion_corrida << "\n";
				}
				int idx_mejor_solucion = mejor_solucion(cromosomas[especie]);
				cout << "Especie|mejor fitness|clusters: " << especie << "|" << cromosomas[especie][idx_mejor_solucion].fitness << "|" << cromosomas[especie][idx_mejor_solucion].clusters << "\n";
				cout << "Promedio/Optima (paso): " << promedio_fitness(cromosomas[especie]) << "/" << cromosoma_objetivo.fitness << " (" << generacion << "/" << generaciones << ")\n";
				archivo_mejores_soluciones << cromosomas[especie][idx_mejor_solucion].to_string() << "\n";
				//archivo_mejores_soluciones_formato_cluster << cromosomas[idx_mejor_solucion].to_cluster() << "\n";
				if(cromosomas[especie][idx_mejor_solucion].fitness == max_fitness[especie]){
					iteraciones_con_max_fitness_igual[especie] = iteraciones_con_max_fitness_igual[especie] + cada_cuantas_generaciones_mostrar_datos;
				}else{
					max_fitness[especie] = cromosomas[especie][idx_mejor_solucion].fitness;
					iteraciones_con_max_fitness_igual[especie] = 1;
				}
				if(iteraciones_con_max_fitness_igual[especie] >= cutoff){
					cout << "Cutoff alcanzado...finalizando...\n";
					break;
				}	
			}
		}
		
		//Mezcla las mejores soluciones de cada especie
		if((generacion+1) % cada_cuantas_generaciones_intercambiar_especies == 0 && especies > 1){
			vector<vector<c_cromosoma>> intercambiar_cromosomas(especies);
			
			//Reserva las soluciones a intercambiar
			for(int especie = 0; especie < especies; especie++){
				reservar_soluciones_de_elite(cromosomas[especie], intercambiar_cromosomas[especie], soluciones_de_elite_a_intercambiar);
			}
			
			//Intercambia las soluciones
			for(int especie = 0; especie < (especies-1); especie++){
				aplicar_soluciones_de_elite(cromosomas[especie], intercambiar_cromosomas[especie+1], generador_aleatorio);
			}
			aplicar_soluciones_de_elite(cromosomas[(especies-1)], intercambiar_cromosomas[0], generador_aleatorio);
			cout << "Intercambio de especies (nichos)\n";
			
		}
		
		//Splitea los clusters mas grandes
		if((generacion+1) % cada_cuantas_generaciones_splitear_clusters == 0){
			vector<vector<c_cromosoma>> intercambiar_cromosomas(especies);
			cout << "Spliteando clusters\n";
			//Reserva las soluciones a intercambiar
			for(int especie = 0; especie < especies; especie++){
				for(int cromosoma = 0; cromosoma < poblacion; cromosoma++){
						cromosomas[especie][cromosoma].splitear_el_cluster_mas_grande();
				}
			}
			
		}		

	}
	
	for(int especie = 0; especie < especies; especie++){
		int idx_mejor_solucion = mejor_solucion(cromosomas[especie]);
		cout << "Especie|mejor fitness|clusters: " << especie << "|" << cromosomas[especie][idx_mejor_solucion].fitness << "|" << cromosomas[especie][idx_mejor_solucion].clusters << "\n";
		cout << "Promedio/Optima  (paso): " << promedio_fitness(cromosomas[especie]) << "/" << cromosoma_objetivo.fitness << " (" << generacion << "/" << generaciones << ")\n";	
	}
	
	std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
	cout << "Duracion real: " << duration/(1000000.0L) << "\n";

 	archivo_mejores_soluciones.close();
 	//archivo_mejores_soluciones_formato_cluster.close();
	
 	for(int p = 0; p < puntos.size(); p++){
		archivo_puntos << puntos[p]->x << "\t" << puntos[p]->y << "\n";
	}
	archivo_puntos.close();

	archivo_soluciones << buffer_de_soluciones;
	archivo_soluciones.close();
	
	//Liberamos la memoria
	for(vector<c_punto*>::const_iterator it = puntos.begin(); it != puntos.end(); it++)
	{
	    delete (*it);
	} 
	puntos.clear();
	centros.clear();
	cromosomas.clear();
	cromosomas_nuevos.clear();
	cromosomas_de_elite.clear();

	return 0;
}


