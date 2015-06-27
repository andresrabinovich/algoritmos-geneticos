#include <iostream>
#include <random>
#include <functional> //Para std::bind
#include <vector>

using namespace std;

/**************
CONFIGURACIONES
**************/

/*-------------------------
Configuraciones del dataset
-------------------------*/
int centros_por_dimension = 2; //los centros son los lugares alrededor de los cuales se van a armar los clusters
int puntos_por_cluster = 10;
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
int generaciones = 2500;

/*--------------------------------
Calcula algunas cantidades previas
--------------------------------*/
int cantidad_de_centros = centros_por_dimension*centros_por_dimension;
int cantidad_de_puntos  = cantidad_de_centros*puntos_por_cluster;

/*-------------------------------------------
Inicializa el generador de numeros aleatorios
-------------------------------------------*/
default_random_engine generador_aleatoreo;

/************************
COMIENZAN LAS ESTRUCTURAS
************************/
struct s_punto{
	float x;
	float y;
};

struct s_rectangulo {
	s_punto inferior_izquierdo;
	s_punto superior_derecho;
};

/**********************
COMIENZAN LAS FUNCIONES
**********************/

/*-------------------------------------
Funcion para generar el dataset inicial
-------------------------------------*/
void generar_dataset(vector<s_punto> &centros, vector<s_punto> &puntos, valores_limite){

	//Habilita el generador de numeros aleatorios para elegir el ancho y el alto
	uniform_real_distribution<double> urd_ancho(-ancho_del_cluster, ancho_del_cluster);
	uniform_real_distribution<double> urd_alto(-alto_del_cluster, alto_del_cluster);

	//Para usos repetidos
	auto runif_ancho = bind(urd_ancho, generador_aleatoreo);
	auto runif_alto = bind(urd_alto, generador_aleatoreo);

	//Resetea los valores limite
	valores_limite.inferior_izquierdo.x = 0;
	valores_limite.inferior_izquierdo.y = 0;
	valores_limite.superior_derecho.x = 0;
	valores_limite.superior_derecho.y = 0; 

	//genera la red equiespaciada
	int centro = 0;	
	for(int x = 1; x <= centros_por_dimension; x++){
		for(int y = 1; y <= centros_por_dimension; y++){
			centros[centro].x = x;
			centros[centro].y = y;
			centro++;
		}
	}	

	//genera los puntos de datos alrededor de los centros
	int punto = 0;
	for(centro = 0; centro < cantidad_de_centros; centro++){
		for(int punto_en_cluster = 0; punto_en_cluster < puntos_por_cluster; punto_en_cluster++){
			puntos[punto].x = centros[centro].x + runif_ancho();
			puntos[punto].y = centros[centro].y + runif_alto();
			punto++;
		}		
	}	
	return;
}

/*******************
COMIENZA EL PROGRAMA
*******************/
int main(){
	vector<s_punto> centros(cantidad_de_centros);
	vector<s_punto> puntos(cantidad_de_puntos);
	s_rectangulo valores_limite;//Va a contener los valores limite que tomaron los puntos, en forma de un rectangulo (o sea, dos puntos, el angulo inferior izquierdo y el superior derecho)
	generar_dataset(centros, puntos, valores_limite);
	for(int centro = 0; centro < cantidad_de_centros; centro++){
		cout << "x: " << centros[centro].x << " y: " << centros[centro].y << "\n";
	}
	for(int punto = 0; punto < cantidad_de_puntos; punto++){
		cout << "x: " << puntos[punto].x << " y: " << puntos[punto].y << "\n";
	}
	return 0;
}


