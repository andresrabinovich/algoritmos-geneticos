#///////////////
#CONFIGURACIONES
#///////////////
poblacion = 100;
tamano_del_cromosoma = 20;
pm = 0.001; #probabilidad de mutacion
pc = 0.1; #probabilidad de single-point crossover
generaciones = 100;
corridas = 1;
cantidad_de_esquemas = 10;

#///////////////////////
#COMIENZAN LAS FUNCIONES
#///////////////////////

#------------------
#Funcion de fitness
#------------------
calcular_fitness <- function(cromosoma){
	return (binario_a_decimal(cromosoma));	
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm){
	#Muta un bit del cromosoma con probabilidad pm
	if(runif(1) <= pm){
		#Elije un locus al azar y lo invierte
		posicion = sample(1:length(cromosoma), 1);
		if(cromosoma[posicion]){
			cromosoma[posicion] = 0;
		}else{
			cromosoma[posicion] = 1;
		}

	}
	return (cromosoma);
}

#--------------------
#Funcion de crossover
#--------------------
cruzar <- function(cromosomas_padres, pc){
	#Creamos los hijos
	cromosomas_hijos = cromosomas_padres; 

	#Hace crossover entre dos cromosomas con probabilidad pc
	if(runif(1) <= pc){
		#Elije un locus desde donde empezar a cruzar
		posicion = sample(1:length(cromosomas_padres[1,]), 1);
		cromosomas_hijos[1, 1:posicion] = cromosomas_padres[2, 1:posicion];
		cromosomas_hijos[2, 1:posicion] = cromosomas_padres[1, 1:posicion];
	}
	return (cromosomas_hijos);
}

#-----------------------------------------------------
#Funcion que elige una pareja en funcion de su fitness
#-----------------------------------------------------
elegir_pareja <- function(fitness){

	#VER COMO IMPLEMENTAR LA FUNCION SAMPLE EN C!
	return (sample(1:length(fitness), 2, replace=FALSE, prob=fitness));
}

#---------------------------------------
#Funcion para pasar de binario a decimal
#---------------------------------------
binario_a_decimal <- function(binario) {
	decimal = 0;
	longitud = length(binario);
	for(i in 1:longitud){
		decimal = decimal + (binario[i] * 2^(longitud - i));
	}
	return(decimal);
}

#--------------------------------------------
#Funcion para machear esquemas con cromosomas
#--------------------------------------------
comparar_cromosomas_con_esquemas <- function(cromosomas, esquemas){
	macheos = matrix(0, ncol=nrow(esquemas));

	#Vamos a chequear todos los cromosomas
	for(cromosoma in 1:nrow(cromosomas)){

		#Contra todos los esquemas
		for(esquema in 1:nrow(esquemas)){

			#Posicion por posicion
			machea = TRUE;#Bandera de comparacion	
			for(posicion in 1:length(esquemas[esquema])){

				if(esquemas[esquema][posicion] != 2){ #2 es el comodin

					if(esquemas[esquema][posicion] != cromosomas[cromosoma][posicion]){
						machea = FALSE;
						break;	
					}

				}

			}
			if(machea) macheos[esquema] = macheos[esquema] + 1;

		}

	}

	return (macheos);
}

#////////////////////
#COMIENZA EL PROGRAMA
#////////////////////

#Genera esquemas al azar
esquemas = matrix(sample(0:2, cantidad_de_esquemas*tamano_del_cromosoma, replace=TRUE), ncol=tamano_del_cromosoma);

#Indica cuantos esquemas se machearon en cada generacion
esquemas_por_generacion = matrix(0, ncol=cantidad_de_esquemas, nrow=generaciones)

#Matriz en blanco que va a guardar los cromosomas de la poblacion nueva despues de cada corrida
cromosomas_nuevos = matrix(0, ncol=tamano_del_cromosoma, nrow=poblacion);

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

#Guarda el mejor fitness (col 1) y el promedio (col 2) de cada generacion
fitness_por_generacion = matrix(0, ncol=2, nrow=generaciones);

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas binarios al azar de la corrida
	cromosomas = matrix(sample(0:1, poblacion*tamano_del_cromosoma, replace=TRUE), ncol=tamano_del_cromosoma);

	#Generando las generaciones
	for(generacion in 1:generaciones){

		esquemas_por_generacion[generacion, ] = comparar_cromosomas_con_esquemas(cromosomas, esquemas);

		#Calcula el fitness de los cromosomas
		for(cromosoma in 1:poblacion){
			fitness[cromosoma] = calcular_fitness(cromosomas[cromosoma, ]);
		}	
		print(mean(fitness));

		#Cruza los cromosomas de acuerdo a su fitness. Cuanto mas fitness mas probabilidad de cruza. 
		#Elige poblacion/2 parejas
		pareja_actual = 1; #Indice de la nueva pareja en cada cruza
		cruzas = c(1: as.integer(poblacion/2));
		for(cruza in cruzas){
			#Elige la pareja a cruzar
			pareja = elegir_pareja(fitness);
			#La cruza y genera dos hijos
			hijos = cruzar(cromosomas[pareja, ], pc);
			#Asigna a la nueva poblacion los dos hijos
			cromosomas_nuevos[pareja_actual, ] = hijos[1, ];
			cromosomas_nuevos[pareja_actual+1, ] = hijos[2, ];
			#Agrega dos al indice de nuevas parejas
			pareja_actual = pareja_actual + 2;
		}
	
		#Asignamos la nueva poblacion como la poblacion actual
		cromosomas = cromosomas_nuevos;

		#Mutamos los nuevos cromosomas
		for(cromosoma in 1:poblacion){
			cromosomas[cromosoma, ] = mutar(cromosomas[cromosoma, ], pm);
		}

		#Guarda el mejor fitness y el promedio
		fitness_por_generacion[generacion, 1] = max(fitness);
		fitness_por_generacion[generacion, 2] = mean(fitness);
		
	}
}

#Muestra los mejores fitness
graphics.off();
plot(fitness_por_generacion[, 2]);
dev.new();
plot(fitness_por_generacion[, 1]);
print(esquemas_por_generacion);
