#///////////////
#CONFIGURACIONES
#///////////////
poblacion = 1;
tamano_del_cromosoma = 20;
pm = 0.001; #probabilidad de mutacion
pc = 0.1; #probabilidad de single-point crossover
generaciones = 10000;
corridas = 1;

#///////////////////////
#COMIENZAN LAS FUNCIONES
#///////////////////////

#------------------
#Funcion de fitness
#------------------
calcular_fitness <- function(cromosoma){
	#return (binario_a_decimal(cromosoma));	
	return (sum(cromosoma));
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm){
	#Elije un locus al azar y lo invierte
	posicion = sample(1:length(cromosoma), 1);
	if(cromosoma[posicion]){
		cromosoma[posicion] = 0;
	}else{
		cromosoma[posicion] = 1;
	}
	return (cromosoma);
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

#////////////////////
#COMIENZA EL PROGRAMA
#////////////////////

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

#Mejores fitness
mejores_fitness = matrix(0, ncol=1, nrow=1);

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas binarios al azar de la corrida
	cromosomas = matrix(sample(0:1, poblacion*tamano_del_cromosoma, replace=TRUE), ncol=tamano_del_cromosoma);

	#Calcula el fitness de los cromosomas
	for(cromosoma in 1:poblacion){
		fitness[cromosoma] = calcular_fitness(cromosomas[cromosoma, ]);
	}	
	#Mejor fitness de la primera corrida
	mejores_fitness[1] = fitness[1];

	#Generando las generaciones
	for(generacion in 1:generaciones){
		
		print(mean(fitness));

		#Mutamos los nuevos cromosomas y los guardamos si su fitness es superior al anterior
		for(cromosoma in 1:poblacion){
			cromosoma_nuevo = mutar(cromosomas[cromosoma, ], pm);
			fitness_nuevo = calcular_fitness(cromosoma_nuevo);
			if(fitness_nuevo > fitness[cromosoma]){
				cromosomas[cromosoma, ] = cromosoma_nuevo;
				fitness[cromosoma] = fitness_nuevo;
			}
		}

		if(generacion%%100 == 0){
			mejores_fitness = rbind(mejores_fitness, fitness[1]);
		}
	}
}
j = 1;
mejores_fitness_promediados = matrix(0, ncol=length(mejores_fitness));
for(i in seq(1, length(mejores_fitness), by=10)){
	mejores_fitness_promediados[j] = mean(mejores_fitness[i:(i+9)]);
	j = j + 1;
}
plot(mejores_fitness_promediados);
