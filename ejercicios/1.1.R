#---------------
#Configuraciones
#---------------
poblacion = 100;
tamano_del_cromosoma = 20;
pm = 0.001; #probabilidad de mutacion
pc = 0.1; #probabilidad de single-point crossover
generaciones = 225;
corridas = 20;

#Matriz en blanco que va a guardar los cromosomas de la poblacion nueva despues de cada corrida
cromosomas_nuevos = matrix(0, ncol=tamano_del_cromosoma, nrow=poblacion);

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

aparicion_del_1 = matrix(0, ncol=1, nrow=corridas);

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas binarios al azar de la corrida
	cromosomas = matrix(round(runif(poblacion*tamano_del_cromosoma)), ncol=tamano_del_cromosoma);

	#Generando las generaciones
	for(generacion in 1:generaciones){

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

		#Nos fijamos si aparecio el string de 1s
		if(any(fitness==20)){
			aparicion_del_1[corrida] = generacion;
			break;	
		}
	}
}
#Promedio de la generacion en la que aparece el string de 1s
print(mean(aparicion_del_1));

#------------------
#Funcion de fitness
#------------------
calcular_fitness <- function(cromosoma){
	return (sum(cromosoma));	
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm){
	#Muta un bit del cromosoma con probabilidad pm
	if(runif(1) <= pm){
		#Elije un locus al azar y lo invierte
		posicion = round(runif(1, 1, length(cromosoma)));
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
		posicion = round(runif(1, 1, length(cromosomas_padres[1,])));
		cromosomas_hijos[1, 1:posicion] = cromosomas_padres[2, 1:posicion];
		cromosomas_hijos[2, 1:posicion] = cromosomas_padres[1, 1:posicion];
	}
	return (cromosomas_hijos);
}

#-----------------------------------------------------
#Funcion que elige una pareja en funcion de su fitness
#-----------------------------------------------------
elegir_pareja <- function(fitness){

	#Trae los cromosomas padre usando el metodo de la ruleta. 
	#La probabilidad de ser elegido es el fitness_individual/fitness_total
	#Es un poco de bricollage
	pareja = c(1, 1);
	poblacion = length(fitness);
	for(padre in 1:2){	
		#No queremos que salga el mismo padre dos veces asi que lo sacamos haciendo su fitness = 0
		if(padre == 2) fitness[pareja[1]] = 0; 
		#Traemos un numero entre 0 y el total del fitness y lo macheamos con el cromosoma correspondiente
		ruleta = round(runif(1, 1, sum(fitness)));
		#Encontramos el cromosoma correspondiente
		suma = 0;
		for(cromosoma in 1:poblacion){
			suma = suma + fitness[cromosoma];
			if(suma >= ruleta) break;
			pareja[padre] = pareja[padre] + 1;
		}
	}
	return (pareja);
}

