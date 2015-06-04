#///////////////
#CONFIGURACIONES
#///////////////
poblacion = 50;
pm = 0.01; #probabilidad de mutacion
pc = 0.8; #probabilidad de single-point crossover
generaciones = 250;
corridas = 1;

dim_red = 2; #los puntos en la red no son reales, son solo los lugares alrededor de los cuales se van a armar los clusters
puntos_por_cluster = 10;
parametro_de_red = 1;
ancho_del_cluster = 0.1; #lo que mide el cluster en x
alto_del_cluster = 0.1; #lo que mide el cluster en y
k = dim_red^2; #TO DO - Ver de adaptar el algoritmo para no tener que decirle cuantos clusters buscar

#////////////////////
#CARGADO DE LIBRERIAS
#////////////////////
library(fpc);

#/////////////////////////////////
#COMIENZAN LAS FUNCIONES GENERICAS
#/////////////////////////////////

#------------------
#Funcion de fitness
#------------------
calcular_fitness <- function(cromosoma, matriz_de_distancias){
	#return (binario_a_decimal(cromosoma));	
	#return (sum(cromosoma));
	#Vamos a probar con silhouette como funcion de fitness
	#s<-silhouette(cromosoma, matriz_de_distancias);
	#si<-summary(s);
	#return (si$si.summary["Mean"]+1); #Le sumamos 1 para que sea positiva
	#Vamos a probar con el indice de Calinski-Harabasz (ch)
	return (cluster.stats(matriz_de_distancias, cromosoma)$ch);
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm, k){
	#Muta un bit del cromosoma con probabilidad pm
	if(runif(1) <= pm){
		#Elije un locus al azar y lo cambia
		posicion = sample(1:length(cromosoma), 1);
		red<-c(1:k);
		red<-red[-cromosoma[posicion]];	#Sacamos como posibilidad que la mutacion lo deje en el mismo lugar	
		cromosoma[posicion] = sample(red, 1);
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
	
		#Elije un locus desde donde empezar a cruzar y los cruza
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

#------------------------------------
#Genera la matriz de distancias entre 
#los puntos y cada elemento de la red
#------------------------------------
#generar_matriz_de_distancias <- function(puntos, redes){
#	matriz_distancia = matrix(0, nrow=nrow(puntos), ncol=nrow(redes));
#	for(punto in 1:nrow(puntos)){
#		for(red in 1:nrow(redes)){
#			matriz_distancia[punto, red] = sqrt((puntos[punto, 1]-redes[red, 1])^2 + (puntos[punto, 2]-redes[red, 2])^2);
#		}
#	}
#	return (matriz_distancia);
#}

#////////////////////
#COMIENZA EL PROGRAMA
#////////////////////

#genero la red equiespaciada
a <- seq(1, dim_red*parametro_de_red, parametro_de_red);
red <- matrix(data=a, nrow=dim_red^2, ncol=2);
red[, 1] <- rep(a, each=dim_red);

#genero los puntos de datos alrededor de la red
puntos_en_la_red <- dim_red^2;
total_de_puntos <- puntos_en_la_red * puntos_por_cluster;
tamano_del_cromosoma = total_de_puntos;

#Los puntos de ruido que vamos a tener
cantidad_de_puntos_de_ruido <- total_de_puntos * nivel_de_ruido;

#Genero los puntos de los clusters
puntos <- matrix(0, nrow=total_de_puntos, ncol=2);
puntos[, 1] <- runif(total_de_puntos, -ancho_del_cluster, ancho_del_cluster) + rep(red[, 1], each=puntos_por_cluster);
puntos[, 2] <- runif(total_de_puntos, -alto_del_cluster, alto_del_cluster) + rep(red[, 2], each=puntos_por_cluster);

#Matriz de distancias entre los puntos y la red
#matriz_de_distancias = generar_matriz_de_distancias(puntos, red);
matriz_de_disimilaridad = dist(puntos);

#Matriz en blanco que va a guardar los cromosomas de la poblacion nueva despues de cada corrida
#Cada cromosoma es una tira ordenada que asigna a cada posicion (cada punto) uno de los clusters posibles
#de la red
cromosomas_nuevos = matrix(0, ncol=total_de_puntos, nrow=poblacion);

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

#Cantidad de cruzas por iteracion			
cruzas = c(1: as.integer(poblacion/2));

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas al azar de la corrida, entre 1 y la cantidad de puntos de la red
	cromosomas = matrix(sample(1:puntos_en_la_red, poblacion*tamano_del_cromosoma, replace=TRUE), ncol=tamano_del_cromosoma);

	#Generando las generaciones
	for(generacion in 1:generaciones){

		#Calcula el fitness de los cromosomas
		for(cromosoma in 1:poblacion){
			fitness[cromosoma] = calcular_fitness(cromosomas[cromosoma, ], matriz_de_disimilaridad);
		}
		
		if(generacion%%1 == 0){
			print(mean(fitness));
			print(generacion);
		}

		#Cruza los cromosomas de acuerdo a su fitness. Cuanto mas fitness mas probabilidad de cruza. 
		#Elige poblacion/2 parejas

		pareja_actual = 1; #Indice de la nueva pareja en cada cruza, es interno a este bucle
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
			cromosoma_nuevo = mutar(cromosomas[cromosoma, ], pm, k);
			#fitness_nuevo = calcular_fitness(cromosoma_nuevo, matriz_de_disimilaridad);
			#if(fitness_nuevo > fitness[cromosoma] || runif(1) > exp(-(fitness_nuevo-fitness[cromosoma]))) {
			#if(fitness_nuevo > fitness[cromosoma]) {
				cromosomas[cromosoma, ] = cromosoma_nuevo;
			#	fitness[cromosoma] = fitness_nuevo;
			#}
		
		}
		
	}
	
}

#Muestra los mejores fitness
graphics.off();
print(cromosomas);
plot(puntos[,1],puntos[,2]);
points(puntos[,1],puntos[,2],col=rainbow(k)[cromosomas[1, ]],pch=20);
dev.new();
plot(silhouette(cromosomas[which.max(fitness), ], matriz_de_disimilaridad));
