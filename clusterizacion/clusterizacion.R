#TO DO
#VER COMO IMPLEMENTAR LA FUNCION SAMPLE EN C!

#///////////////
#CONFIGURACIONES
#///////////////

#---------------------------
#Configuraciones del dataset
#---------------------------
poblacion = 50;
pm = 0.5; #probabilidad de mutacion
pc = 0.5; #probabilidad de single-point crossover
generaciones = 2500;
corridas = 1;

#---------------------
#Configuraciones de AG
#---------------------
dim_red = 2; #los puntos en la red no son reales, son solo los lugares alrededor de los cuales se van a armar los clusters
puntos_por_cluster = 10;
parametro_de_red = 1;
ancho_del_cluster = 0.1; #lo que mide el cluster en x
alto_del_cluster = 0.1; #lo que mide el cluster en y
k_max = 25; #Maxima cantidad de clusters a buscar
k_min = 2; #Minima cantidad de clusters a buscar
soluciones_de_elite = 4; #Las mejores soluciones pasan sin alteraciones a la proxima generacion

#Setea la semilla aleatoria para tener resultados reproducibles
set.seed(123457)

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
calcular_fitness <- function(puntos, cromosoma){
	#Vamos a probar con el indice de Calinski-Harabasz (ch)
	return (calinhara(puntos, cromosoma));
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm, k_max, k_min){
	#Tres operadores de mutacion: mutar, mergear, splitear
	
	#Muta el cromosoma con probabilidad pm
	if(runif(1) <= pm){

		#Elije un locus al azar y lo cambia
		posicion = sample(1:length(cromosoma), 1);
		red<-c(1:k_max);
		red<-red[-cromosoma[posicion]];	#Sacamos como posibilidad que la mutacion lo deje en el mismo lugar	
		cluster_anterior_a_la_mutacion = cromosoma[posicion]; #Guardamos la mutacion anterior
		cromosoma[posicion] = sample(red, 1);
		#Si la mutacion provoco que la solucion tenga menos clusters que el minimo, volvemos el cambio para atras
		if(length(unique(cromosoma)) < k_min) cromosoma[posicion] = cluster_anterior_a_la_mutacion;
		
	
	}
	#Mergea dos clusters dentro del cromosomas con probabilidad pm
	if(runif(1) <= pm){
		#Cuantos clusters hay
		clusters_en_cromosoma = unique(cromosoma);
		
		#Como minimo tiene que tener k_min clusters
		if(length(clusters_en_cromosoma) > k_min){
			#Elije dos clusters al azar y los junta
			clusters_a_mergear = sample(clusters_en_cromosoma, 2, replace=FALSE);
			cromosoma[which(cromosoma==clusters_a_mergear[2])]=clusters_a_mergear[1];	
		}
	}
	#Splitea un cluster dentro del cromosomas con probabilidad pm
	if(runif(1) <= pm){
		#Elije un cluster al azar
		clusters_en_cromosoma = unique(cromosoma);
		
		#Elije un cluster al azar y lo divide (elige en realidad dos clusters al azar,
		#uno el que va a dividir, el otro el que va a usar para asignar 
		#a la mitad de los elementos del primero).
		clusters_a_dividir = sample(clusters_en_cromosoma, 2, replace=FALSE);
		elementos_del_cluster_a_dividir = which(cromosoma==clusters_a_dividir[1]);
		
		cromosoma[elementos_del_cluster_a_dividir[1:ceiling(length(elementos_del_cluster_a_dividir)/2)]]=clusters_a_dividir[2];	
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
	return (sample(1:length(fitness), 2, replace=FALSE, prob=(fitness/sum(fitness))));
}

#---------------------------------------
#Funcion para generar el dataset inicial
#---------------------------------------
generar_dataset <- function(dim_red, puntos_por_cluster, parametro_de_red, ancho_del_clustero, alto_del_cluster){
	#genero la red equiespaciada
	a <- seq(1, dim_red*parametro_de_red, parametro_de_red);
	red <- matrix(data=a, nrow=dim_red^2, ncol=2);
	red[, 1] <- rep(a, each=dim_red);

	#genero los puntos de datos alrededor de la red
	puntos_en_la_red <- dim_red^2;
	total_de_puntos <- puntos_en_la_red * puntos_por_cluster;

	#Genero los puntos de los clusters
	puntos <- matrix(0, nrow=total_de_puntos, ncol=2);
	puntos[, 1] <- runif(total_de_puntos, -ancho_del_cluster, ancho_del_cluster) + rep(red[, 1], each=puntos_por_cluster);
	puntos[, 2] <- runif(total_de_puntos, -alto_del_cluster, alto_del_cluster) + rep(red[, 2], each=puntos_por_cluster);
	return (puntos);
}

#////////////////////
#COMIENZA EL PROGRAMA
#////////////////////

#Genera el dataset
puntos <- generar_dataset(dim_red, puntos_por_cluster, parametro_de_red, ancho_del_clustero, alto_del_cluster);
total_de_puntos = nrow(puntos);

#Matriz de distancias entre los puntos
matriz_de_disimilaridad = dist(puntos);

#Matriz en blanco que va a guardar los cromosomas de la poblacion nueva despues de cada corrida
#Cada cromosoma es una tira ordenada que asigna a cada posicion (cada punto) uno de los clusters posibles
#de la red
cromosomas_nuevos = matrix(0, ncol=total_de_puntos, nrow=poblacion);

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

#Cantidad de cruzas por iteracion			
cruzas = c(1: as.integer(poblacion/2));

#Registro de fitness
registro_de_fitness = matrix(0, ncol=1, nrow=generaciones);

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas al azar de la corrida, entre 1 y la cantidad de puntos de la red
	cromosomas = matrix(sample(1:k_max, poblacion*total_de_puntos, replace=TRUE), ncol=total_de_puntos);

	#Generando las generaciones
	for(generacion in 1:generaciones){

		#Calcula el fitness de los cromosomas
		for(cromosoma in 1:poblacion){
			fitness[cromosoma] = calcular_fitness(puntos, cromosomas[cromosoma, ]);
		}
		registro_de_fitness[generacion] = mean(fitness);
		
		#Las soluciones con mejor fitness pasan inalteradas
		indice_mejores_soluciones = sort(fitness, index.return=TRUE)$ix[length(fitness):(length(fitness)-soluciones_de_elite + 1)];
		mejores_soluciones = cromosomas[indice_mejores_soluciones, ];		
		
		if(generacion%%1 == 0){
                        ibestf<-which.max(fitness)
                        nn <- apply(cromosomas,1,function(x){ return(length(unique(x)))})
			cat(paste("generacion:",generacion,"- fitness mean:sd:max", 
				  signif(mean(fitness),2),signif(sd(fitness),2),
			          signif(fitness[ibestf],2),"\n"))
			cat(paste("               - N mean:sd:max",
			          mean(nn),sd(nn),nn[ibestf],"\n\n"))
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
			cromosomas[cromosoma, ] = mutar(cromosomas[cromosoma, ], pm, k_max, k_min);
		
		}
		
		#Descartamos los cambios a las soluciones de elite
		cromosomas[indice_mejores_soluciones, ] = mejores_soluciones;			
		
	}
	
}

#Muestra los mejores fitness
graphics.off();
ibestf<-which.max(fitness);
plot(puntos[,1],puntos[,2]);
points(puntos[,1],puntos[,2],col=rainbow(length(unique(cromosomas[ibestf, ])))[cromosomas[ibestf, ]],pch=20);
dev.new();
plot(silhouette(cromosomas[ibestf, ], matriz_de_disimilaridad));
dev.new()
plot(registro_de_fitness);
