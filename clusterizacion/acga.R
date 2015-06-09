#TO DO
#VER COMO IMPLEMENTAR LA FUNCION SAMPLE EN C!

#///////////////
#CONFIGURACIONES
#///////////////

#---------------------------
#Configuraciones del dataset
#---------------------------
dim_red = 2; #los puntos en la red no son reales, son solo los lugares alrededor de los cuales se van a armar los clusters
puntos_por_cluster = 10;
parametro_de_red = 1;
ancho_del_cluster = 0.1; #lo que mide el cluster en x
alto_del_cluster = 0.1; #lo que mide el cluster en y

#---------------------
#Configuraciones de AG
#---------------------
poblacion = 10;
pm = 0.1; #probabilidad de mutacion
pc = 0.3; #probabilidad de single-point crossover
pp = 3; #Cromosomas a seleccionar aleatoreamente en cada busqueda de padres
k_max = 10; #Maxima cantidad de clusters a buscar
alfa = 0.1; #Amplitud de mutacion de valores de activacion
epsilon = 0.3; #Amplitud de mutacion de centroides
soluciones_de_elite = 4; #Las mejores soluciones pasan sin alteraciones a la proxima generacion
generaciones = 1;
corridas = 1;


#////////////////////
#CARGADO DE LIBRERIAS
#////////////////////
#library(fpc);

#///////////////////////
#COMIENZAN LAS FUNCIONES
#///////////////////////

#------------------
#Funcion de fitness
#------------------
calcular_fitness <- function(cromosoma, promedio, k_max, puntos){

	#Va a usar como funcion de fitness el indice Calinski-Harabasz 
	#Guarda la distancia de cada punto a su correspondiente cluster
	distancias = matrix(0, ncol=1, nrow=nrow(puntos));

	#Cuantos clusters hay habilitados
	clusters = 0;
	for(i in 1:k_max){
		#Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
		if(cromosoma[i] >= 0.5) clusters = clusters + 1;
	}
	
	#Si solo hay un cluster, devuelve 0
	if(clusters < 2) return (0);

	#Asigna primero cada punto a un cluster por proximidad
	for(punto in 1:nrow(puntos)){
		
		#Inicializa la distancia al cluster al que pertenece (inf = no pertenece a ninguno)		
		distancia_a_k = Inf;
		for(i in 1:k_max){
			#Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
			if(cromosoma[i] >= 0.5){
				#Si la distancia a este centroide es la mas chica, se lo asigno al centroide
				distancia_a_i = (puntos[punto, 1] - cromosoma[(k_max + (2*i - 1))])^2 + (puntos[punto, 2] - cromosoma[(k_max + (2*i))])^2
				if(distancia_a_i < distancia_a_k) {
					distancia_a_k = distancia_a_i;
				}
			}		
			distancias[punto] = distancia_a_k; #La distancia cuadriatica del punto al cluster que pertenece
		
		}
	}
	#Suma todas las distancias intra-clusters
	distancia_intra_cluster = sum(distancias);
		
	#Teniendo todo solo resta calcular el indice CH
	ch = (((promedio - distancia_intra_cluster)/distancia_intra_cluster)*(nrow(puntos)-clusters)/(clusters-1));
	return (ch);
}

#-------------------
#Funcion de mutacion
#-------------------
mutar <- function(cromosoma, pm, k_max, alfa, epsilon, valores_limite){
	
	#Muta primero los valores de activacion
	mascara = sample(0:1, k_max, replace=TRUE, c((1-pm),pm));
	for(i in 1:k_max){
		if(mascara[i]){
	
			#muta
			nuevo_va = cromosoma[i] + runif(min=-alfa, max=alfa, 1);
			if(nuevo_va > 1){			
				cromosoma[i] = 1;
			}else if(nuevo_va < 0){
				cromosoma[i] = 0;
			}else{
				cromosoma[i] = nuevo_va;
			}
		}
		
	}

	#Muta ahora los valores de los centroides
	mascara = sample(0:1, k_max, replace=TRUE, c((1-pm),pm));
	t = c(0,0);
	t[1] = epsilon*(valores_limite[1,2] - valores_limite[1,1]); #Amplitud de la mutacion en x
	t[2] = epsilon*(valores_limite[2,2] - valores_limite[2,1]); #Amplitud de la mutacion en y
	for(i in 1:k_max){
		if(mascara[i]){
	
			#muta
			nuevo_x = cromosoma[(k_max+(2*i - 1))] + runif(min=-t[1], max=t[1], 1);
			nuevo_y = cromosoma[(k_max+(2*i))] + runif(min=-t[2], max=t[2], 1);
			if(nuevo_x > valores_limite[1,2]){			
				cromosoma[(k_max+(2*i - 1))] = valores_limite[1,2];
			}else if(nuevo_x < valores_limite[1,1]){
				cromosoma[(k_max+(2*i - 1))] = valores_limite[1,1];
			}else{
				cromosoma[(k_max+(2*i - 1))] = nuevo_x;
			}

			if(nuevo_y > valores_limite[2,2]){			
				cromosoma[(k_max+(2*i))] = valores_limite[2,2];
			}else if(nuevo_y < valores_limite[2,1]){
				cromosoma[(k_max+(2*i))] = valores_limite[2,1];
			}else{
				cromosoma[(k_max+(2*i))] = nuevo_y;
			}

		}
		
	}
	return (cromosoma);
}

#--------------------
#Funcion de crossover
#--------------------
cruzar <- function(cromosomas_padres, pc, k_max){
	#Creamos los hijos
	cromosomas_hijos = cromosomas_padres; 

	#Hace crossover entre dos cromosomas con probabilidad pc para cada elemento
	mascara = sample(0:1, k_max, replace=TRUE, c((1-pc),pc));
	for(i in 1:k_max){
		if(mascara[i]){
	
			#intercambiamos los locus
			cromosomas_hijos[1, (k_max+(2*i - 1)):(k_max+(2*i))] = cromosomas_padres[2, (k_max+(2*i - 1)):(k_max+(2*i))];
			cromosomas_hijos[2, (k_max+(2*i - 1)):(k_max+(2*i))] = cromosomas_padres[1, (k_max+(2*i - 1)):(k_max+(2*i))];
		}
		
	}
	return (cromosomas_hijos);
}

#-----------------------------------------------------
#Funcion que elige una pareja en funcion de su fitness
#-----------------------------------------------------
elegir_pareja <- function(fitness, pp){

	#Toma pp soluciones aleatoreas y nos quedamos con las dos de mejor fitness	
	cromosomas <- sample(1:length(fitness), pp, replace=FALSE);
	pareja = c(0,0);
	pareja[1] = which.max(fitness[cromosomas]);
	fitness = fitness[-pareja[1]];
	pareja[2] = which.max(fitness[cromosomas]);	
	return (pareja);
	
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

#Trae los valores maximos y minimos para cada dimension de los puntos
valores_limite = matrix(0, ncol=2, nrow=2);
valores_limite[1,1] = min(puntos[, 1]); #[1,1] es el minimo en x
valores_limite[1,2] = max(puntos[, 1]); #[1,2] es el maximo en x
valores_limite[2,1] = min(puntos[, 2]); #[2,1] es el minimo en y
valores_limite[2,2] = max(puntos[, 2]); #[2,2] es el maximo en y

#Matriz en blanco que va a guardar los cromosomas de la poblacion nueva despues de cada corrida
#Cada cromosoma es una tira ordenada de k_max valores de activacion + 2k_max de la informacion de cada centroide
cromosomas_nuevos = matrix(0, ncol=(3*k_max), nrow=poblacion);

#Matriz que guarda el fitness de cada cromosoma
fitness = matrix(0, ncol=1, nrow=poblacion);

#Cantidad de cruzas por iteracion			
cruzas = c(1: as.integer(poblacion/2));

#Para el fitness de los cromosomas vamos a necesitar el promedio 
#de la distancia de los puntos al centroide. Calculamos esto una sola vez.
centroide = t(c(mean(puntos[, 1]), mean(puntos[, 2])));
promedio = sum((puntos[, 1]-centroide[1])^2+(puntos[, 2]-centroide[2])^2);

#Comienzan las corridas
for(corrida in 1:corridas){

	#Genera los cromosomas al azar de la corrida. Primero los valores de activacion
	cromosomas = matrix(runif(k_max*poblacion), ncol=k_max);
	#Revisa que haya al menos dos con más de 0.5
	for(i in 1:poblacion){
		
		#Si no hay dos con mas de 0.5 elige dos posiciones al azar y le pone dos numeros mayores a 0.5
		if(length(which(cromosomas[i, ] > 0.5) < 2)){
			va_a_cambiar = sample(1:k_max, 2, replace=FALSE);
			cromosomas[i, va_a_cambiar[1]] = runif(1, min=0.5);
			cromosomas[i, va_a_cambiar[2]] = runif(1, min=0.5);
		}
	}
	#Ahora genera los centroides. Los genera dentro de los rangos maximos y minimos de los valores limite
	cromosomas = cbind(cromosomas, matrix(runif(k_max*poblacion, min=valores_limite[1,1], max=valores_limite[1,2]), ncol=k_max));
	cromosomas = cbind(cromosomas, matrix(runif(k_max*poblacion, min=valores_limite[2,1], max=valores_limite[2,2]), ncol=k_max));

	#Generacionando
	for(generacion in 1:generaciones){

		#Calcula el fitness de los cromosomas
		for(cromosoma in 1:poblacion){
			fitness[cromosoma] = calcular_fitness(cromosomas[cromosoma, ], promedio, k_max, puntos);
		}
		
		if(generacion%%1 == 0){
                        ibestf<-which.max(fitness);
                        nn <- apply(cromosomas,1,function(x){ return(length(unique(x)))})
			cat(paste("generacion:",generacion,"- fitness mean:sd:max:optimo", 
				  signif(mean(fitness),2),
				  signif(sd(fitness),2),
			          signif(fitness[ibestf],2),
				  signif(fitness_objetivo, 2),
				  "\n"))
			cat(paste("               - N mean:sd:max",
			          mean(nn),sd(nn),nn[ibestf],"\n\n"))
		}

		#Las soluciones con mejor fitness pasan inalteradas
		indice_mejores_soluciones = sort(fitness, index.return=TRUE)$ix[length(fitness):(length(fitness)-soluciones_de_elite + 1)];
		mejores_soluciones = cromosomas[indice_mejores_soluciones, ];

		#Cruza los cromosomas de acuerdo a su fitness. Cuanto mas fitness mas probabilidad de cruza. 
		#Elige poblacion/2 parejas

		pareja_actual = 1; #Indice de la nueva pareja en cada cruza, es interno a este bucle
		for(cruza in cruzas){
			#Elige la pareja a cruzar
			pareja = elegir_pareja(fitness, pp);
			#La cruza y genera dos hijos
			hijos = cruzar(cromosomas[pareja, ], pc, k_max);
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
			cromosomas[cromosoma, ] = mutar(cromosomas[cromosoma, ], pm, k_max, alfa, epsilon, valores_limite);
		
		}
		
		#Descartamos los cambios a las soluciones de elite
		cromosomas[indice_mejores_soluciones, ] = mejores_soluciones;		
	
	}
	
}

#Muestra los mejores fitness
graphics.off();
plot(puntos[,1],puntos[,2]);
mejor_cromosoma = cromosomas[which.max(fitness), ];
ks_habilitados = which(mejor_cromosoma[1:8] >= 0.5);
print(length(ks_habilitados));
points(mejor_cromosoma[(k_max+(ks_habilitados*2)-1)],mejor_cromosoma[(k_max+(ks_habilitados*2))],col=rainbow(k_max),pch=20);


