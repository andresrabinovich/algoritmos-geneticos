k_max = 20;
puntos <- read.delim("~/algoritmos-geneticos/clusterizacion/puntos.txt", header=FALSE)
soluciones <- data.matrix(read.delim("~/algoritmos-geneticos/clusterizacion/soluciones.txt", sep=",", header=FALSE))
plot(puntos)
for(i in nrow(soluciones):nrow(soluciones)){
   a = which(soluciones[i, 1:k_max] >= 0.5);
   b = matrix(0, ncol=2, nrow=length(a))
   k = 1;
   for(j in a){
     b[k, 1] = soluciones[i, (k_max+(2*j-1))]
     b[k, 2] = soluciones[i, (k_max+(2*j))]
     k = k + 1;
   }
   points(b, col="red")
}
cromosoma = soluciones[nrow(soluciones), ]

#Incializa la tira que representa la solucion en el formato de clusters
tira = rep(0, nrow(puntos));

#Si el valor de activacion es mayor a 0.5 (i.e. esta activado) busca los puntos que pertenecen a ese cluster		
k_activos = which(cromosoma[1:k_max] >= 0.5);
#Asigna cada punto a un cluster por proximidad
for(punto in 1:nrow(puntos)){
  
  #Inicializa la distancia al cluster al que pertenece (inf = no pertenece a ninguno)		
  distancia_a_k = Inf;
  k = 0; #El cluster al que va a pertenecer
  for(i in k_activos){
    #Si la distancia a este centroide es la mas chica, se lo asigno al centroide
    distancia_a_i = (puntos[punto, 1] - cromosoma[(k_max + (2*i - 1))])^2 + (puntos[punto, 2] - cromosoma[(k_max + (2*i))])^2
    if(distancia_a_i < distancia_a_k) {
      distancia_a_k = distancia_a_i;
      k = i;
    }
  }
  #Asigna el punto al cluster mas cercano
  tira[punto] = k; 
}
#Ordena la lista para que no queden agujeros entre los clusters (el 111133335555 tiene que ser 111122223333)
k_actual = 1
for(i in k_activos){
  elementos_a_modificar = which(tira == i);
  if(length(elementos_a_modificar)){
    tira[elementos_a_modificar] = k_actual;
    k_actual = k_actual + 1;
  }
}
dev.new();
plot(puntos)
points(puntos[,1],puntos[,2],col=rainbow(length(unique(tira)))[tira],pch=20);
