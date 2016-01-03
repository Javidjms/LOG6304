##########################################################################
#Projet
##########################################################################

#Chargement de la librairie Rcurl
library(RCurl)

#Chargement des données
u.data <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.data.csv', userpwd='20113:20113'), sep='|', header=T)
u.user <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.user.csv', userpwd='20113:20113'), sep='|', header=T)
u.item <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.item.csv', userpwd='20113:20113'), sep='|', header=T)

#Jointure du tableau u.data avec u.user sur la colonne user.id
v.u <- merge(u.data, u.user, by.x='user.id', by.y='id')

#Chargement de la librairie Matrix
library(Matrix)

#Transformation en matrice (pour certaines opérations) avec les entetes des lignes u et colonnes i
m.sparse <- sparseMatrix(u.data[,1],u.data[,2],x=u.data[,3])
rownames(m.sparse) <- paste('u', 1:nrow(m.sparse), sep='')
colnames(m.sparse) <- paste('i', 1:ncol(m.sparse), sep='')
m <- as.matrix(m.sparse)

lambda=0.0000005

#Algorithme SGD with RMSE Criteria
algo.sgd = function(R,K,rmse.threshold = 0.1,gamma=0.0002,lambda=0.0000005) {
  H = dim(R)[1]
  W = dim(R)[2]
  P = matrix(rnorm(1:(H*K),mean = 0, sd = 0.01),H,K)
  Q = matrix(rnorm(1:(W*K),mean = 0, sd = 0.01),W,K)
  step = 1
  r = rmse(R,(P %*% t(Q)))
  while (r >= rmse.threshold){
    
    E = R - (P %*% t(Q))
    P = P + gamma * (E %*% Q - lambda * P)
    Q = Q + gamma * (t(E) %*% P - lambda * Q)
    print(step)
    r = rmse(R,(P %*% t(Q)))
    print(r)
    step = step +1
  }
  
  return(list(P,Q))
}


#Algorithme SGD with Regularized RMSE Criteria
algo.sgd.reg = function(R,K,rmse.reg.threshold = 0.1,gamma=0.0002,lambda=0.0000005) {
  H = dim(R)[1]
  W = dim(R)[2]
  P = matrix(rnorm(1:(H*K),mean = 0, sd = 0.01),H,K)
  Q = matrix(rnorm(1:(W*K),mean = 0, sd = 0.01),W,K)
  step = 1
  r = rmse.reg(R,P,Q,lambda)
  while (r >= rmse.reg.threshold){
    
    E = R - (P %*% t(Q))
    P = P + gamma * (E %*% Q - lambda * P)
    Q = Q + gamma * (t(E) %*% P - lambda * Q)
    print(step)
    r = rmse.reg(R,P,Q,lambda)
    print(r)
    step = step +1
  }
  
  return(list(P,Q))
}


#Recuperation des index de categorie film
u.factor = u.item[,c(1,6:24)]
colSums(u.factor)

#RMSE
rmse = function(m.real,m.pred) {n = sum(m.real!=0)  ;sqrt(sum((m.pred[m.real!=0] - m.real[m.real!=0])^2)/n)}

#RMSE Regularised : RSE  1
rmse.reg = function(m.real,P,Q,lambda) { m.pred= P %*% t(Q);sum((m.real[m.real!=0]-m.pred[m.real!=0])^2) + lambda * (norm(P,type = "F") + norm(Q,type = "F")) }

#Test
R = matrix(round(abs(rnorm(1:20,mean = 0, sd = 0.01)*1000)),5,4)

resultat = algo.sgd(R,19)
nP = resultat[[1]]
nQ = resultat[[2]]
nR= nP %*% t(nQ)

rmse(R,nR)
rmse.reg(R,nP,nQ,lambda)


#Algorithe SGD avec Biais
algo.sgd.bias = function(R,K,steps = 5000,gamma=0.0002,lambda=0.02) {
  H = dim(R)[1]
  W = dim(R)[2]
  P = matrix(rnorm(1:(H*K),mean = 0, sd = 0.01),H,K)
  Q = matrix(rnorm(1:(W*K),mean = 0, sd = 0.01),W,K)
  R.indices=which(R>=0,arr.ind = TRUE,useNames = TRUE)
  R.random = sample(length(R),length(R))
  R.na = R
  R.na[R==0] = NA
  average.mean = mean(R.na,na.rm = TRUE)
  bias.user = (rowSums(R.na - average.mean,na.rm = TRUE)/(rowSums(!is.na(R.na))+5))
  bias.movie = (colMeans(R.na - average.mean,na.rm = TRUE)/(colSums(!is.na(R.na))+10))
  step = 1
  while (step !=(steps*H*W)){
    r =R.indices[R.random[(step %% length(R))+1 ],]
    i = r[1]
    j = r[2]
    p = R[i,j]
    predict.value =  average.mean + bias.user[i] +bias.movie[j] + P[i,] %*% Q[j,]
    eij = p - predict.value
    average.mean = average.mean + gamma * eij
    bias.user[i] = bias.user[i] + gamma * (eij - lambda *bias.user[i] )
    bias.movie[j] = bias.movie[j] + gamma * (eij - lambda *bias.movie[j] )
    P[i,] = P[i,] + gamma * (eij * Q[j,] - lambda * P[i,])
    Q[j,] = Q[j,] + gamma * (eij * P[i,] - lambda * Q[j,])
    
    step = step +1
  }
  
  return(list(P,Q,average.mean,bias.user,bias.movie))
}

#RMSE Regularised : RSE  2 with bias
rse.2 = function(m.real,P,Q,average.mean,bias.user,bias.movie,lambda) { 
  m.pred= P %*% t(Q)
  bu = matrix(bias.user, nrow(m.real), ncol(m.real))
  bm = matrix(bias.movie, nrow(m.real), ncol(m.real),byrow = TRUE)
  average.mean = matrix(average.mean, nrow(m.real), ncol(m.real))
  sum((m.real[m.real!=0]-( bu + bm + average.mean )[m.real!=0] -m.pred[m.real!=0])^2 + lambda * (norm(P,type = "F") + norm(Q,type = "F") + (bu^2)[m.real!=0] + (bm^2)[m.real!=0] ) ) 
}

#Test
R = matrix(round(abs(rnorm(1:20,mean = 0, sd = 0.01)*1000)),5,4)

resultat = algo.sgd.bias(R,19,steps = 10000)
nP = resultat[[1]]
nQ = resultat[[2]]
naverage.mean = resultat[[3]]
nbias.user = resultat[[4]]
nbias.movie = resultat[[5]]
(nR= nP %*% t(nQ))

rmse(R,nR)
rse.2(R,nP,nQ,naverage.mean,nbias.user,nbias.movie,lambda)
