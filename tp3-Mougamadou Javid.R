##########################################################################
#Question 1
##########################################################################

#Chargement de la librairie Rcurl
library(RCurl)
library(parallel)
library(parallelize.dynamic)

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

#Fonction Erreur Absolue Moyenne MAE
mae <- function(m.pred, m.real) mean(abs(m.pred[!is.na(m.real)]-m.real[!is.na(m.real)]))

#Fonction Erreur Quadratique  RMSE
rmse = function(m.pred, m.real) {n = sum(!is.na(m.real))  ;sqrt(sum((m.pred[!is.na(m.real)] - m.real[!is.na(m.real)])^2)/n)}

#Dimension de la matrice m
H = dim(m)[1]
W = dim(m)[2]

#Valeurs NA
m.sparse[m.sparse==0] <- NA

m.na <- m
m.na[m==0] <- NA

#Calcul de la moyenne des utilisateurs et des items
m.useraverage = rowMeans(m.sparse,na.rm = TRUE)
m.itemaverage = colMeans(m.sparse,na.rm = TRUE)

#Creation de la matrice moyenne user/item
m.useraverage2 = cbind(m.useraverage,rep(1,H))
m.itemaverage2 = cbind(rep(1,W),m.itemaverage)
m.bothaverage = (m.useraverage2 %*% t(m.itemaverage2))/2

m.bothaverage[1:5,1:5]
m[1:5,1:5]

#Check values
#all(apply(as.matrix(1:943),1,function(i) all((m.useraverage[i]+m.itemaverage)/2==m.bothaverage[i,])))

#Calcul du  MAE
(mae(m.bothaverage,m.na)) # = 0.78

#Calcul du RMSE
(rmse(m.bothaverage,m.na)) #N = 0.96

##########################################################################
#Question 2
##########################################################################

#Normalization
m.itemaverage.matrix = matrix(m.itemaverage,H,W,byrow = TRUE)
m.filled = m
m.filled[m==0] = m.itemaverage.matrix[m==0]
m.normalized = m.filled - m.useraverage

#Calcul SVD  (M = U . D . V)
m.svd = svd(m.normalized)

##########################################################################
#Question 3
##########################################################################

#Nombre de dimension k = 10
k = 10

#Matrice Diagonale D - Conservation des 10 premieres valeurs
d = m.svd$d
d.k = d[1:k]

#Matrice D Squarred 
d.sq = sqrt(d.k)

#Matrice U reduit à 10 dimensions => Uk
u = m.svd$u
u.k = u[,1:k]

#Matrice V reduit à 10 dimensions => Vk
v = m.svd$v
v.k = v[,1:k]

#Matrice des moyennes des utilisateurs
C = matrix(rowMeans(m.na,na.rm = TRUE),H,W)

#Matrice de Prediction P = C + Uk.tSk.Sk.tVk
P = C + (u.k %*% diag(d.sq,k)) %*% (diag(d.sq,k) %*% t(v.k))

##########################################################################
#Question 4
##########################################################################

#Calcul du MAE
(mae(P,m.na)) # = 0.70

#Calcul du RMSE
(rmse(P,m.na)) # = 0.88


##########################################################################
#Question 5
##########################################################################


#Algortithme pour le calcul du SVD + RMSE en fonction de la matrice m et du nombre de dimension k
algo.svd = function(k,m.na,m.svd) {
  
  d = m.svd$d
  d.k = d[1:k]
  
  d.sq = sqrt(d.k)
  
  u = m.svd$u
  u.k = u[,1:k]
  
  v = m.svd$v
  v.k = v[,1:k]
  
  C = matrix(rowMeans(m.na,na.rm = TRUE),H,W)
  
  P = C + (u.k %*% diag(d.sq,k)) %*% (diag(d.sq,k) %*% t(v.k))
  
  return(P)
  
}


#Generation des test des dimensions de 1 à 100
m.na <- m
m.na[m==0] <- NA
testrange = c(1,5,10,14,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,943)
predict=apply(as.matrix(testrange),1,function(i){ 
    P= algo.svd(i,m.na,m.svd)
    return(cbind(mae(P,m.na),rmse(P,m.na)))
    } )
colnames(predict) = testrange

#Affichage de la courbe du MAE en fonction du nombre de dimensions
plot(testrange,predict[1,],type = "l")

#Affichage de la courbe du RMSE en fonction du nombre de dimensions
plot(testrange,predict[2,],type = "l")

#Recuperation du nombre de dimension optimal
k.optimal = sort(predict[1,])[1] 
k.optimal = which(predict[1,] <= k.optimal) # K Optimal = 943
(predict[1,k.optimal])    # MAE K Optimal = 5.006641e-15

k.optimal = sort(predict[2,])[1] 
k.optimal = which(predict[2,] <= k.optimal) # K Optimal = 943
(predict[2,k.optimal])    # RMSE K Optimal = 1.096776e-14

##########################################################################
#Question 6   k=943
##########################################################################

## Index aléatoire des données de tests
i.observed <- which(m > 0)
i.hasard <- sample(i.observed, length(i.observed))
length(i.hasard)
fold.size <- round(length(i.hasard) / 10)
i.false <- rep(FALSE, length(m))

predict = lapply(as.matrix(1:10),function(i){


fold.number <- i

## Index booléen pour les cellules de test et d'entraînement
i.test.b <- i.false
## Les cellules indexées du replis correspondant sont fixées à TRUE pour le test...
i.test.b[ i.hasard[((fold.number-1) * fold.size):((fold.number) * fold.size)] ] <- TRUE
## ...et à FALSE pour l'entraînement
i.train.b <-  !i.test.b
m.na.train <- m.na
m.na.train[i.test.b] <- NA                # remove test data from training
table(m.na.train)
votes.films.moyens <- colMeans(m.na.train, na.rm=T)
mean(votes.films.moyens)                # des NaN pourraient être créés car certains films n'ont plus aucun vote
## Il faudrait alors remplacer ces colonnes par une valeur correspondant à la moyenne générale.
moy.globale <- mean(m.na.train, na.rm=T)
films.sans.votes <- colSums(m.na.train, na.rm=T) == 0
sum(films.sans.votes)                   # si 0 alors pas besoin de faire l'ajustement suivant
m.na.train[,films.sans.votes] <- moy.globale
votes.films.moyen <- colMeans(m.na.train, na.rm=T)
## fin de l'ajustement
hist(votes.films.moyens)
## votes moyens des utilisateurs de test
votes.utilisateurs.moyen <- rowMeans(m.na.train, na.rm=T)
## Normalisation

votes.films.moyen.matrix = matrix(votes.films.moyen,H,W,byrow = TRUE)
m.na.train.filled = m.na.train
m.na.train.filled[is.na(m.na.train.filled)] = votes.films.moyen.matrix[is.na(m.na.train.filled)]
m.na.train.normalized = m.na.train.filled - votes.utilisateurs.moyen

m.na.train.svd = svd(m.na.train.normalized)
  
P= algo.svd(k.optimal,m.na,m.na.train.svd)
  
## Erreur absolue moyenne
e1=mean(abs(P[i.test.b] - m.na[i.test.b]), na.rm=T)
## Racine carrée de erreur quadratique moyenne
e2=sqrt(mean((P[i.test.b] - m.na[i.test.b])^2, na.rm=T))
  
return(cbind(e1,e2))

})

#cl = makeCluster(detectCores())
#system.time(parLapply(cl,1:100,algo.svd,m.na,m.svd, H, W))
#system.time(lapply(1:100,algo.svd,m.na,m.svd, H, W))

#a= lapply(lapply(predict,'[[',1),sum)

#Mean of the 10 mae values
mae.mean=(predict[[1]][1]+predict[[2]][1]+predict[[3]][1]+predict[[4]][1]+predict[[5]][1]+
  predict[[6]][1]+predict[[7]][1]+predict[[8]][1]+predict[[9]][1]+predict[[10]][1])/10

#Mean of the 10 rmse values
rmse.mean=(predict[[1]][2]+predict[[2]][2]+predict[[3]][2]+predict[[4]][2]+predict[[5]][2]+
            predict[[6]][2]+predict[[7]][2]+predict[[8]][2]+predict[[9]][2]+predict[[10]][2])/10


##########################################################################
#Question 7   Find Optimal k by Cross Validation
##########################################################################

## Index aléatoire des données de tests
i.observed <- which(m > 0)
i.hasard <- sample(i.observed, length(i.observed))
length(i.hasard)
fold.size <- round(length(i.hasard) / 10)
i.false <- rep(FALSE, length(m))

predict = lapply(as.matrix(1:10),function(i){
  
  
    fold.number <- i
    
    ## Index booléen pour les cellules de test et d'entraînement
    i.test.b <- i.false
    ## Les cellules indexées du replis correspondant sont fixées à TRUE pour le test...
    i.test.b[ i.hasard[((fold.number-1) * fold.size):((fold.number) * fold.size)] ] <- TRUE
    ## ...et à FALSE pour l'entraînement
    i.train.b <-  !i.test.b
    m.na.train <- m.na
    m.na.train[i.test.b] <- NA                # remove test data from training
    table(m.na.train)
    votes.films.moyens <- colMeans(m.na.train, na.rm=T)
    mean(votes.films.moyens)                # des NaN pourraient être créés car certains films n'ont plus aucun vote
    ## Il faudrait alors remplacer ces colonnes par une valeur correspondant à la moyenne générale.
    moy.globale <- mean(m.na.train, na.rm=T)
    films.sans.votes <- colSums(m.na.train, na.rm=T) == 0
    sum(films.sans.votes)                   # si 0 alors pas besoin de faire l'ajustement suivant
    m.na.train[,films.sans.votes] <- moy.globale
    votes.films.moyen <- colMeans(m.na.train, na.rm=T)
    ## fin de l'ajustement
    hist(votes.films.moyens)
    ## votes moyens des utilisateurs de test
    votes.utilisateurs.moyen <- rowMeans(m.na.train, na.rm=T)
    ## Normalisation
    
    votes.films.moyen.matrix = matrix(votes.films.moyen,H,W,byrow = TRUE)
    m.na.train.filled = m.na.train
    m.na.train.filled[is.na(m.na.train.filled)] = votes.films.moyen.matrix[is.na(m.na.train.filled)]
    m.na.train.normalized = m.na.train.filled - votes.utilisateurs.moyen
    
    m.na.train.svd = svd(m.na.train.normalized)
    
    comparaison=apply(as.matrix(testrange),1,function(j) {
      P= algo.svd(j,m.na,m.na.train.svd)
      
      ## Erreur absolue moyenne
      e1=mean(abs(P[i.test.b] - m.na[i.test.b]), na.rm=T)
      
      return(e1)
      
      })
    names(comparaison)=testrange
    k.optimal = sort(comparaison)[1] 
    k.optimal = which(comparaison <= k.optimal) 
    (comparaison[k.optimal]) 
    return(cbind(k.optimal,comparaison[k.optimal]))
    
  
})

predict

#  Best Dimension k = 14
#  Best MAE for k=14 : 0.7377675 
