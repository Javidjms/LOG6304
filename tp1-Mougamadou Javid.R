##########################################################################
#Question 1
##########################################################################

#Chargement de la librairie Rcurl
library(RCurl)

#Chargement des données
u.data <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.data.csv', userpwd='20113:20113'), sep='|', header=T)
u.user <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.user.csv', userpwd='20113:20113'), sep='|', header=T)
u.item <- read.csv(text=getURL('http://www.groupes.polymtl.ca/log6304/Tp/20153/u.item.csv', userpwd='20113:20113'), sep='|', header=T)

#Jointure du tableau u.data avec u.user sur la colonne user.id
v.u <- merge(u.data, u.user, by.x='user.id', by.y='id')
#head(v.u)

#Moyenne des votes par Age
(aggregate(v.u$rating,list(v.u$age),mean))
#plot(agg.age$Age,agg.age$rating)

#Moyenne des votes par Profession
(aggregate(v.u$rating,list(v.u$job),mean))
#plot(agg.job$Job,agg.job$rating)

##########################################################################
#Question 2
##########################################################################

# Star Trek V: The Final Frontier (1989)  =====>  i450

#Chargement de la librairie Matrix
library(Matrix)

#Transformation en matrice (pour certaines opérations) avec les entetes des lignes u et colonnes i
m <- sparseMatrix(u.data[,1],u.data[,2],x=u.data[,3])
m.sparse <- sparseMatrix(u.data[,1],u.data[,2],x=u.data[,3])
rownames(m.sparse) <- paste('u', 1:nrow(m.sparse), sep='')
colnames(m.sparse) <- paste('i', 1:ncol(m.sparse), sep='')
m <- as.matrix(m.sparse)

m.sparse<-m
m.sparse[m.sparse==0] <- NA

#Tranposition de la matrice
tm = t(m)

## Calculs avec le cosinus
(num <- (tm[450,] %*% t(tm[-450,]))) # numerateurs de chaque poids
(denom <- sqrt(sum(tm[450,]^2)) * sqrt(rowSums(tm[-450,]^2))) # denominateurs
(wcos <- num/denom) #Cosinus entre i450 et les autres i

#Determination de la valeur seuil entre les 10 premier film les plus similaires
(seuil <- sort(wcos,decreasing = TRUE)[10])# Tri Decroissant + Filtre les 10 valeurs
#Recupere les 10 premiers indices des 10 films les plus similaires
(similaire <- which(wcos >= seuil)) 
#Recuperation des noms des 10 films les plus similaires
u.item[similaire,2]


##Calcul avec la corrélation
#Determination de la corelation entre i450 et les autres
tm.na = tm
tm.na[tm==0]=NA
wcov = cor(t(tm.na),tm.na[450,],use="pairwise.complete.obs")

#Determination de la valeur seuil entre les 10 premier film les plus similaires
(seuil <- sort(wcov,decreasing = TRUE)[10])

#Recupere les 10 premiers indices des 10 films les plus similaires
(similaire <- which(wcov >= seuil))

#Recuperation des noms des 10 films les plus similaires
u.item[similaire,2][1:10]

##########################################################################
#Question 3
##########################################################################

#Fonctions utiles 

min.nindex = function(m,n=5){
  i = order(m)
  return (i[1:n])
}

max.nindex = function(m,n=5){
  i = order(m,decreasing = TRUE)
  return (i[1:n])
}

#Calcul de la distance euclidienne
distance.450 = sqrt(colSums((m[,450] - m)^2))
# Calcul des votes communs
votes.communs <- (tm[450,]>0) %*% t(tm>0)

# Utilisateurs qui ont voté pour le film 450
uservoted = tm[450,]!="0"

# Utilisateurs qui n'ont pas voté pour le film 450
usernotvoted = tm[450,]=="0"

# On recupère les 40 voisins les plus proches de i450
i.distance.450 = min.nindex(distance.450,40)

# Verification des votes communs entre voisins
votes.communs[i.distance.450]

# Selection des voisins qui n'ont pas de votes en communs (0 ou 1)
i.nonfiable <- which(votes.communs[i.distance.450] <= 1)

# Selection des 20 voisins fiable qui ont aux moins 2 votes en communs
i.distance.450.fiable = i.distance.450[-i.nonfiable]
i.voisins.450 = i.distance.450.fiable[-1]
votes.communs[i.voisins.450]

# Calcul de la matrice centrée
(v.i <- rowMeans(tm[i.voisins.450,]))
(tm.centre <- tm[i.voisins.450,] - v.i)

# Calcul du cosinus
(num <- (tm[450,] %*% t(tm[i.voisins.450,]))) # numerateurs de chaque poids
(denom <- sqrt(sum(tm[450,]^2)) * sqrt(rowSums(tm[i.voisins.450,]^2))) # denominateurs
(wcos <- num/denom) #Cosinus entre i450 et les autres i

# Estimation du vote sur la base des voisins choisis
vote.estimation = mean(tm[450,uservoted]) + 1/sum(abs(wcos)) * (wcos %*% tm.centre)

# Visualisation des estimations de votes

table(vote.estimation)

##########################################################################
#Question 4
##########################################################################

# Nombre d'utilisateur ayant voté pour le film i450
n = sum(uservoted)

# Calcul de l'erreur quadratique RMSE
rmse = sqrt(sum((tm[450,uservoted] - vote.estimation[uservoted])^2)/n)


##########################################################################
##########################################################################

#Question 5  (Version Revisé du chargé de TP Peng Xu)

##########################################################################
##########################################################################

# Recuperation des indices des films de Star Trek et de Star Wars
indices.star.trek = grep("trek", as.character(u.item$movie.title), ignore.case=T)
u.item$movie.title[indices.star.trek]
indices.star.wars = c(172,181)
u.item$movie.title[indices.star.wars]
indices.star = c(indices.star.trek,indices.star.wars)

# Creation du nouvel utilisateur
newuser = rep.int(0, 1682)
newuser[indices.star.trek] = 5
newuser[indices.star.wars] = 1

# Calcul de la distance euclidienne

# First, calculate the distance by using m.sparse, thus only 10 columns are used for calculation
newuser.distance = sqrt(rowSums(t(t(m.sparse[,indices.star])-newuser[indices.star])^2,na.rm=TRUE))
# Find the number of common votes for each user
votes.communs <- (newuser>0) %*% t(m > 0)
# Set all the distance of users with no common votes to be NA
newuser.distance[votes.communs==0] <- NA

# Obtain the 20 neighbours for newUser
newuser.voisins <- min.nindex(newuser.distance,20)

# Calculate the cosine value
# This part remains the same
num <- newuser %*% t(m[newuser.voisins, ])
denom <- sqrt(sum(newuser ^ 2)) * sqrt(rowSums(m[newuser.voisins, ]^2))
wcos <- num/denom

# Calculate the predicted votes
# Pay attention here, this part is where the errors come from.
# Obtain the rowMeans, notice that we should use m.sparse
v.i <- rowMeans(m.sparse[newuser.voisins, ], na.rm=T)
# Centering
m.centre <- (m.sparse[newuser.voisins,]) - v.i
# Assign all the NA values to be 0
m.centre[is.na(m.centre)] <- 0
# Average vote of the newUser
u.average <- mean(newuser[indices.star])
# Apply the user-user formula
vote.estimation.newuser <- u.average + 1/sum(abs(wcos)) * (wcos %*% m.centre)

# Have a look at the result
hist(vote.estimation.newuser)
# Remove the dim attribute of vote.estimation.newUser
dim(vote.estimation.newuser)<-NULL
# Combine relavant information in a dataframe
result<-data.frame("id" = u.item$movie.id,
                   "title" = u.item$movie.title,
                   "Prediction" = vote.estimation.newuser)
# Remove the items which the newUser has already rated
result<-result[-indices.star,]
# Final recommendation
recommend<-head(result[order(result$Prediction,decreasing=T),],10)
recommend
# Check the common votes of these movies
m.sparse[newuser.voisins,recommend$id]

##########################################################################
#Question 6
##########################################################################

# Fonctions utiles

neg = function(p) {1-p }
Odds = function(p) { p/(1-p) }
OddsToP = function(o) { o/(1+o) }

ratio.chances = function(rating.vec,seuil=3){
  sum(rating.vec > seuil)/(1+sum(rating.vec <= seuil && rating.vec >0))
}

# Jointure entre u.user et u.date
mx = merge(u.user,u.data,1)

# Utilisation de notre matrice tm
tm

# I : Vecteur avec le nombre des utilisateurs qui aime le film par film 
i = rowSums(tm > 3)

# NI : Vecteur avec le nombre des utilisateurs qui n'aime le film par film 
ni = rowSums(tm <= 3 & tm >0)

# O(H) Chances Initiales avec Correction pour eviter les valeurs infini 0/0 ou 1/0
O.h = (1+i) /(1+ni)

# P(H)
P.h = OddsToP(O.h)

# Creation de notre algorithme pour un nouvel utilisateur

algo.newuser.recommandation = function(job,gender,age){
  
  # Recuperation des utilisateurs ayant pour job celui donnée en paramètre
  E1.id = mx[mx$job==job,'id']
  
  # Recuperation des utilisateurs sans doublons
  E1.id = unique(E1.id)
  
  # P(E_1 |H) : probabilité qu'il aime le film sachant qu'il a le job "job"
  E1.like = rowSums(tm[,E1.id] > 3)/(rowSums(tm[,E1.id])+1)
  P.E1.like = (E1.like)/(P.h)
  
  # P(E_1 |~H) : probabilité qu'il n'aime le film sachant qu'il a le job "job"
  E1.dislike = rowSums((tm[,E1.id] <= 3) & (tm[,E1.id] > 0))/(rowSums(tm[,E1.id])+1)
  P.E1.dislike = (E1.dislike)/(neg(P.h))
  
  # Recuperation des utilisateurs ayant pour sexe celui donnée en paramètre
  E2.id = mx[mx$gender==gender,'id']
  
  # Recuperation des utilisateurs sans doublons
  E2.id = unique(E2.id)
  
  # P(E_2 |H) : probabilité qu'il aime le film sachant qu'il a le sexe "gender"
  E2.like = rowSums(tm[,E2.id] > 3)/(rowSums(tm[,E2.id])+1)
  P.E2.like = (E2.like)/(P.h)
  
  # P(E_2 |~H) : probabilité qu'il n'aime le film sachant qu'il a le sexe "gender"
  E2.dislike = rowSums((tm[,E2.id] <= 3) & (tm[,E2.id] > 0))/(rowSums(tm[,E2.id])+1)
  P.E2.dislike = (E2.dislike)/(neg(P.h))
  
  # Recuperation par tranche d'age
  # Recuperation des utilisateurs ayant pour age celui donnée en paramètre
  if (age=='child') {
    E3.id = mx[mx$age<20,'id']
  }
  
  else if (age=='young') {
    E3.id = mx[mx$age>=20 & mx$age<50,'id']
  }

  else if (age=='old') {
    E3.id = mx[mx$age>=50,'id']
  }
  
  # Recuperation des utilisateurs sans doublons
  E3.id = unique(E3.id)
  
  # P(E_3 |H) : probabilité qu'il aime le film sachant qu'il a l'age "age"
  E3.like = rowSums(tm[,E3.id] > 3)/(rowSums(tm[,E3.id])+1)
 
  P.E3.like = (E3.like)/(P.h)
  
  # P(E_3 |~H) : probabilité qu'il n'aime le film sachant qu'il a l'age "age"
  E3.dislike = rowSums((tm[,E3.id] <= 3) & (tm[,E3.id] > 0))/(rowSums(tm[,E3.id])+1)
  P.E3.dislike = (E3.like)/(neg(P.h))
  
  # O(H|E_1, E_2) = O(H) \frac{P(E_1 | H)}{P(E_1 | ~H)} * \frac{P(E_2 | H)}{P(E_2 | ~H)}
  O.All_E = O.h * ((1 +P.E1.like) /(1+P.E1.dislike)) * ((1 +P.E2.like) /(1+P.E2.dislike)) * ((1 +P.E3.like) /(1+P.E3.dislike)) 

  # Conversion Chance (Odds) to Probabilité
  P.All_E = OddsToP(O.All_E)
  
  # Liste des 10 films recommandés par notre algorithme bahésien
  (u.item[max.nindex(P.All_E,10),2])
  
  
}

# Example to Use

algo.newuser.recommandation('engineer','M','old')

