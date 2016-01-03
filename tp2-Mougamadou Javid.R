#######################################################################################
#  TP2
#######################################################################################

# Chargement de la liste
m = read.table("http://www.cours.polymtl.ca/inf6304/Public/citeseer.rtable")

# Conversion en matrice 
m = as.matrix(m)

# Fonction de distance euclidienne pour vecteur
distance.euclidian = function(vector1,vector2){
  d = sqrt(sum((vector1 - vector2)^2))
  return(d)
}

#######################################################################################
#  Question 1
#######################################################################################

pagerank.algo = function(m,d=0.85,threeshold = 0.0001) {
  
  # Number of Article
  N <- dim(m)[1]
  # Initialisation de PR a 1/N pour tous les articles
  pr <- rep(1/N,N)
  # Application de l'algorithime avec correction
  # Correction (Google Pagerank Algorithm Correction)
  
  # Calcul de la matrice U Unitaire
  U = matrix(1,N,N)
  
  # Calcul de la matrice C ColSums
  C = 1 / (colSums(m))
  C = t(matrix(C,N,N))
  # Matrice A (Dangling Nodes Matrix)
  A = C
  A[A!=Inf]=0
  A[A==Inf]= 1/N
  
  # Calcul de la matrice H (Hyperlink Matrix)
  C[C==Inf]=0
  H = m * C
  
  # Formule du PageRank
  (G <- ( d * (H +A) ) +  (1-d)/N * U)
  
  repeat{
    
    prnew = G %*% pr
    
    # Calcul de la difference prnew - pr en valeur absolue
    difference = distance.euclidian(prnew,pr)
    # Re init
    pr = prnew
    # Condition de threeshold (par defaut 0.0001)
    if(difference <= threeshold) {
      break
      }
    
    }
    return(pr)
}

# Appel de la fonction - M doit etre transpose pour que la formule reste valide !!!!
resultat = pagerank.algo(t(m))
# Rappel des references de 422908
references =  which(m["422908",] ==1)

# Trie des references par rapport au pagerank
tri = sort(resultat[references,],decreasing = TRUE)[1:3]

# Recuperation des 3 pages recommandés par l'algorithme
which(tri>0)

### Recommandation : X311874 puis X19422 puis X17094

# Variante : References des References

m2 = m %*% m  # Adjacent Matrix Squarred (In order to determine Reference linked to References)

# Recuperation des references de longueur 2
references2 =  which(m2["422908",] >=1)

# Concatenation des references 
referencesGlobales = c(names(references),names(references2))
# Suppressions des doublons des references 
referencesGlobales = unique(referencesGlobales)

# Trie des references par rapport au pagerank
tri = sort(resultat[referencesGlobales,],decreasing = TRUE)[1:3]

# Recuperation des 3 pages recommandés par l'algorithme
which(tri>0)

### Recommandation : X311874 puis X19422 puis X17094



#######################################################################################
#  Question 2 (Item - Item)
#######################################################################################

# Calcul de la distance sur les lignes
m.distance.row = sqrt(colSums((t(m)[,"422908"] - t(m))^2))

# Calcul de la distance sur les colonnes
m.distance.column = sqrt(colSums((m[,"X422908"] - m)^2))

# Moyenne de la distance sur les lignes et sur les colonnes
m.distance = (m.distance.row + m.distance.column)/2

# Recuperation des 20 voisins les plus proches
m.voisins = names(which(sort(m.distance)[2:21]>0))

# Calcul du cosinus 
(num <- (m["422908",] %*% t(m[m.voisins,]))) # numerateurs de chaque poids
(denom <- sqrt(sum(m["422908",]^2)) * sqrt(rowSums(m[m.voisins,]^2))) # denominateurs
(wcos <- num/denom) 

(seuil <- sort(wcos[,],decreasing = TRUE)[3])# Tri Decroissant + Filtre les 3 valeurs
#Recupere les 3 articles les plus similaires
(similaire <- names(which(sort(wcos[,],decreasing = TRUE) >= seuil)))

### Recommandation : 96767 puis 149673 puis 466838