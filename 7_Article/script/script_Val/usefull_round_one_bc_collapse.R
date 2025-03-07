# Type: function

round_one_bc_collapse <- function(data, threshold, csv_file_filter){
  # Lecture du fichier de filtrage
  r6pt <- read.csv(csv_file_filter)
  r6pt <- r6pt$x
  
  # Sélection et combinaison des colonnes
  temp <- data[,r6pt]
  temp <- temp[,which(c(1:length(r6pt))%%2==1)] + temp[,which(c(1:length(r6pt))%%2==0)]
  colnames(temp) <- r6pt[which(c(1:length(r6pt))%%2==1)]
  rownames(temp) <- rownames(data)
  
  # Suppression des colonnes vides
  temp <- temp[,which(colSums(temp) > 0)]
  
  # Analyse de Zipf
  y = as.numeric(log(colSums(temp)[order(colSums(temp), decreasing = TRUE)]))
  x = log(c(1:length(y)))
  
  # Premier graphique: Distribution log-log (Zipf)
  plot(x, y, main="Distribution log-log (Zipf)", 
       xlab="Log(Rang)", ylab="Log(Somme des colonnes)",
       pch=19, col="blue")
  
  # Régression linéaire
  slope = as.numeric((y[length(y)]-y[1])/(x[length(x)]-x[1]))
  yline = slope*x + y[1]
  lines(x, yline, col="red", lwd=2)
  
  # Ajout de légende au premier graphique
  legend("topright", legend=c("Données", "Régression linéaire"),
         col=c("blue", "red"), pch=c(19, NA), lty=c(NA, 1), lwd=c(NA, 2))
  
  # Calcul des distances
  d = sqrt((y-yline)^2)
  
  # Deuxième graphique: Distances à la régression
  plot(x, d, main="Distances à la régression linéaire",
       xlab="Log(Rang)", ylab="Distance euclidienne",
       pch=19, col="black")
  
  # Points dépassant le seuil
  d2 = d[which(d >= d[which.max(d)]*threshold)]
  x2 = x[which(d >= d[which.max(d)]*threshold)]
  points(x2, d2, col="green", pch=19)
  
  # Ajout du seuil comme ligne horizontale
  abline(h=d[which.max(d)]*threshold, col="orange", lty=2)
  
  # Ajout de légende au deuxième graphique
  legend("topright", 
         legend=c("Toutes distances", "Distances > seuil", "Seuil"),
         col=c("black", "green", "orange"), 
         pch=c(19, 19, NA), 
         lty=c(NA, NA, 2))
  
  # Filtrage des colonnes selon le seuil calculé
  temp <- temp[,which(colSums(temp) >= exp(y[which(x == x2[which.max(x2)])]))]
  
  return(temp)
}

#'7_Article/script/utile_bact/r6ptorderedbcs.csv' pour csv_file path