# usefull R fonction



# Ce code convertit les noms de gènes bactériens en utilisant une table de conversion stockée dans un fichier CSV. 

convert_gene_names_bacteria <- function(genes, conversion_file){
  # genes: vecteur de noms de gènes bactériens
  # conversion_file: chemin vers le fichier CSV contenant la table de conversion


  # Lire le fichier CSV passé en argument
  conversion_table <- read.csv(conversion_file)
  new.features <- c()
  
  for (i in 1:length(genes)){
    
    if (genes[i] %in% conversion_table$genecode){
    
        genename <- conversion_table$genename[which(conversion_table$genecode == genes[i])]
        
    } else {
      
        genename <- genes[i]
        
    }
    
    new.features <- c(new.features, genename)
  }
  
  return(new.features)
}
#'7_Article/script/utile_bact/bacteria_gene_conversion.csv' pour conversion file path 



