convert_gene_names_bacteria <- function(genes){
  
  conversion_table <- read.csv('C:/Users/lbrettne/Desktop/bacteria_gene_conversion.csv')
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


# Ce code convertit les noms de gènes bactériens en utilisant une table de conversion stockée dans un fichier CSV. Voici une analyse détaillée de son fonctionnement :

# 🔍 Fonctionnement du code

# 1️⃣ Chargement de la table de conversion

# conversion_table <- read.csv('C:/Users/lbrettne/Desktop/bacteria_gene_conversion.csv')

# 	•	La table de conversion bacteria_gene_conversion.csv est chargée dans conversion_table.
# 	•	On suppose qu’elle contient au moins deux colonnes :
# 	•	genecode : le code du gène (ancien identifiant).
# 	•	genename : le nom du gène (nouvel identifiant).

# Exemple de table bacteria_gene_conversion.csv :

# genecode	genename
# b0001	thrL
# b0002	thrA
# b0003	thrB
# b0004	thrC

# 2️⃣ Initialisation d’un vecteur vide

# new.features <- c()

# 	•	new.features est un vecteur vide qui stockera les nouveaux noms de gènes.

# 3️⃣ Boucle de conversion

# for (i in 1:length(genes)){

# 	•	On parcourt chaque élément du vecteur genes, qui contient la liste des identifiants à convertir.

# 3.1 Vérification si le gène est dans la table

# if (genes[i] %in% conversion_table$genecode){

# 	•	On vérifie si le gène actuel (genes[i]) se trouve dans la colonne genecode de la table de conversion.

# 3.2 Conversion du nom du gène

# genename <- conversion_table$genename[which(conversion_table$genecode == genes[i])]

# 	•	Si le gène est trouvé, on récupère le nom correspondant dans genename.

# 3.3 Conservation du nom d’origine si non trouvé

# } else {
#     genename <- genes[i]
# }

# 	•	Si genes[i] n’existe pas dans la table de conversion, on le garde tel quel.

# 3.4 Ajout du nom dans le vecteur new.features

# new.features <- c(new.features, genename)

# 	•	Le nom (converti ou non) est ajouté à new.features.

# 4️⃣ Retour des gènes convertis

# return(new.features)

# 	•	La fonction retourne la liste des noms de gènes convertis.

# 🧐 Exemple d’utilisation

# Cas 1 : Liste de gènes avec conversions disponibles

# genes <- c("b0001", "b0002", "b9999")  # "b9999" n'existe pas dans la table
# convert_gene_names_bacteria(genes)

# Table de conversion

# genecode	genename
# b0001	thrL
# b0002	thrA

# Résultat attendu

# [1] "thrL" "thrA" "b9999"

# 	•	b0001 → thrL
# 	•	b0002 → thrA
# 	•	b9999 → pas dans la table, donc inchangé

# 📌 Résumé : Que fait cette fonction ?

# ✅ Charge une table de conversion des gènes bactériens.
# ✅ Remplace chaque identifiant de gène (genecode) par son nom (genename).
# ✅ Conserve les gènes non trouvés tels quels.
# ✅ Retourne la liste des noms de gènes convertis.

# 📌 Contexte d’application :
# 	•	Annotation génomique : conversion des codes de gènes en noms lisibles.
# 	•	Analyse de transcriptome/metagenome : standardisation des identifiants de gènes bactériens.