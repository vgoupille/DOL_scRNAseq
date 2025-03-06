
# Ce code R définit une fonction nommée `assign_cell_wells` qui associe des codes-barres (barcodes) à des positions spécifiques dans une plaque de laboratoire (puits ou "wells").

# Voici ce que fait cette fonction :

# 1. Elle prend en entrée un vecteur de codes-barres (`barcodes`)
# 2. Elle initialise un dataframe vide `cell_wells` avec deux colonnes (barcode et well)
# 3. Pour chaque code-barre dans le vecteur d'entrée :
#    - Elle extrait la dernière partie du code-barre en le divisant par le caractère "_"
#    - Elle vérifie cette séquence (BC1) par rapport à une longue liste de séquences connues
#    - Selon la séquence identifiée, elle attribue le code-barre à un puits spécifique (A1, A2, B1, etc.)
#    - Si aucune correspondance n'est trouvée, elle attribue "unknown" au puits
#    - Elle affiche le numéro d'itération (i) pour suivre la progression
# 4. À la fin, elle renvoie le dataframe complet avec les associations code-barre/puits

# Cette fonction est typiquement utilisée dans des analyses de séquençage à haut débit où les échantillons sont identifiés par des codes-barres et disposés dans une plaque de 96 puits (format 8×12, d'où les désignations A1-H12).




assign_cell_wells <- function(barcodes){

  cell_wells <- data.frame(barcode <- c(), well <- c())
  
  for (i in 1:length(barcodes)){
    
    BC1 <- last(unlist(strsplit(barcodes[i], "_")))
    
    if (BC1 == 'ACTCGTAA' | BC1 == 'CTGCTTTG'){
      temp <- data.frame(barcode = barcodes[i], well = "A1")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'AAACGATA' | BC1 == 'CATGATCA'){
      temp <- data.frame(barcode = barcodes[i], well = "A2")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTACCTCG' | BC1 == 'GGGTAGCG'){
      temp <- data.frame(barcode = barcodes[i], well = "A3")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GCCTGCAA' | BC1 == 'CCGAGAAA'){
      temp <- data.frame(barcode = barcodes[i], well = "A4")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TGGTATAC' | BC1 == 'ACGGACTC'){
      temp <- data.frame(barcode = barcodes[i], well = "A5")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CGTTCGAG' | BC1 == 'ACTTACGA'){
      temp <- data.frame(barcode = barcodes[i], well = "A6")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TCTATTAC' | BC1 == 'TATTTAAG'){
      temp <- data.frame(barcode = barcodes[i], well = "A7")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ATAAGCTC' | BC1 == 'ACCGTACG'){
      temp <- data.frame(barcode = barcodes[i], well = "A8")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ATTCATGG' | BC1 == 'TATAGTCG'){
      temp <- data.frame(barcode = barcodes[i], well = "A9")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ATCCGCGA' | BC1 == 'TGGGCATC'){
      temp <- data.frame(barcode = barcodes[i], well = "A10")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ATCGCATA' | BC1 == 'TACCTAGA'){
      temp <- data.frame(barcode = barcodes[i], well = "A11")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CCGTTCTA' | BC1 == 'GCTGCATG'){
      temp <- data.frame(barcode = barcodes[i], well = "A12")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TGGCGCGC' | BC1 == 'GTCATATG'){
      temp <- data.frame(barcode = barcodes[i], well = "B1")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TGTCTGAA' | BC1 == 'ATATTGGC'){
      temp <- data.frame(barcode = barcodes[i], well = "B2")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CTGTCCCG' | BC1 == 'CTAAGGGA'){
      temp <- data.frame(barcode = barcodes[i], well = "B3")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'AATTTCTC' | BC1 == 'TCGTTTCG'){
      temp <- data.frame(barcode = barcodes[i], well = "B4")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CGCGACTA' | BC1 == 'GAATAATG'){
      temp <- data.frame(barcode = barcodes[i], well = "B5")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GGGATCGG' | BC1 == 'ACTGCGCA'){
      temp <- data.frame(barcode = barcodes[i], well = "B6")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTATTCTG' | BC1 == 'GCTTATAG'){
      temp <- data.frame(barcode = barcodes[i], well = "B7")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'AGGCGGCA' | BC1 == 'ATCATGCA'){
      temp <- data.frame(barcode = barcodes[i], well = "B8")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ACGCCGGC' | BC1 == 'ACGTTAAC'){
      temp <- data.frame(barcode = barcodes[i], well = "B9")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTGTCTTA' | BC1 == 'CCATCTTG'){
      temp <- data.frame(barcode = barcodes[i], well = "B10")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TACGGTTA' | BC1 == 'CATAGCTA'){
      temp <- data.frame(barcode = barcodes[i], well = "B11")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTGGGAGA' | BC1 == 'GAGGTTGA'){
      temp <- data.frame(barcode = barcodes[i], well = "B12")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TGCTTGGG' | BC1 == 'GCACTGAC'){
      temp <- data.frame(barcode = barcodes[i], well = "C1")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TAAATATC' | BC1 == 'TTCATCGC'){
      temp <- data.frame(barcode = barcodes[i], well = "C2")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CACAATTG' | BC1 == 'GAAATTAG'){
      temp <- data.frame(barcode = barcodes[i], well = "C3")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GTGCTAGC' | BC1 == 'AGGATTAA'){
      temp <- data.frame(barcode = barcodes[i], well = "C4")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CGCCCGGA' | BC1 == 'AATAGAAC'){
      temp <- data.frame(barcode = barcodes[i], well = "C5")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GCTCGCGG' | BC1 == 'TCTTAATC'){
      temp <- data.frame(barcode = barcodes[i], well = "C6")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CTTTGGTC' | BC1 == 'TAATACGC'){
      temp <- data.frame(barcode = barcodes[i], well = "C7")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTCCGATC' | BC1 == 'GTTTGTGA'){
      temp <- data.frame(barcode = barcodes[i], well = "C8")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TTCGCTAC' | BC1 == 'CGAACGTC'){
      temp <- data.frame(barcode = barcodes[i], well = "C9")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'AGCGAAAC' | BC1 == 'GGTTCTTC'){
      temp <- data.frame(barcode = barcodes[i], well = "C10")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'AAATAGCA' | BC1 == 'GCAAATTC'){
      temp <- data.frame(barcode = barcodes[i], well = "C11")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CGTCTAGG' | BC1 == 'GCTATGCG'){
      temp <- data.frame(barcode = barcodes[i], well = "C12")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GCCGTGTA' | BC1 == 'CTACCCTA'){
      temp <- data.frame(barcode = barcodes[i], well = "D1")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CGCTTAAA' | BC1 == 'GTGGGTTC'){
      temp <- data.frame(barcode = barcodes[i], well = "D2")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GACCTTTC' | BC1 == 'GTCCGTAG'){
      temp <- data.frame(barcode = barcodes[i], well = "D3")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GGTGGAGC' | BC1 == 'TGCGATCG'){
      temp <- data.frame(barcode = barcodes[i], well = "D4")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TACTCGAA' | BC1 == 'TATCCGGG'){
      temp <- data.frame(barcode = barcodes[i], well = "D5")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CATTTGGA' | BC1 == 'AGGTAATA'){
      temp <- data.frame(barcode = barcodes[i], well = "D6")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GACGGGAC' | BC1 == 'CGTGGTTG'){
      temp <- data.frame(barcode = barcodes[i], well = "D7")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GTCGCGCG' | BC1 == 'GACAAAGC'){
      temp <- data.frame(barcode = barcodes[i], well = "D8")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'GTTACGTA' | BC1 == 'GGGCGATG'){
      temp <- data.frame(barcode = barcodes[i], well = "D9")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'CTATTTCA' | BC1 == 'ATCTATAA'){
      temp <- data.frame(barcode = barcodes[i], well = "D10")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'ACTATATA' | BC1 == 'GCCCATGA'){
      temp <- data.frame(barcode = barcodes[i], well = "D11")
      cell_wells <- rbind(cell_wells, temp)
    } else if (BC1 == 'TCACTTTA' | BC1 == 'CTGAAAGG'){
      temp <- data.frame(barcode = barcodes[i], well = "D12")
      cell_wells <- rbind(cell_wells, temp)
    } else {
      temp <- data.frame(barcode = barcodes[i], well = "unknown")
      cell_wells <- rbind(cell_wells, temp) 
    }
    #print(i)
  }

  return(cell_wells)
}