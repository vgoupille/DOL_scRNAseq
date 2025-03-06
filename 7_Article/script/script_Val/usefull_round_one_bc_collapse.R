# Type: function

round_one_bc_collapse <- function(data, threshold, csv_file_filter){
  
  r6pt <- read.csv(csv_file_filter)
  r6pt <- r6pt$x
  
  temp <- data[,r6pt]
 
  temp <- temp[,which(c(1:length(r6pt))%%2==1)] + temp[,which(c(1:length(r6pt))%%2==0)]
 
 colnames(temp) <- r6pt[which(c(1:length(r6pt))%%2==1)]
 rownames(temp) <- rownames(data)
 
 temp <- temp[,which(colSums(temp) > 0)]
 
 y = as.numeric(log(colSums(temp)[order(colSums(temp), decreasing = TRUE)]))
 x = log(c(1:length(y)))
 plot(x,y)

 slope = as.numeric((y[length(y)]-y[1])/(x[length(x)]-x[1]))
 yline = slope*x + y[1]
 lines(x,yline, col = 2)

 d = sqrt((y-yline)^2)
 plot(x,d)
 d2 = d[which(d >= d[which.max(d)]*threshold)]
 x2 = x[which(d >= d[which.max(d)]*threshold)]
 points(x2,d2, col = 3)
 
 temp <- temp[,which(colSums(temp) >= exp(y[which(x == x2[which.max(x2)])]))]
 #dim(temp)
 
 return(temp)
}

#'7_Article/script/utile_bact/r6ptorderedbcs.csv' pour csv_file path