library(ggplot2)

hour  <- c(0.00,1.10,2.13,2.45,2.75,2.90,3.18,3.28,3.60,3.93,4.0,4.42,4.67,4.95,6.12,6.97)
OD600 <- c(0.02,0.04,0.27,0.49,0.94,0.93,1.16,1.34,1.58,1.76,2.0,2.87,3.20,3.56,5.31,6.34)

M14hour <- c(4.0, 4.67)
M14samp <- c(2.0, 3.2)
M15hour <- c(3.28,4.95,6.12,6.97)
M15samp <- c(1.34,3.56,5.31,6.34)
M45hour <- c(2.45,2.75,3.93,4.42)
M45samp <- c(0.49,0.94,1.76,2.87)

m <- c(1.34,3.56,5.31,6.34,0.49,0.94,1.76,2.87)

plot(hour, log(OD600), type = "l"#, main = "B. subtilis growth curve (M15)"
     )
points(hsamp,log(ODsamp),col = 2)

fit <- nls(OD600 ~ L/(1+exp(-k*(hour - x0))), start = list(L = 6, k = 0.5, x0 = 3))

x <- seq(0,10,by=0.01)
y <- 7.2/(1+exp(-0.95*(x-4.96)))
plot(x,y, xlab = "hour", ylab = "OD600", type = "l")
points(M14hour,M14samp)
points(M15hour,M15samp,col = 2)
points(M45hour,M45samp,col = 3)

fit <- data.frame(x = x, y = y)

time <- c(4.0,4.67,3.28,4.95,6.12,6.97,2.45,2.75,3.93,4.42)
OD   <- c(2.0,3.2,1.34,3.56,5.31,6.34,0.49,0.94,1.76,2.87)
exp  <- c("experiment 1","experiment 1","experiment 2","experiment 2","experiment 2","experiment 2","combined","combined","combined","combined")

data <- data.frame(hour = time, OD600 = OD, experiment = exp)

#pink  #e39fc0 experiment 2
#blue  #383a73 experiment 1
#green #5eccad combined

p <- ggplot() +
      geom_line(data = fit, aes(x = x, y = y), size = 1) +
      geom_point(data = data, aes(x = hour, y = OD600, fill = "black"), colour = "black", pch = 21, size = 5) +
      #scale_fill_manual(values = c("#5eccad", "#383a73", "#e39fc0")) +
      theme_bw() +
      theme(axis.text = element_text(size = 16),
        text = element_text(size = 16),
        legend.text = element_text(size = 16),
        aspect.ratio = 1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
      xlab("hour") +
      ylab("OD600") +
      NoLegend() +
      title("bacillus subtilis")
      

p

ggsave(
  "ggtest.png",
  p,
  dpi = 1200
)

p <- ggplot() +
  geom_line(data = fit, aes(x = x, y = y), size = 1) +
  geom_point(data = data, aes(x = hour, y = OD600, fill = exp), colour = "black", pch = 21, size = 5) +
  scale_fill_manual(values = c("#5eccad", "#383a73", "#e39fc0")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16),
        text = element_text(size = 16),
        legend.text = element_text(size = 16),
        aspect.ratio = 1,
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.8,0.15),
        legend.title = element_blank()) +
  xlab("hour") +
  ylab("OD600")


p

doublingrate <- function(fit, x, t){
  drs <- c()
  for (i in 2:length(fit)){
    dr <- log2(fit[i]/fit[i-1])/(x[i]-x[i-1])
    drs <- c(drs,dr)
  }
  
  grs <- c()
  for (j in 1:length(t)){
    k <- which.min(abs(x - t[j]))
    gr <- drs[k]
    grs <- c(grs, gr)
  }
  return(grs)
}

dt <- doublingrate(y,x,hour)