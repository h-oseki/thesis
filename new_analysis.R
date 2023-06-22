setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 12),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)
##
#### ABIERTA 1 ####
ab1_res190 <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_res190_dist.xvg")
ab1_res190 <- as.data.frame(sapply(ab1_res190, as.numeric))
colnames(ab1_res190) <- c("time", "distance")
ab1_a3a <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_alpha3A.xvg")
ab1_a3a <- as.data.frame(sapply(ab1_a3a, as.numeric))
colnames(ab1_a3a) <- c("time", "x", "y", "z")
ab1_a3b <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_alpha3B.xvg")
ab1_a3b <- as.data.frame(sapply(ab1_a3b, as.numeric))
colnames(ab1_a3b) <- c("time", "x", "y", "z")
ab1_a4a <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_alpha4A.xvg")
ab1_a4a <- as.data.frame(sapply(ab1_a4a, as.numeric))
colnames(ab1_a4a) <- c("time", "x", "y", "z")
ab1_a4b <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_alpha4B.xvg")
ab1_a4b <- as.data.frame(sapply(ab1_a4b, as.numeric))
colnames(ab1_a4b) <- c("time", "x", "y", "z")
ab1_prot <- readXVG("abierta-rep1/crop-200/new_analysis/abierta1_prot.xvg")
ab1_prot <- as.data.frame(sapply(ab1_prot, as.numeric))
colnames(ab1_prot) <- c("time", "x", "y", "z")
ab1_dists <- data.frame("time" = ab1_a3a$time)
ab1_dists$a3a_com <- (sqrt((ab1_a3a$x - ab1_prot$x)^2 + 
                            (ab1_a3a$y - ab1_prot$y)^2 + 
                            (ab1_a3a$z - ab1_prot$z)^2)) - sqrt((ab1_a3a$x[1] - ab1_prot$x[1])^2 + 
                                                                 (ab1_a3a$y[1] - ab1_prot$y[1])^2 + 
                                                                 (ab1_a3a$z[1] - ab1_prot$z[1])^2)


ab1_dists$a3b_com <-( sqrt((ab1_a3b$x - ab1_prot$x)^2 + 
                            (ab1_a3b$y - ab1_prot$y)^2 + 
                            (ab1_a3b$z - ab1_prot$z)^2)) -  sqrt((ab1_a3b$x[1] - ab1_prot$x[1])^2 + 
                                                                   (ab1_a3b$y[1] - ab1_prot$y[1])^2 + 
                                                                   (ab1_a3b$z[1] - ab1_prot$z[1])^2)
ab1_dists$a4a_com <- (sqrt((ab1_a4a$x - ab1_prot$x)^2 + 
                            (ab1_a4a$y - ab1_prot$y)^2 + 
                            (ab1_a4a$z - ab1_prot$z)^2)) - sqrt((ab1_a4a$x[1] - ab1_prot$x[1])^2 + 
                                                                  (ab1_a4a$y[1] - ab1_prot$y[1])^2 + 
                                                                  (ab1_a4a$z[1] - ab1_prot$z[1])^2)
ab1_dists$a4b_com <- (sqrt((ab1_a4b$x - ab1_prot$x)^2 + 
                            (ab1_a4b$y - ab1_prot$y)^2 + 
                            (ab1_a4b$z - ab1_prot$z)^2)) - sqrt((ab1_a4b$x[1] - ab1_prot$x[1])^2 + 
                                                                  (ab1_a4b$y[1] - ab1_prot$y[1])^2 + 
                                                                  (ab1_a4b$z[1] - ab1_prot$z[1])^2)
ab1_dists$a3a_a4a <- (sqrt((ab1_a3a$x - ab1_a4a$x)^2 + 
                            (ab1_a3a$y - ab1_a4a$y)^2 + 
                            (ab1_a3a$z - ab1_a4a$z)^2)) - sqrt((ab1_a3a$x[1] - ab1_a4a$x[1])^2 + 
                                                                 (ab1_a3a$y[1] - ab1_a4a$y[1])^2 + 
                                                                 (ab1_a3a$z[1] - ab1_a4a$z[1])^2)
ab1_dists$a3b_a4b <-( sqrt((ab1_a3b$x - ab1_a4b$x)^2 + 
                            (ab1_a3b$y - ab1_a4b$y)^2 + 
                            (ab1_a3b$z - ab1_a4b$z)^2) ) -  sqrt((ab1_a3b$x[1] - ab1_a4b$x[1])^2 + 
                                                                   (ab1_a3b$y[1] - ab1_a4b$y[1])^2 + 
                                                                   (ab1_a3b$z[1] - ab1_a4b$z[1])^2) 

ab1_fig1 <- ggplot(ab1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") + ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

ab1_fig2 <- ggplot(ab1_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

ab1_fig3 <- ggplot(ab1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/abierta1_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(ab1_fig1, ab1_fig2, ab1_fig3, ncol = 1, nrow = 3)
dev.off()

#### ABIERTA 2 ####
ab2_res190 <- readXVG("abierta-rep2/new_analysis/abierta_res190.xvg")
ab2_res190 <- as.data.frame(sapply(ab2_res190, as.numeric))
colnames(ab2_res190) <- c("time", "distance")
ab2_a3a <- readXVG("abierta-rep2/new_analysis/abierta2_alpha3A.xvg")
ab2_a3a <- as.data.frame(sapply(ab2_a3a, as.numeric))
colnames(ab2_a3a) <- c("time", "x", "y", "z")
ab2_a3b <- readXVG("abierta-rep2/new_analysis/abierta2_alpha3B.xvg")
ab2_a3b <- as.data.frame(sapply(ab2_a3b, as.numeric))
colnames(ab2_a3b) <- c("time", "x", "y", "z")
ab2_a4a <- readXVG("abierta-rep2/new_analysis/abierta2_alpha4A.xvg")
ab2_a4a <- as.data.frame(sapply(ab2_a4a, as.numeric))
colnames(ab2_a4a) <- c("time", "x", "y", "z")
ab2_a4b <- readXVG("abierta-rep2/new_analysis/abierta2_alpha4B.xvg")
ab2_a4b <- as.data.frame(sapply(ab2_a4b, as.numeric))
colnames(ab2_a4b) <- c("time", "x", "y", "z")
ab2_prot <- readXVG("abierta-rep2/new_analysis/abierta2_prot.xvg")
ab2_prot <- as.data.frame(sapply(ab2_prot, as.numeric))
colnames(ab2_prot) <- c("time", "x", "y", "z")
ab2_dists <- data.frame("time" = ab2_a3a$time)
ab2_dists$a3a_com <- (sqrt((ab2_a3a$x - ab2_prot$x)^2 + 
                             (ab2_a3a$y - ab2_prot$y)^2 + 
                             (ab2_a3a$z - ab2_prot$z)^2)) - sqrt((ab2_a3a$x[1] - ab2_prot$x[1])^2 + 
                                                                   (ab2_a3a$y[1] - ab2_prot$y[1])^2 + 
                                                                   (ab2_a3a$z[1] - ab2_prot$z[1])^2)


ab2_dists$a3b_com <-( sqrt((ab2_a3b$x - ab2_prot$x)^2 + 
                             (ab2_a3b$y - ab2_prot$y)^2 + 
                             (ab2_a3b$z - ab2_prot$z)^2)) -  sqrt((ab2_a3b$x[1] - ab2_prot$x[1])^2 + 
                                                                    (ab2_a3b$y[1] - ab2_prot$y[1])^2 + 
                                                                    (ab2_a3b$z[1] - ab2_prot$z[1])^2)
ab2_dists$a4a_com <- (sqrt((ab2_a4a$x - ab2_prot$x)^2 + 
                             (ab2_a4a$y - ab2_prot$y)^2 + 
                             (ab2_a4a$z - ab2_prot$z)^2)) - sqrt((ab2_a4a$x[1] - ab2_prot$x[1])^2 + 
                                                                   (ab2_a4a$y[1] - ab2_prot$y[1])^2 + 
                                                                   (ab2_a4a$z[1] - ab2_prot$z[1])^2)
ab2_dists$a4b_com <- (sqrt((ab2_a4b$x - ab2_prot$x)^2 + 
                             (ab2_a4b$y - ab2_prot$y)^2 + 
                             (ab2_a4b$z - ab2_prot$z)^2)) - sqrt((ab2_a4b$x[1] - ab2_prot$x[1])^2 + 
                                                                   (ab2_a4b$y[1] - ab2_prot$y[1])^2 + 
                                                                   (ab2_a4b$z[1] - ab2_prot$z[1])^2)
ab2_dists$a3a_a4a <- (sqrt((ab2_a3a$x - ab2_a4a$x)^2 + 
                             (ab2_a3a$y - ab2_a4a$y)^2 + 
                             (ab2_a3a$z - ab2_a4a$z)^2)) - sqrt((ab2_a3a$x[1] - ab2_a4a$x[1])^2 + 
                                                                  (ab2_a3a$y[1] - ab2_a4a$y[1])^2 + 
                                                                  (ab2_a3a$z[1] - ab2_a4a$z[1])^2)
ab2_dists$a3b_a4b <-( sqrt((ab2_a3b$x - ab2_a4b$x)^2 + 
                             (ab2_a3b$y - ab2_a4b$y)^2 + 
                             (ab2_a3b$z - ab2_a4b$z)^2) ) -  sqrt((ab2_a3b$x[1] - ab2_a4b$x[1])^2 + 
                                                                    (ab2_a3b$y[1] - ab2_a4b$y[1])^2 + 
                                                                    (ab2_a3b$z[1] - ab2_a4b$z[1])^2) 

ab2_fig1 <- ggplot(ab2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

ab2_fig2 <- ggplot(ab2_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

ab2_fig3 <- ggplot(ab2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "Darkgreen", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/abierta2_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(ab2_fig1, ab2_fig2, ab2_fig3, ncol = 1, nrow = 3)
dev.off()
#### CERRADA 1 ####
ce1_res190 <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_res190.xvg")
ce1_res190 <- as.data.frame(sapply(ce1_res190, as.numeric))
colnames(ce1_res190) <- c("time", "distance")
ce1_a3a <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_alpha3A.xvg")
ce1_a3a <- as.data.frame(sapply(ce1_a3a, as.numeric))
colnames(ce1_a3a) <- c("time", "x", "y", "z")
ce1_a3b <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_alpha3B.xvg")
ce1_a3b <- as.data.frame(sapply(ce1_a3b, as.numeric))
colnames(ce1_a3b) <- c("time", "x", "y", "z")
ce1_a4a <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_alpha4A.xvg")
ce1_a4a <- as.data.frame(sapply(ce1_a4a, as.numeric))
colnames(ce1_a4a) <- c("time", "x", "y", "z")
ce1_a4b <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_alpha4B.xvg")
ce1_a4b <- as.data.frame(sapply(ce1_a4b, as.numeric))
colnames(ce1_a4b) <- c("time", "x", "y", "z")
ce1_prot <- readXVG("cerrada-rep1/prod/new_analysis/cerrada1_prot.xvg")
ce1_prot <- as.data.frame(sapply(ce1_prot, as.numeric))
colnames(ce1_prot) <- c("time", "x", "y", "z")
ce1_dists <- data.frame("time" = ce1_a3a$time)
ce1_dists$a3a_com <- (sqrt((ce1_a3a$x - ce1_prot$x)^2 + 
                             (ce1_a3a$y - ce1_prot$y)^2 + 
                             (ce1_a3a$z - ce1_prot$z)^2)) - sqrt((ce1_a3a$x[1] - ce1_prot$x[1])^2 + 
                                                                   (ce1_a3a$y[1] - ce1_prot$y[1])^2 + 
                                                                   (ce1_a3a$z[1] - ce1_prot$z[1])^2)


ce1_dists$a3b_com <-( sqrt((ce1_a3b$x - ce1_prot$x)^2 + 
                             (ce1_a3b$y - ce1_prot$y)^2 + 
                             (ce1_a3b$z - ce1_prot$z)^2)) -  sqrt((ce1_a3b$x[1] - ce1_prot$x[1])^2 + 
                                                                    (ce1_a3b$y[1] - ce1_prot$y[1])^2 + 
                                                                    (ce1_a3b$z[1] - ce1_prot$z[1])^2)
ce1_dists$a4a_com <- (sqrt((ce1_a4a$x - ce1_prot$x)^2 + 
                             (ce1_a4a$y - ce1_prot$y)^2 + 
                             (ce1_a4a$z - ce1_prot$z)^2)) - sqrt((ce1_a4a$x[1] - ce1_prot$x[1])^2 + 
                                                                   (ce1_a4a$y[1] - ce1_prot$y[1])^2 + 
                                                                   (ce1_a4a$z[1] - ce1_prot$z[1])^2)
ce1_dists$a4b_com <- (sqrt((ce1_a4b$x - ce1_prot$x)^2 + 
                             (ce1_a4b$y - ce1_prot$y)^2 + 
                             (ce1_a4b$z - ce1_prot$z)^2)) - sqrt((ce1_a4b$x[1] - ce1_prot$x[1])^2 + 
                                                                   (ce1_a4b$y[1] - ce1_prot$y[1])^2 + 
                                                                   (ce1_a4b$z[1] - ce1_prot$z[1])^2)
ce1_dists$a3a_a4a <- (sqrt((ce1_a3a$x - ce1_a4a$x)^2 + 
                             (ce1_a3a$y - ce1_a4a$y)^2 + 
                             (ce1_a3a$z - ce1_a4a$z)^2)) - sqrt((ce1_a3a$x[1] - ce1_a4a$x[1])^2 + 
                                                                  (ce1_a3a$y[1] - ce1_a4a$y[1])^2 + 
                                                                  (ce1_a3a$z[1] - ce1_a4a$z[1])^2)
ce1_dists$a3b_a4b <-( sqrt((ce1_a3b$x - ce1_a4b$x)^2 + 
                             (ce1_a3b$y - ce1_a4b$y)^2 + 
                             (ce1_a3b$z - ce1_a4b$z)^2) ) -  sqrt((ce1_a3b$x[1] - ce1_a4b$x[1])^2 + 
                                                                    (ce1_a3b$y[1] - ce1_a4b$y[1])^2 + 
                                                                    (ce1_a3b$z[1] - ce1_a4b$z[1])^2) 

ce1_fig1 <- ggplot(ce1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

ce1_fig2 <- ggplot(ce1_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

ce1_fig3 <- ggplot(ce1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/cerrada1_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(ce1_fig1, ce1_fig2, ce1_fig3, ncol = 1, nrow = 3)
dev.off()

#### CERRADA 2 ####
ce2_res190 <- readXVG("cerrada-rep2/new_analysis/cerrada2_res190.xvg")
ce2_res190 <- as.data.frame(sapply(ce2_res190, as.numeric))
colnames(ce2_res190) <- c("time", "distance")
ce2_a3a <- readXVG("cerrada-rep2/new_analysis/cerrada2-alpha3A.xvg")
ce2_a3a <- as.data.frame(sapply(ce2_a3a, as.numeric))
colnames(ce2_a3a) <- c("time", "x", "y", "z")
ce2_a3b <- readXVG("cerrada-rep2/new_analysis/cerrada2-alpha3B.xvg")
ce2_a3b <- as.data.frame(sapply(ce2_a3b, as.numeric))
colnames(ce2_a3b) <- c("time", "x", "y", "z")
ce2_a4a <- readXVG("cerrada-rep2/new_analysis/cerrada2-alpha4A.xvg")
ce2_a4a <- as.data.frame(sapply(ce2_a4a, as.numeric))
colnames(ce2_a4a) <- c("time", "x", "y", "z")
ce2_a4b <- readXVG("cerrada-rep2/new_analysis/cerrada2-alpha4B.xvg")
ce2_a4b <- as.data.frame(sapply(ce2_a4b, as.numeric))
colnames(ce2_a4b) <- c("time", "x", "y", "z")
ce2_prot <- readXVG("cerrada-rep2/new_analysis/cerrada2-prot.xvg")
ce2_prot <- as.data.frame(sapply(ce2_prot, as.numeric))
colnames(ce2_prot) <- c("time", "x", "y", "z")
ce2_dists <- data.frame("time" = ce2_a3a$time)
ce2_dists$a3a_com <- (sqrt((ce2_a3a$x - ce2_prot$x)^2 + 
                             (ce2_a3a$y - ce2_prot$y)^2 + 
                             (ce2_a3a$z - ce2_prot$z)^2)) - sqrt((ce2_a3a$x[1] - ce2_prot$x[1])^2 + 
                                                                   (ce2_a3a$y[1] - ce2_prot$y[1])^2 + 
                                                                   (ce2_a3a$z[1] - ce2_prot$z[1])^2)


ce2_dists$a3b_com <-( sqrt((ce2_a3b$x - ce2_prot$x)^2 + 
                             (ce2_a3b$y - ce2_prot$y)^2 + 
                             (ce2_a3b$z - ce2_prot$z)^2)) -  sqrt((ce2_a3b$x[1] - ce2_prot$x[1])^2 + 
                                                                    (ce2_a3b$y[1] - ce2_prot$y[1])^2 + 
                                                                    (ce2_a3b$z[1] - ce2_prot$z[1])^2)
ce2_dists$a4a_com <- (sqrt((ce2_a4a$x - ce2_prot$x)^2 + 
                             (ce2_a4a$y - ce2_prot$y)^2 + 
                             (ce2_a4a$z - ce2_prot$z)^2)) - sqrt((ce2_a4a$x[1] - ce2_prot$x[1])^2 + 
                                                                   (ce2_a4a$y[1] - ce2_prot$y[1])^2 + 
                                                                   (ce2_a4a$z[1] - ce2_prot$z[1])^2)
ce2_dists$a4b_com <- (sqrt((ce2_a4b$x - ce2_prot$x)^2 + 
                             (ce2_a4b$y - ce2_prot$y)^2 + 
                             (ce2_a4b$z - ce2_prot$z)^2)) - sqrt((ce2_a4b$x[1] - ce2_prot$x[1])^2 + 
                                                                   (ce2_a4b$y[1] - ce2_prot$y[1])^2 + 
                                                                   (ce2_a4b$z[1] - ce2_prot$z[1])^2)
ce2_dists$a3a_a4a <- (sqrt((ce2_a3a$x - ce2_a4a$x)^2 + 
                             (ce2_a3a$y - ce2_a4a$y)^2 + 
                             (ce2_a3a$z - ce2_a4a$z)^2)) - sqrt((ce2_a3a$x[1] - ce2_a4a$x[1])^2 + 
                                                                  (ce2_a3a$y[1] - ce2_a4a$y[1])^2 + 
                                                                  (ce2_a3a$z[1] - ce2_a4a$z[1])^2)
ce2_dists$a3b_a4b <-( sqrt((ce2_a3b$x - ce2_a4b$x)^2 + 
                             (ce2_a3b$y - ce2_a4b$y)^2 + 
                             (ce2_a3b$z - ce2_a4b$z)^2) ) -  sqrt((ce2_a3b$x[1] - ce2_a4b$x[1])^2 + 
                                                                    (ce2_a3b$y[1] - ce2_a4b$y[1])^2 + 
                                                                    (ce2_a3b$z[1] - ce2_a4b$z[1])^2) 

ce2_fig1 <- ggplot(ce2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

ce2_fig2 <- ggplot(ce2_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

ce2_fig3 <- ggplot(ce2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "dodgerblue3", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/cerrada2_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(ce2_fig1, ce2_fig2, ce2_fig3, ncol = 1, nrow = 3)
dev.off()

#### MUTANTE 1 ####
mu1_res190 <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_res190.xvg")
mu1_res190 <- as.data.frame(sapply(mu1_res190, as.numeric))
colnames(mu1_res190) <- c("time", "distance")
mu1_a3a <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_dist_alpha3A.xvg")
mu1_a3a <- as.data.frame(sapply(mu1_a3a, as.numeric))
colnames(mu1_a3a) <- c("time", "x", "y", "z")
mu1_a3b <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_dist_alpha3B.xvg")
mu1_a3b <- as.data.frame(sapply(mu1_a3b, as.numeric))
colnames(mu1_a3b) <- c("time", "x", "y", "z")
mu1_a4a <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_dist_alpha4A.xvg")
mu1_a4a <- as.data.frame(sapply(mu1_a4a, as.numeric))
colnames(mu1_a4a) <- c("time", "x", "y", "z")
mu1_a4b <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_dist_alpha4B.xvg")
mu1_a4b <- as.data.frame(sapply(mu1_a4b, as.numeric))
colnames(mu1_a4b) <- c("time", "x", "y", "z")
mu1_prot <- readXVG("v160m-rep1/prod_200ns/new_analysis/v160m1_dist_prot.xvg")
mu1_prot <- as.data.frame(sapply(mu1_prot, as.numeric))
colnames(mu1_prot) <- c("time", "x", "y", "z")
mu1_dists <- data.frame("time" = mu1_a3a$time)
mu1_dists$a3a_com <- (sqrt((mu1_a3a$x - mu1_prot$x)^2 + 
                             (mu1_a3a$y - mu1_prot$y)^2 + 
                             (mu1_a3a$z - mu1_prot$z)^2)) - sqrt((mu1_a3a$x[1] - mu1_prot$x[1])^2 + 
                                                                   (mu1_a3a$y[1] - mu1_prot$y[1])^2 + 
                                                                   (mu1_a3a$z[1] - mu1_prot$z[1])^2)


mu1_dists$a3b_com <-( sqrt((mu1_a3b$x - mu1_prot$x)^2 + 
                             (mu1_a3b$y - mu1_prot$y)^2 + 
                             (mu1_a3b$z - mu1_prot$z)^2)) -  sqrt((mu1_a3b$x[1] - mu1_prot$x[1])^2 + 
                                                                    (mu1_a3b$y[1] - mu1_prot$y[1])^2 + 
                                                                    (mu1_a3b$z[1] - mu1_prot$z[1])^2)
mu1_dists$a4a_com <- (sqrt((mu1_a4a$x - mu1_prot$x)^2 + 
                             (mu1_a4a$y - mu1_prot$y)^2 + 
                             (mu1_a4a$z - mu1_prot$z)^2)) - sqrt((mu1_a4a$x[1] - mu1_prot$x[1])^2 + 
                                                                   (mu1_a4a$y[1] - mu1_prot$y[1])^2 + 
                                                                   (mu1_a4a$z[1] - mu1_prot$z[1])^2)
mu1_dists$a4b_com <- (sqrt((mu1_a4b$x - mu1_prot$x)^2 + 
                             (mu1_a4b$y - mu1_prot$y)^2 + 
                             (mu1_a4b$z - mu1_prot$z)^2)) - sqrt((mu1_a4b$x[1] - mu1_prot$x[1])^2 + 
                                                                   (mu1_a4b$y[1] - mu1_prot$y[1])^2 + 
                                                                   (mu1_a4b$z[1] - mu1_prot$z[1])^2)
mu1_dists$a3a_a4a <- (sqrt((mu1_a3a$x - mu1_a4a$x)^2 + 
                             (mu1_a3a$y - mu1_a4a$y)^2 + 
                             (mu1_a3a$z - mu1_a4a$z)^2)) - sqrt((mu1_a3a$x[1] - mu1_a4a$x[1])^2 + 
                                                                  (mu1_a3a$y[1] - mu1_a4a$y[1])^2 + 
                                                                  (mu1_a3a$z[1] - mu1_a4a$z[1])^2)
mu1_dists$a3b_a4b <-( sqrt((mu1_a3b$x - mu1_a4b$x)^2 + 
                             (mu1_a3b$y - mu1_a4b$y)^2 + 
                             (mu1_a3b$z - mu1_a4b$z)^2) ) -  sqrt((mu1_a3b$x[1] - mu1_a4b$x[1])^2 + 
                                                                    (mu1_a3b$y[1] - mu1_a4b$y[1])^2 + 
                                                                    (mu1_a3b$z[1] - mu1_a4b$z[1])^2) 

mu1_fig1 <- ggplot(mu1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "purple", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

mu1_fig2 <- ggplot(mu1_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "purple", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

mu1_fig3 <- ggplot(mu1_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "purple", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/mutante1_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(mu1_fig1, mu1_fig2, mu1_fig3, ncol = 1, nrow = 3)
dev.off()

#### MUTANTE 2 ####
mu2_res190 <- readXVG("v160m-rep2/new_analysis/v160m2_res190.xvg")
mu2_res190 <- as.data.frame(sapply(mu2_res190, as.numeric))
colnames(mu2_res190) <- c("time", "distance")
mu2_a3a <- readXVG("v160m-rep2/new_analysis/v160m2_dist_alpha3A.xvg")
mu2_a3a <- as.data.frame(sapply(mu2_a3a, as.numeric))
colnames(mu2_a3a) <- c("time", "x", "y", "z")
mu2_a3b <- readXVG("v160m-rep2/new_analysis/v160m2_dist_alpha3B.xvg")
mu2_a3b <- as.data.frame(sapply(mu2_a3b, as.numeric))
colnames(mu2_a3b) <- c("time", "x", "y", "z")
mu2_a4a <- readXVG("v160m-rep2/new_analysis/v160m2_dist_alpha4A.xvg")
mu2_a4a <- as.data.frame(sapply(mu2_a4a, as.numeric))
colnames(mu2_a4a) <- c("time", "x", "y", "z")
mu2_a4b <- readXVG("v160m-rep2/new_analysis/v160m2_dist_alpha4B.xvg")
mu2_a4b <- as.data.frame(sapply(mu2_a4b, as.numeric))
colnames(mu2_a4b) <- c("time", "x", "y", "z")
mu2_prot <- readXVG("v160m-rep2/new_analysis/v160m2_dist_prot.xvg")
mu2_prot <- as.data.frame(sapply(mu2_prot, as.numeric))
colnames(mu2_prot) <- c("time", "x", "y", "z")
mu2_dists <- data.frame("time" = mu2_a3a$time)
mu2_dists$a3a_com <- (sqrt((mu2_a3a$x - mu2_prot$x)^2 + 
                             (mu2_a3a$y - mu2_prot$y)^2 + 
                             (mu2_a3a$z - mu2_prot$z)^2)) - sqrt((mu2_a3a$x[1] - mu2_prot$x[1])^2 + 
                                                                   (mu2_a3a$y[1] - mu2_prot$y[1])^2 + 
                                                                   (mu2_a3a$z[1] - mu2_prot$z[1])^2)


mu2_dists$a3b_com <-( sqrt((mu2_a3b$x - mu2_prot$x)^2 + 
                             (mu2_a3b$y - mu2_prot$y)^2 + 
                             (mu2_a3b$z - mu2_prot$z)^2)) -  sqrt((mu2_a3b$x[1] - mu2_prot$x[1])^2 + 
                                                                    (mu2_a3b$y[1] - mu2_prot$y[1])^2 + 
                                                                    (mu2_a3b$z[1] - mu2_prot$z[1])^2)
mu2_dists$a4a_com <- (sqrt((mu2_a4a$x - mu2_prot$x)^2 + 
                             (mu2_a4a$y - mu2_prot$y)^2 + 
                             (mu2_a4a$z - mu2_prot$z)^2)) - sqrt((mu2_a4a$x[1] - mu2_prot$x[1])^2 + 
                                                                   (mu2_a4a$y[1] - mu2_prot$y[1])^2 + 
                                                                   (mu2_a4a$z[1] - mu2_prot$z[1])^2)
mu2_dists$a4b_com <- (sqrt((mu2_a4b$x - mu2_prot$x)^2 + 
                             (mu2_a4b$y - mu2_prot$y)^2 + 
                             (mu2_a4b$z - mu2_prot$z)^2)) - sqrt((mu2_a4b$x[1] - mu2_prot$x[1])^2 + 
                                                                   (mu2_a4b$y[1] - mu2_prot$y[1])^2 + 
                                                                   (mu2_a4b$z[1] - mu2_prot$z[1])^2)
mu2_dists$a3a_a4a <- (sqrt((mu2_a3a$x - mu2_a4a$x)^2 + 
                             (mu2_a3a$y - mu2_a4a$y)^2 + 
                             (mu2_a3a$z - mu2_a4a$z)^2)) - sqrt((mu2_a3a$x[1] - mu2_a4a$x[1])^2 + 
                                                                  (mu2_a3a$y[1] - mu2_a4a$y[1])^2 + 
                                                                  (mu2_a3a$z[1] - mu2_a4a$z[1])^2)
mu2_dists$a3b_a4b <-( sqrt((mu2_a3b$x - mu2_a4b$x)^2 + 
                             (mu2_a3b$y - mu2_a4b$y)^2 + 
                             (mu2_a3b$z - mu2_a4b$z)^2) ) -  sqrt((mu2_a3b$x[1] - mu2_a4b$x[1])^2 + 
                                                                    (mu2_a3b$y[1] - mu2_a4b$y[1])^2 + 
                                                                    (mu2_a3b$z[1] - mu2_a4b$z[1])^2) 

mu2_fig1 <- ggplot(mu2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_com), color = "purple", alpha = 1) + 
  geom_line(aes(y = a3b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 3 y el COM")

mu2_fig2 <- ggplot(mu2_dists, aes(x = time)) + 
  geom_line(aes(y = a4a_com), color = "purple", alpha = 1) + 
  geom_line(aes(y = a4b_com), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia de las helices alpha 4 y el COM")

mu2_fig3 <- ggplot(mu2_dists, aes(x = time)) + 
  geom_line(aes(y = a3a_a4a), color = "purple", alpha = 1) + 
  geom_line(aes(y = a3b_a4b), color = "gray60", alpha = 1) + a + 
  ylab("Distance (nm)") + xlab("Time (ns)") +  ylim(-0.3, 0.3) +
  ggtitle("Distancia entre helices alpha 3 y 4")

png("new_analysis/mutante2_poldomain_distances.png", height = 200, width = 184, units = "mm", res = 300)
ggarrange(mu2_fig1, mu2_fig2, mu2_fig3, ncol = 1, nrow = 3)
dev.off()

#### RESIDUO 190 ####

res190 <- data.frame("time" = c(ab1_res190$time, ab2_res190$time,
                                ce1_res190$time, ce2_res190$time,
                                mu1_res190$time, mu2_res190$time),
                     "dist" = c(ab1_res190$distance, ab2_res190$distance,
                                ce1_res190$distance, ce2_res190$distance,
                                mu1_res190$distance, mu2_res190$distance), 
                     "Structure" = "s")

res190$Structure[1:1001] <- "Open STING 1"
res190$Structure[1002:2002] <- "Open STING 2"
res190$Structure[2003:3003] <- "Closed STING 1"
res190$Structure[3004:4004] <- "Closed STING 2"
res190$Structure[4005:5005] <- "V160M STING 1"
res190$Structure[5006:6006] <- "V160M STING 2"

res190$Structure <- as.factor(res190$Structure)
res190$Structure <- ordered(res190$Structure, levels=c("Open STING 1", "Open STING 2",
                                                     "Closed STING 1","Closed STING 2",
                                                     "V160M STING 1", "V160M STING 2"))

png("new_analysis/residue190_distances.png", height = 90, width = 184, units = "mm", res = 300)
ggplot(res190, aes(x = time, y = dist, color = Structure)) +
  geom_line(size = .4, alpha = .6) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("palegreen", "Darkgreen", 
                                "dodgerblue", "blue4",
                                "magenta", "Purple")) + a +
  ggtitle("Distancia de los residuos THR190 (HIS185 en hSTING)") +
  geom_point(aes(x = 192.6, y = 4.191), color = "palegreen", size = 3) + 
  geom_point(aes(x = 132.4, y = 4.938), color = "Darkgreen", size = 3) + 
  geom_point(aes(x = 179.4, y = 3.424), color = "dodgerblue", size = 3) + 
  geom_point(aes(x = 169.4, y = 3.54), color = "blue4", size = 3) + 
  geom_point(aes(x = 108.6, y = 5.161), color = "magenta", size = 3) + 
  geom_point(aes(x = 146.2, y = 5.014), color = "Purple", size = 3) 
dev.off()

  
#### NEW ALT

ab2_polA <- readXVG("abierta-rep2/new_analysis/abierta2_pol_chainsA.xvg")
ab2_polA <- as.data.frame(sapply(ab2_polA, as.numeric))
colnames(ab2_polA) <- c("time", "x", "y", "z")
ab2_polB <- readXVG("abierta-rep2/new_analysis/abierta2_pol_chainsB.xvg")
ab2_polB <- as.data.frame(sapply(ab2_polB, as.numeric))
colnames(ab2_polB) <- c("time", "x", "y", "z")
ab2_prot <- readXVG("abierta-rep2/new_analysis/abierta2_prot.xvg")
ab2_prot <- as.data.frame(sapply(ab2_prot, as.numeric))
colnames(ab2_prot) <- c("time", "x", "y", "z")

ce2_polA <- readXVG("cerrada-rep2/new_analysis/cerrada2_pol_chainA.xvg")
ce2_polA <- as.data.frame(sapply(ce2_polA, as.numeric))
colnames(ce2_polA) <- c("time", "x", "y", "z")
ce2_polB <- readXVG("cerrada-rep2/new_analysis/cerrada2_pol_chainB.xvg")
ce2_polB <- as.data.frame(sapply(ce2_polB, as.numeric))
colnames(ce2_polB) <- c("time", "x", "y", "z")
ce2_prot <- readXVG("cerrada-rep2/new_analysis/cerrada2-prot.xvg")
ce2_prot <- as.data.frame(sapply(ce2_prot, as.numeric))
colnames(ce2_prot) <- c("time", "x", "y", "z")

dists <- data.frame("time" = ce2_polA$time)

dists$ab_polA_com <- (sqrt((ab2_polA$x - ab2_prot$x)^2 + 
                             (ab2_polA$y - ab2_prot$y)^2 + 
                             (ab2_polA$z - ab2_prot$z)^2)) - sqrt((ab2_polA$x[1] - ab2_prot$x[1])^2 + 
                                                                   (ab2_polA$y[1] - ab2_prot$y[1])^2 + 
                                                                   (ab2_polA$z[1] - ab2_prot$z[1])^2)
dists$ab_polB_com <- (sqrt((ab2_polB$x - ab2_prot$x)^2 + 
                             (ab2_polB$y - ab2_prot$y)^2 + 
                             (ab2_polB$z - ab2_prot$z)^2)) - sqrt((ab2_polB$x[1] - ab2_prot$x[1])^2 + 
                                                                   (ab2_polB$y[1] - ab2_prot$y[1])^2 + 
                                                                   (ab2_polB$z[1] - ab2_prot$z[1])^2)

dists$ce_polA_com <- (sqrt((ce2_polA$x - ab2_prot$x)^2 + 
                             (ce2_polA$y - ab2_prot$y)^2 + 
                             (ce2_polA$z - ab2_prot$z)^2)) - sqrt((ce2_polA$x[1] - ab2_prot$x[1])^2 + 
                                                                    (ce2_polA$y[1] - ab2_prot$y[1])^2 + 
                                                                    (ce2_polA$z[1] - ab2_prot$z[1])^2)
dists$ce_polB_com <- (sqrt((ce2_polB$x - ce2_prot$x)^2 + 
                             (ce2_polB$y - ce2_prot$y)^2 + 
                             (ce2_polB$z - ce2_prot$z)^2)) - sqrt((ce2_polB$x[1] - ce2_prot$x[1])^2 + 
                                                                    (ce2_polB$y[1] - ce2_prot$y[1])^2 + 
                                                                    (ce2_polB$z[1] - ce2_prot$z[1])^2)


dists2 <- data.frame("time" = c(dists$time, dists$time, dists$time, dists$time),
                     "dist" = c(dists$ab_polA_com, dists$ab_polB_com, dists$ce_polA_com, dists$ce_polB_com),
                     "Structure" = "Structure") 

dists2$Structure[1:1001] <- "Open STING - Chain A"
dists2$Structure[1002:2002] <- "Open STING - Chain B"
dists2$Structure[2003:3003] <- "Closed STING - Chain A"
dists2$Structure[3004:4004] <- "Closed STING - Chain B"

dists2$Structure <- as.factor(dists2$Structure)
dists2$Structure <- ordered(dists2$Structure, levels=c("Open STING - Chain A",
                                                     "Open STING - Chain B",
                                                     "Closed STING - Chain A", 
                                                     "Closed STING - Chain B"))
dist_both <- ggplot(dists2, aes(x = time, y = dist, color = Structure)) + 
  scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) + 
  geom_hline(yintercept = 0) +
  #geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank()) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)")  +  a + ylim(-0.5, 0.5) 

figd <- ggplot(dists2[c(501:1001, 1502:2002, 2503:3003,3504:4004),], 
               aes(dist, fill = Structure, color = Structure)) + coord_flip() +
  geom_density(alpha = .7) + a + ylab("Density") + theme(axis.title.y = element_blank()) + 
  scale_fill_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) +
  theme(axis.title = element_blank()) + 
  scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) + xlim(-0.5, 0.5)

png("figuras_paper/figura_4/distancias_oswald_alt_with_chains.png", height = 60, width = 184, units = "mm", res = 300)

ggarrange(dist_both, figd, ncol = 2, nrow = 1, widths = c(2,1),
          common.legend = TRUE, legend = "bottom")

dev.off()
