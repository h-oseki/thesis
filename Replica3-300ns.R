rmsd3 <- read.csv("cerrada-rep3/prod-300/cerrada-rep3-300ns-rmsd.xvg", sep="")
gr3 <- read.csv("cerrada-rep3/prod-300/cerrada-rep3-300ns-rg.xvg", sep="")
gr3$time <- gr3$time/1000
bfac3 <- read.pdb("cerrada-rep3/prod-300/cerrada-rep3-300ns-bfac.pdb")
bfac3_db <- data.frame("NUM" = bfac3$atom$resno, 
                       "RES" = bfac3$atom$resid, 
                       "BFAC" = bfac3$atom$b)
prom3 <- mean(bfac3_db$BFAC)
sd3 <- sd(bfac3_db$BFAC)
bfac3_db$NBFAC <- (bfac3_db$BFAC - prom3)/sd3
j = 146
for (i in 1:394){
  bfac3_db$ATOM[i] = j
  j = j+1
}
rm(i, j, prom3, sd3)
png("cerrada-rep3/prod-300/rep3-300ns-rmsd.png", height = 200, width = 300, units = "mm", res = 300)
rep3_ext_rmsd <- ggplot(rmsd3, aes(x = time, y = rmsd)) + 
  geom_line(color = "deeppink1", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + 
  ggtitle("Root mean square deviation")
dev.off()
png("cerrada-rep3/prod-300/rep3-300ns-rg.png", height = 200, width = 300, units = "mm", res = 300)
rep3_ext_gr <- ggplot(gr3, aes(x = time, y = all)) + 
  geom_line(color = "deeppink1", size = 1) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  a + ggtitle("Radius of gyration")
dev.off()


bfacs_rep3_ext_A <- ggplot(bfac3_db[1:197,], aes(x = ATOM, y = NBFAC)) + ylim(-0.6, 9) + 
  geom_line(color = "deeppink3", size = 1) + ylab("Normalized B-factor") + a +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.6, ymax=9, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_rep3_ext_B <- ggplot(bfac3_db[198:394,], aes(x = ATOM, y = NBFAC)) + 
  geom_line(color = "deeppink3", size = 1) + ylab("Normalized B-factor") + a +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.6, ymax=9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.6, ymax=9, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 
png("cerrada-rep3/prod-300/rep3-300ns-bfacs.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(bfacs_reps_A, bfacs_reps_B,  labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom") 
dev.off()


rep3_ext_chainA <- read.delim("cerrada-rep3/prod-300/distance/chainA-rep3-ext.xvg")
rep3_ext_chainB <- read.delim("cerrada-rep3/prod-300/distance/chainB-rep3-ext.xvg")
rep3_ext_protein <- read.delim("cerrada-rep3/prod-300/distance/protein-rep3-ext.xvg")

rep3_ext <- data.frame("time" = rep3_ext_chainA$time/1000, 
                       "chainA_x" = rep3_ext_chainA$x, "chainA_y" = rep3_ext_chainA$y, "chainA_z" = rep3_ext_chainA$z,
                       "chainB_x" = rep3_ext_chainB$x, "chainB_y" = rep3_ext_chainB$y, "chainB_z" = rep3_ext_chainB$z, 
                       "protein_x" = rep3_ext_protein$x,"protein_y" = rep3_ext_protein$y, "protein_z" = rep3_ext_protein$z)
rep3_ext$chainA_dist <- sqrt((rep3_ext$chainA_x - rep3_ext$protein_x)^2 + 
                               (rep3_ext$chainA_y - rep3_ext$protein_y)^2 + 
                               (rep3_ext$chainA_z - rep3_ext$protein_z)^2)
rep3_ext$chainB_dist <- sqrt((rep3_ext$chainB_x - rep3_ext$protein_x)^2 + 
                               (rep3_ext$chainB_y - rep3_ext$protein_y)^2 + 
                               (rep3_ext$chainB_z - rep3_ext$protein_z)^2)  
rep3_ext$dist_normA <- rep3_ext$chainA_dist - rep3_ext$chainA_dist[1]
rep3_ext$dist_normB <- rep3_ext$chainB_dist - rep3_ext$chainB_dist[1]

distances_rep3_ext <- data.frame("Time" = c(rep3_ext_chainA$time, rep3_ext_chainA$time),
                                 "Distance" = c(rep3_ext$dist_normA, rep3_ext$dist_normB),
                                 "Label" = "label")
rm(rep3_ext_chainA, rep3_ext_chainB, rep3_ext_protein, rep3_ext)

distances_rep3_ext$Label[1:1501] <- "Rep3 ext - Chain A" 
distances_rep3_ext$Label[1502:3002] <- "Rep3 ext - Chain B" 

rep3_ext_dist <- ggplot(distances_rep3_ext[1:3002,], aes(x = Time/1000, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("royalblue", "orangered")) + a + ylim(-0.25, 0.25) +
  ggtitle("Replica 3 extendida") 

