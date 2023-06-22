rmsd_33cgamp <- read.csv("cerrada+33cGAMP/prod_200/cerrada+33cgamp-200ns-rmsd.xvg", sep="")
gr_33cgamp <- read.csv("cerrada+33cGAMP/prod_200/cerrada+33cgamp-200ns-rg.xvg", sep="")
gr_33cgamp$time <- gr_33cgamp$time/1000
bfac_33cgamp <- read.pdb("cerrada+33cGAMP/prod_200/cerrada+33cgamp-200ns-bfac.pdb")
bfac_33cgamp_db <- data.frame("NUM" = bfac_33cgamp$atom$resno, 
                       "RES" = bfac_33cgamp$atom$resid, 
                       "BFAC" = bfac_33cgamp$atom$b)
prom_33cgamp <- mean(bfac_33cgamp_db$BFAC)
sd_33cgamp <- sd(bfac_33cgamp_db$BFAC)
bfac_33cgamp_db$NBFAC <- (bfac_33cgamp_db$BFAC - prom_33cgamp)/sd_33cgamp
j = 146
for (i in 1:394){
  bfac_33cgamp_db$ATOM[i] = j
  j = j+1
}
rm(i, j, prom_33cgamp, sd_33cgamp)
png("cerrada+33cGAMP/prod_200/holo+33cgamp-rmsd.png", height = 200, width = 300, units = "mm", res = 300)
rmsd_33cgamp <- ggplot(rmsd_33cgamp, aes(x = time, y = rmsd)) + 
  geom_line(color = "turquoise4", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
dev.off()
png("cerrada+33cGAMP/prod_200/holo+33cgamp-rg.png", height = 200, width = 300, units = "mm", res = 300)
gr_33cgamp <- ggplot(gr_33cgamp, aes(x = time, y = all)) + 
  geom_line(color = "turquoise4", size = 1) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  a + ggtitle("Radius of gyration")
dev.off()
bfacs_33cgamp_A <- ggplot(bfac_33cgamp_db[1:197,], aes(x = ATOM, y = NBFAC))  + 
  geom_line(color = "turquoise4", size = 1) + ylab("Normalized B-factor") + a +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.6, ymax=7.9, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_33cgamp_B <- ggplot(bfac_33cgamp_db[198:394,], aes(x = ATOM, y = NBFAC)) + 
  geom_line(color = "turquoise4", size = 1) + ylab("Normalized B-factor") + a +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.6, ymax=7.9, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.6, ymax=7.9, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 
png("cerrada+33cGAMP/prod_200/holo+33cgamp-rmsf.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(bfacs_33cgamp_A, bfacs_33cgamp_B,
          common.legend = TRUE, legend = "bottom") 
dev.off()

rep1_33cgamp_chainA <- read.delim("cerrada+33cGAMP/prod_200/distances/chainA-33cgamp-rep1.xvg")
rep1_33cgamp_chainB <- read.delim("cerrada+33cGAMP/prod_200/distances/chainB-33cgamp-rep1.xvg")
rep1_33cgamp_protein <- read.delim("cerrada+33cGAMP/prod_200/distances/protein-33cgamp-rep1.xvg")

rep1_33cgamp <- data.frame("time" = rep1_33cgamp_chainA$time/1000, 
                       "chainA_x" = rep1_33cgamp_chainA$x, "chainA_y" = rep1_33cgamp_chainA$y, "chainA_z" = rep1_33cgamp_chainA$z,
                       "chainB_x" = rep1_33cgamp_chainB$x, "chainB_y" = rep1_33cgamp_chainB$y, "chainB_z" = rep1_33cgamp_chainB$z, 
                       "protein_x" = rep1_33cgamp_protein$x,"protein_y" = rep1_33cgamp_protein$y, "protein_z" = rep1_33cgamp_protein$z)
rep1_33cgamp$chainA_dist <- sqrt((rep1_33cgamp$chainA_x - rep1_33cgamp$protein_x)^2 + 
                               (rep1_33cgamp$chainA_y - rep1_33cgamp$protein_y)^2 + 
                               (rep1_33cgamp$chainA_z - rep1_33cgamp$protein_z)^2)
rep1_33cgamp$chainB_dist <- sqrt((rep1_33cgamp$chainB_x - rep1_33cgamp$protein_x)^2 + 
                               (rep1_33cgamp$chainB_y - rep1_33cgamp$protein_y)^2 + 
                               (rep1_33cgamp$chainB_z - rep1_33cgamp$protein_z)^2)  
rep1_33cgamp$dist_normA <- rep1_33cgamp$chainA_dist - rep1_33cgamp$chainA_dist[1]
rep1_33cgamp$dist_normB <- rep1_33cgamp$chainB_dist - rep1_33cgamp$chainB_dist[1]

distances_33cgamp <- data.frame("Time" = c(rep1_33cgamp_chainA$time, rep1_33cgamp_chainA$time),
                         "Distance" = c(rep1_33cgamp$dist_normA, rep1_33cgamp$dist_normB),
                         "Label" = "label")
rm(rep1_33cgamp_chainA, rep1_33cgamp_chainB, rep1_33cgamp_protein, rep1_33cgamp)

distances_33cgamp$Label[1:1001] <- "33cGAMP Rep1 - Chain A" 
distances_33cgamp$Label[1002:2002] <- "33cGAMP Rep1 - Chain B" 

rep1_33cgamp_dist <- ggplot(distances_33cgamp[1:2002,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("darkslateblue", "olivedrab3")) + a + ylim(-0.25, 0.25) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of STING") 
