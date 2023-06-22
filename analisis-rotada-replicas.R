setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

rmsd_rot <- read.csv("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-rmsd.xvg", sep="")
rmsd_rot2 <- read.csv("cerrada-ligrot-rep2/cerrada-ligrot-rep2-rmsd.xvg", sep="")
rmsd_rots <- data.frame("rmsd" = c(rmsd_rot$rmsd, rmsd_rot2$rmsd),
                   "time" = c(rmsd_rot$time, rmsd_rot2$time),
                   "Structure" = "label")
rmsd_rots$Structure[1:1001] <- "replica 1"
rmsd_rots$Structure[1002:2002] <- "replica 2"
rm(rmsd_rot, rmsd_rot2)

gr_rot <- read.csv("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-gr.xvg", sep="")
gr_rot2 <- read.csv("cerrada-ligrot-rep2/cerrada-ligrot-rep2-gr.xvg", sep="")
gr_rots <- data.frame("time" = c(gr_rot$time, gr_rot2$time),
                 "all" = c(gr_rot$all, gr_rot2$all),
                 "x" = c(gr_rot$x, gr_rot2$x),
                 "y" = c(gr_rot$y, gr_rot2$y),
                 "z" = c(gr_rot$z, gr_rot2$z),
                 "Structure" = "label")

gr_rots$time <- gr_rots$time/1000
gr_rots$Structure[1:1001] <- "replica 1"
gr_rots$Structure[1002:2002] <- "replica 2"
rm(gr_rot, gr_rot2)

bfac_rot <- read.pdb("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-bfac.pdb")
bfac_rot2 <- read.pdb("cerrada-ligrot-rep2/cerrada-ligrot-rep2-bfac.pdb")
bfac_db_rot <- data.frame("NUM" = bfac_rot$atom$resno, 
                       "RES" = bfac_rot$atom$resid, 
                       "BFAC" = bfac_rot$atom$b)
prom_rot <- mean(bfac_db_rot$BFAC)
sd_rot <- sd(bfac_db_rot$BFAC)
bfac_db_rot$NBFAC <- (bfac_db_rot$BFAC - prom_rot)/sd_rot
j = 146
for (i in 1:394){
  bfac_db_rot$ATOM[i] = j
  j = j+1
}

bfac_db_rot2 <- data.frame("NUM" = bfac_rot2$atom$resno, 
                          "RES" = bfac_rot2$atom$resid, 
                          "BFAC" = bfac_rot2$atom$b)
prom_rot2 <- mean(bfac_db_rot2$BFAC)
sd_rot2 <- sd(bfac_db_rot2$BFAC)
bfac_db_rot2$NBFAC <- (bfac_db_rot2$BFAC - prom_rot2)/sd_rot2
j = 146
for (i in 1:394){
  bfac_db_rot2$ATOM[i] = j
  j = j+1
}

bfacs_rots <- data.frame("res" = c(bfac_db_rot$ATOM, bfac_db_rot2$ATOM),
                    "nbfac" = c(bfac_db_rot$NBFAC, bfac_db_rot2$NBFAC),
                    "Structure" = "label")
bfacs_rots$Structure[1:394] <- "replica 1"
bfacs_rots$Structure[395:788] <- "replica 2"
rm(i, j, prom_rot, prom_rot2, sd_rot, sd_rot2, bfac_rot, bfac_db_rot, bfac_rot2, bfac_db_rot2)


rot_rmsd <- ggplot(rmsd_rots, aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("dodgerblue3", "deeppink")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("ligrot-replicas-rmsds.png", height = 150, width = 300, units = "mm", res = 300)
rot_rmsd
dev.off()


rot_gr <- ggplot(gr_rots, aes(x = time, y = all)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("dodgerblue3", "deeppink")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  scale_linetype_manual(values=c("dashed", "solid")) + a + ggtitle("Radius of gyration")
png("ligrot-replicas-grs.png", height = 150, width = 300, units = "mm", res = 300)
rot_gr
dev.off()

rot_bfacA <- ggplot(bfacs_rots[c(1:197, 395:591),], aes(x = res, y = nbfac)) + ylim(-0.6, 9) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("dodgerblue3", "deeppink")) +
  scale_linetype_manual(values=c("twodash", "solid")) + a +
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
rot_bfacB <- ggplot(bfacs_rots[c(198:394, 592:788),], aes(x = res, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("dodgerblue3", "deeppink")) +
  scale_linetype_manual(values=c("twodash", "solid")) + a + ylim(-0.6, 9) +
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
png("ligrot-replicas-bfacs.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(rot_bfacA, rot_bfacB,  labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom") 
dev.off()

rot_rep1_chainA <- read.delim("cerrada-ligrot-rep1/prod/distances/chainA-rot-rep1.xvg")
rot_rep1_chainB <- read.delim("cerrada-ligrot-rep1/prod/distances/chainB-rot-rep1.xvg")
rot_rep1_protein <- read.delim("cerrada-ligrot-rep1/prod/distances/protein-rot-rep1.xvg")
rot_rep2_chainA <- read.delim("cerrada-ligrot-rep2/distance/chainA-rot-rep2.xvg")
rot_rep2_chainB <- read.delim("cerrada-ligrot-rep2/distance/chainB-rot-rep2.xvg")
rot_rep2_protein <- read.delim("cerrada-ligrot-rep2/distance/protein-rot-rep2.xvg")

rot_rep1_hbondA <- read.table("cerrada-ligrot-rep1/prod/distances/ser268a-lig.dat", quote="\"", comment.char="")
rot_rep1_hbondA <- rot_rep1_hbondA[-1,]
rot_rep1_hbondB <- read.table("cerrada-ligrot-rep1/prod/distances/ser268b-lig.dat", quote="\"", comment.char="")
rot_rep1_hbondB <- rot_rep1_hbondB[-1,]

rot_rep1 <- data.frame("time" = rot_rep1_chainA$time/1000, 
                   "chainA_x" = rot_rep1_chainA$x, "chainA_y" = rot_rep1_chainA$y, "chainA_z" = rot_rep1_chainA$z,
                   "chainB_x" = rot_rep1_chainB$x, "chainB_y" = rot_rep1_chainB$y, "chainB_z" = rot_rep1_chainB$z, 
                   "protein_x" = rot_rep1_protein$x,"protein_y" = rot_rep1_protein$y, "protein_z" = rot_rep1_protein$z,
                   "hbondA" = rot_rep1_hbondA$V2, "hbondB" = rot_rep1_hbondB$V2)
rot_rep1$chainA_dist <- sqrt((rot_rep1$chainA_x - rot_rep1$protein_x)^2 + 
                           (rot_rep1$chainA_y - rot_rep1$protein_y)^2 + 
                           (rot_rep1$chainA_z - rot_rep1$protein_z)^2)
rot_rep1$chainB_dist <- sqrt((rot_rep1$chainB_x - rot_rep1$protein_x)^2 + 
                           (rot_rep1$chainB_y - rot_rep1$protein_y)^2 + 
                           (rot_rep1$chainB_z - rot_rep1$protein_z)^2)
rot_rep1$ser268A_norm <- rot_rep1$chainA_dist - rot_rep1$chainA_dist[1]
rot_rep1$ser268B_norm <- rot_rep1$chainB_dist - rot_rep1$chainB_dist[1]

rot_rep2 <- data.frame("time" = rot_rep2_chainA$time/1000, 
                       "chainA_x" = rot_rep2_chainA$x, "chainA_y" = rot_rep2_chainA$y, "chainA_z" = rot_rep2_chainA$z,
                       "chainB_x" = rot_rep2_chainB$x, "chainB_y" = rot_rep2_chainB$y, "chainB_z" = rot_rep2_chainB$z, 
                       "protein_x" = rot_rep2_protein$x,"protein_y" = rot_rep2_protein$y, "protein_z" = rot_rep2_protein$z)
rot_rep2$chainA_dist <- sqrt((rot_rep2$chainA_x - rot_rep2$protein_x)^2 + 
                               (rot_rep2$chainA_y - rot_rep2$protein_y)^2 + 
                               (rot_rep2$chainA_z - rot_rep2$protein_z)^2)
rot_rep2$chainB_dist <- sqrt((rot_rep2$chainB_x - rot_rep2$protein_x)^2 + 
                               (rot_rep2$chainB_y - rot_rep2$protein_y)^2 + 
                               (rot_rep2$chainB_z - rot_rep2$protein_z)^2)
rot_rep2$ser268A_norm <- rot_rep2$chainA_dist - rot_rep2$chainA_dist[1]
rot_rep2$ser268B_norm <- rot_rep2$chainB_dist - rot_rep2$chainB_dist[1]

distances3 <- data.frame("Time" = c(rep1$time, rep1$time, rep1$time, rep1$time),
                        "Distance" = c(rot_rep1$ser268A_norm, rot_rep1$ser268B_norm, 
                                       rot_rep2$ser268A_norm, rot_rep2$ser268B_norm),
                        "Label" = "label")
rm(rot_rep1, rot_rep1_chainA, rot_rep1_chainB, rot_rep1_hbondA, rot_rep1_hbondB, rot_rep1_protein, 
   rot_rep2, rot_rep2_chainA, rot_rep2_chainB, rot_rep2_protein)

distances3$Label[1:1001] <- "Rot Rep1 - Chain A" 
distances3$Label[1002:2002] <- "Rot Rep1 - Chain B" 
distances3$Label[2003:3003] <- "Rot Rep2 - Chain A" 
distances3$Label[3004:4004] <- "Rot Rep2 - Chain B" 

rotrep1a <- ggplot(distances3[1:2002,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("steelblue", "seagreen3")) + a + ylim(-0.25, 0.25) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of STING")

rotrep2a <- ggplot(distances3[2003:4004,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("darkgoldenrod1", "springgreen")) + a + ylim(-0.25, 0.25) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of STING")

rotrep1c <- ggplot(distances[12013:13013,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "steelblue") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 A")

rotrep1d <- ggplot(distances[13014:14014,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "seagreen3") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 B")
