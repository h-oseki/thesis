setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

rmsd1 <- read.csv("HOLO-WTSTING-control/prod/control_rmsd.xvg", sep="")
rmsd2 <- read.csv("cerrada-rep2/cerrada-rep2-rmsd.xvg", sep="")
rmsd3 <- read.csv("cerrada-rep3/cerrada-rep3-rmsd.xvg", sep="")
rmsd <- data.frame("rmsd" = c(rmsd1$rmsd, rmsd2$rmsd, rmsd3$rmsd),
                   "time" = c(rmsd1$time, rmsd2$time, rmsd3$time),
                   "Structure" = "label")
rmsd$Structure[1:1001] <- "replica 1"
rmsd$Structure[1002:2002] <- "replica 2"
rmsd$Structure[2003:3003] <- "replica 3"
rm(rmsd1, rmsd2, rmsd3)

gr1 <- read.csv("HOLO-WTSTING-control/prod/control_rg.xvg", sep="")
gr2 <- read.csv("cerrada-rep2/cerrada-rep2-rg.xvg", sep="")
gr3 <- read.csv("cerrada-rep3/cerrada-rep3-rg.xvg", sep="")
gr <- data.frame("time" = c(gr1$time, gr2$time, gr3$time),
                 "all" = c(gr1$all, gr2$all, gr3$all),
                 "x" = c(gr1$x, gr2$x, gr3$x),
                 "y" = c(gr1$y, gr2$y, gr3$y),
                 "z" = c(gr1$z, gr2$z, gr3$z),
                 "Structure" = "label")

gr$time <- gr$time/1000
gr$Structure[1:1001] <- "replica 1"
gr$Structure[1002:2002] <- "replica 2"
gr$Structure[2003:3003] <- "replica 3"
rm(gr1, gr2, gr3)

bfac1 <- read.pdb("HOLO-WTSTING-control/prod/control_bfac.pdb")
bfac2 <- read.pdb("cerrada-rep2/cerrada-rep2-bfac.pdb")
bfac3 <- read.pdb("cerrada-rep3/cerrada-rep3-bfac.pdb")
bfac1_db <- data.frame("NUM" = bfac1$atom$resno, 
                              "RES" = bfac1$atom$resid, 
                              "BFAC" = bfac1$atom$b)
prom1 <- mean(bfac1_db$BFAC)
sd1 <- sd(bfac1_db$BFAC)
bfac1_db$NBFAC <- (bfac1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  bfac1_db$ATOM[i] = j
  j = j+1
}
bfac2_db <- data.frame("NUM" = bfac2$atom$resno, 
                       "RES" = bfac2$atom$resid, 
                       "BFAC" = bfac2$atom$b)
prom2 <- mean(bfac2_db$BFAC)
sd2 <- sd(bfac2_db$BFAC)
bfac2_db$NBFAC <- (bfac2_db$BFAC - prom2)/sd2
j = 146
for (i in 1:394){
  bfac2_db$ATOM[i] = j
  j = j+1
}
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
bfacs <- data.frame("res" = c(bfac1_db$ATOM, bfac2_db$ATOM, bfac3_db$ATOM),
                    "nbfac" = c(bfac1_db$NBFAC, bfac2_db$NBFAC, bfac3_db$NBFAC),
                    "Structure" = "label")
bfacs$Structure[1:394] <- "replica 1"
bfacs$Structure[395:788] <- "replica 2"
bfacs$Structure[789:1182] <- "replica 3"
rm(i, j, prom1, prom2, sd1, sd2, bfac1, bfac1_db, bfac2, bfac2_db, bfac3, bfac3_db)

png("replicas-rmsd.png", height = 200, width = 300, units = "mm", res = 300)
rmsd_reps <- ggplot(rmsd, aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("limegreen", "mediumorchid2", "firebrick2")) +
  scale_linetype_manual(values=c("dashed", "solid", "twodash")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
dev.off()
png("replicas-gr.png", height = 200, width = 300, units = "mm", res = 300)
gr_reps <- ggplot(gr, aes(x = time, y = all)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("limegreen", "mediumorchid2", "firebrick2")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  scale_linetype_manual(values=c("dashed", "solid", "twodash")) + a + ggtitle("Radius of gyration")
dev.off()
bfacs_reps_A <- ggplot(bfacs[c(1:197, 395:591, 789:985),], aes(x = res, y = nbfac)) + ylim(-0.6, 9) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("limegreen", "mediumorchid2", "firebrick2")) +
  scale_linetype_manual(values=c("twodash", "solid", "twodash")) + a +
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
bfacs_reps_B <- ggplot(bfacs[c(198:394, 592:788, 986:1182),], aes(x = res, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("limegreen", "mediumorchid2", "firebrick2")) +
  scale_linetype_manual(values=c("twodash", "solid", "twodash")) + a + ylim(-0.6, 9) +
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
png("replicas-bfacs.png", height = 150, width = 300, units = "mm", res = 300)
ggarrange(A, B,  labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom") 
dev.off()


rep1_chainA <- read.delim("HOLO-WTSTING-control/prod/distances/chainA-rep1.xvg")
rep1_chainB <- read.delim("HOLO-WTSTING-control/prod/distances/chainB-rep1.xvg")
rep1_protein <- read.delim("HOLO-WTSTING-control/prod/distances/protein-rep1.xvg")
rep1_ser268a <- read.delim("HOLO-WTSTING-control/prod/distances/ser268a-rep1.xvg")
rep1_ser268b <- read.delim("HOLO-WTSTING-control/prod/distances/ser268b-rep1.xvg")
rep1_lig <- read.delim("HOLO-WTSTING-control/prod/distances/lig-rep1.xvg")
rep1_hbondA <- read.table("HOLO-WTSTING-control/prod/distances/ser268a-lig.dat", quote="\"", comment.char="")
rep1_hbondA <- rep1_hbondA[-1,]
rep1_hbondB <- read.table("HOLO-WTSTING-control/prod/distances/ser268b-lig.dat", quote="\"", comment.char="")
rep1_hbondB <- rep1_hbondB[-1,]

rep2_chainA <- read.delim("cerrada-rep2/distances/chainA-rep2.xvg")
rep2_chainB <- read.delim("cerrada-rep2/distances/chainB-rep2.xvg")
rep2_protein <- read.delim("cerrada-rep2/distances/protein-rep2.xvg")
rep2_ser268a <- read.delim("cerrada-rep2/distances/ser268a-rep2.xvg")
rep2_ser268b <- read.delim("cerrada-rep2/distances/ser268b-rep2.xvg")
rep2_lig <- read.delim("cerrada-rep2/distances/lig-rep2.xvg")
rep2_hbondA <- read.table("cerrada-rep2/distances/ser268a-lig.dat", quote="\"", comment.char="")
rep2_hbondA <- rep2_hbondA[-1,]
rep2_hbondB <- read.table("cerrada-rep2/distances/ser268b-lig.dat", quote="\"", comment.char="")
rep2_hbondB <- rep2_hbondB[-1,]


rep3_chainA <- read.delim("cerrada-rep3/distances/chainA-rep3.xvg")
rep3_chainB <- read.delim("cerrada-rep3/distances/chainB-rep3.xvg")
rep3_protein <- read.delim("cerrada-rep3/distances/protein-rep3.xvg")
rep3_ser268a <- read.delim("cerrada-rep3/distances/ser268a-rep3.xvg")
rep3_ser268b <- read.delim("cerrada-rep3/distances/ser268b-rep3.xvg")
rep3_lig <- read.delim("cerrada-rep3/distances/lig-rep3.xvg")
rep3_hbondA <- read.table("cerrada-rep3/distances/ser268a-lig.dat", quote="\"", comment.char="")
rep3_hbondA <- rep3_hbondA[-1,]
rep3_hbondB <- read.table("cerrada-rep3/distances/ser268b-lig.dat", quote="\"", comment.char="")
rep3_hbondB <- rep3_hbondB[-1,]

rep1 <- data.frame("time" = rep1_chainA$time/1000, 
                   "chainA_x" = rep1_chainA$x,"chainA_y" = rep1_chainA$y, "chainA_z" = rep1_chainA$z,
                   "chainB_x" = rep1_chainB$x, "chainB_y" = rep1_chainB$y,"chainB_z" = rep1_chainB$z, 
                   "protein_x" = rep1_protein$x,"protein_y" = rep1_protein$y, "protein_z" = rep1_protein$z,
                   "ser268a_x" = rep1_ser268a$x,"ser268a_y" = rep1_ser268a$y, "ser268a_z" = rep1_ser268a$z,
                   "ser268b_x" = rep1_ser268b$x,"ser268b_y" = rep1_ser268b$y, "ser268b_z" = rep1_ser268b$z,
                   "lig_x" = rep1_lig$x, "lig_y" = rep1_lig$y, "lig_z" = rep1_lig$z,
                   "hbondA" = rep1_hbondA$V2, "hbondB" = rep1_hbondB$V2)


rep2 <- data.frame("time" = rep2_chainA$time/1000, 
                   "chainA_x" = rep2_chainA$x,"chainA_y" = rep2_chainA$y, "chainA_z" = rep2_chainA$z,
                   "chainB_x" = rep2_chainB$x, "chainB_y" = rep2_chainB$y,"chainB_z" = rep2_chainB$z, 
                   "protein_x" = rep2_protein$x,"protein_y" = rep2_protein$y, "protein_z" = rep2_protein$z,
                   "ser268a_x" = rep2_ser268a$x,"ser268a_y" = rep2_ser268a$y, "ser268a_z" = rep2_ser268a$z,
                   "ser268b_x" = rep2_ser268b$x,"ser268b_y" = rep2_ser268b$y, "ser268b_z" = rep2_ser268b$z,
                   "lig_x" = rep2_lig$x, "lig_y" = rep2_lig$y, "lig_z" = rep2_lig$z,
                   "hbondA" = rep2_hbondA$V2, "hbondB" = rep2_hbondB$V2)

rep3 <- data.frame("time" = rep3_chainA$time/1000, 
                   "chainA_x" = rep3_chainA$x,"chainA_y" = rep3_chainA$y, "chainA_z" = rep3_chainA$z,
                   "chainB_x" = rep3_chainB$x, "chainB_y" = rep3_chainB$y,"chainB_z" = rep3_chainB$z, 
                   "protein_x" = rep3_protein$x,"protein_y" = rep3_protein$y, "protein_z" = rep3_protein$z,
                   "ser268a_x" = rep3_ser268a$x,"ser268a_y" = rep3_ser268a$y, "ser268a_z" = rep3_ser268a$z,
                   "ser268b_x" = rep3_ser268b$x,"ser268b_y" = rep3_ser268b$y, "ser268b_z" = rep3_ser268b$z,
                   "lig_x" = rep3_lig$x, "lig_y" = rep3_lig$y, "lig_z" = rep3_lig$z,
                   "hbondA" = rep3_hbondA$V2, "hbondB" = rep3_hbondB$V2)

rep1$chainA_dist <- sqrt((rep1$chainA_x - rep1$protein_x)^2 + 
                           (rep1$chainA_y - rep1$protein_y)^2 + 
                           (rep1$chainA_z - rep1$protein_z)^2)
rep1$chainB_dist <- sqrt((rep1$chainB_x - rep1$protein_x)^2 + 
                           (rep1$chainB_y - rep1$protein_y)^2 + 
                           (rep1$chainB_z - rep1$protein_z)^2)
rep1$ser268A_dist <- sqrt((rep1$ser268a_x - rep1$lig_x)^2 + 
                            (rep1$ser268a_y - rep1$lig_y)^2 +
                            (rep1$ser268a_z - rep1$lig_z)^2)
rep1$ser268B_dist <- sqrt((rep1$ser268b_x - rep1$lig_x)^2 + 
                            (rep1$ser268b_y - rep1$lig_y)^2 +
                            (rep1$ser268b_z - rep1$lig_z)^2)
rep1$ser268A_norm <- rep1$ser268A_dist - rep1$ser268A_dist[1]
rep1$ser268B_norm <- rep1$ser268B_dist - rep1$ser268B_dist[1]

rep2$chainA_dist <- sqrt((rep2$chainA_x - rep2$protein_x)^2 + 
                           (rep2$chainA_y - rep2$protein_y)^2 + 
                           (rep2$chainA_z - rep2$protein_z)^2)
rep2$chainB_dist <- sqrt((rep2$chainB_x - rep2$protein_x)^2 + 
                           (rep2$chainB_y - rep2$protein_y)^2 + 
                           (rep2$chainB_z - rep2$protein_z)^2)
rep2$ser268A_dist <- sqrt((rep2$ser268a_x - rep2$lig_x)^2 + 
                            (rep2$ser268a_y - rep2$lig_y)^2 +
                            (rep2$ser268a_z - rep2$lig_z)^2)
rep2$ser268B_dist <- sqrt((rep2$ser268b_x - rep2$lig_x)^2 + 
                            (rep2$ser268b_y - rep2$lig_y)^2 +
                            (rep2$ser268b_z - rep2$lig_z)^2)
rep2$ser268A_norm <- rep2$ser268A_dist - rep2$ser268A_dist[1]
rep2$ser268B_norm <- rep2$ser268B_dist - rep2$ser268B_dist[1]

rep3$chainA_dist <- sqrt((rep3$chainA_x - rep3$protein_x)^2 + 
                           (rep3$chainA_y - rep3$protein_y)^2 + 
                           (rep3$chainA_z - rep3$protein_z)^2)
rep3$chainB_dist <- sqrt((rep3$chainB_x - rep3$protein_x)^2 + 
                           (rep3$chainB_y - rep3$protein_y)^2 + 
                           (rep3$chainB_z - rep3$protein_z)^2)
rep3$ser268A_dist <- sqrt((rep3$ser268a_x - rep3$lig_x)^2 + 
                            (rep3$ser268a_y - rep3$lig_y)^2 +
                            (rep3$ser268a_z - rep3$lig_z)^2)
rep3$ser268B_dist <- sqrt((rep3$ser268b_x - rep3$lig_x)^2 + 
                            (rep3$ser268b_y - rep3$lig_y)^2 +
                            (rep3$ser268b_z - rep3$lig_z)^2)
rep3$ser268A_norm <- rep3$ser268A_dist - rep3$ser268A_dist[1]
rep3$ser268B_norm <- rep3$ser268B_dist - rep3$ser268B_dist[1]

distances <- data.frame("Time" = c(rep1$time, rep1$time, rep1$time, rep1$time, rep1$time, rep1$time,
                                   rep1$time, rep1$time, rep1$time, rep1$time, rep1$time, rep1$time,
                                   rep1$time, rep1$time, rep1$time, rep1$time, rep1$time, rep1$time),
                        "Distance" = c(rep1$ser268A_norm, rep1$ser268B_norm, 
                                       rep2$ser268A_norm, rep2$ser268B_norm,
                                       rep3$ser268A_norm, rep3$ser268B_norm,
                                       rep1$ser268A_dist, rep1$ser268B_dist,
                                       rep2$ser268A_dist, rep2$ser268B_dist,
                                       rep3$ser268A_dist, rep3$ser268B_dist,
                                       rep1_hbondA$V2, rep1_hbondB$V2,
                                       rep2_hbondA$V2, rep2_hbondB$V2,
                                       rep3_hbondA$V2, rep3_hbondB$V2),
                        "Label" = "label")


rm(rep1_hbondA, rep1_hbondB, rep2_hbondA, rep2_hbondB, rep3_hbondA, rep3_hbondB)
rm(rep1_chainA, rep1_chainB, rep1_protein, rep1_lig, rep1_ser268a, rep1_ser268b, 
   rep2_chainA, rep2_chainB, rep2_protein, rep2_lig, rep2_ser268a, rep2_ser268b,
   rep3_chainA, rep3_chainB, rep3_protein, rep3_lig, rep3_ser268a, rep3_ser268b)

distances$Label[1:1001] <- "Rep1 - Chain A" #deeppink
distances$Label[1002:2002] <- "Rep1 - Chain B" #chartreuse3
distances$Label[2003:3003] <- "Rep2 - Chain A" #red3
distances$Label[3004:4004] <- "Rep2 - Chain B" #purple
distances$Label[4005:5005] <- "Rep3 - Chain A" #royalblue
distances$Label[5006:6006] <- "Rep3 - Chain B" #orangered
distances$Label[6007:7007] <- "Rep1 - SER268A" #deeppink
distances$Label[7008:8008] <- "Rep1 - SER268B" #chartreuse3
distances$Label[8009:9009] <- "Rep2 - SER268A" #red3
distances$Label[9010:10010] <- "Rep2 - SER268B" #purple
distances$Label[10011:11011] <- "Rep3 - SER268A" #royalblue
distances$Label[11012:12012] <- "Rep3 - SER268B" #orangered
distances$Label[12013:13013] <- "Rep1 - HB SER268A - cGAMP" #deeppink
distances$Label[13014:14014] <- "Rep1 - HB SER268B - cGAMP" #chartreuse3
distances$Label[14015:15015] <- "Rep2 - HB SER268A - cGAMP" #red3
distances$Label[15016:16016] <- "Rep2 - HB SER268B - cGAMP" #purple
distances$Label[16017:17017] <- "Rep3 - HB SER268A - cGAMP" #royalblue
distances$Label[17018:18018] <- "Rep3 - HB SER268B - cGAMP" #orangered

A <- ggplot(distances[1:6006,], aes(x = Time, y = Distance, color = Label)) + geom_line(size = 0.8) + a

B <- ggplot(distances[6007:12012,], aes(x = Time, y = Distance, color = Label)) + geom_line(size = 0.8) + a

ggarrange(A, B)

#### Replica 1 ####

rep1a <- ggplot(distances[1:2002,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "chartreuse3")) + a + 
  ggtitle("Replic 1") + ylim(-0.25, 0.25) 
  

rep1b <- ggplot(distances[6007:8008,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "chartreuse3")) + a + 
  ggtitle("Distance between cGAMP and Serine 268")

rep1c <- ggplot(distances[12013:13013,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "deeppink") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 A")

rep1d <- ggplot(distances[13014:14014,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "chartreuse3") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 B")

rep1e <- ggplot(rep1, aes(x = chainA_dist, y = ser268A_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Time (ns)")

rep1f <- ggplot(rep1, aes(x = chainA_dist, y = ser268A_dist, color = as.factor(hbondA))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

rep1g <- ggplot(rep1, aes(x = chainB_dist, y = ser268B_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Time (ns)")

rep1h <- ggplot(rep1, aes(x = chainB_dist, y = ser268B_dist, color = as.factor(hbondB))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

png("rep1-dist1.png", height = 500, width = 300, units = "mm", res = 300)
ggarrange(rep1a, rep1b, rep1c, rep1d, 
          common.legend = TRUE, legend = "bottom",
          ncol = 1, nrow = 4)
dev.off()
png("rep1-dist2.png", height = 500, width = 500, units = "mm", res = 300)
ggarrange(ggarrange(rep1e, rep1g, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ggarrange(rep1f, rep1h, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ncol = 2, nrow = 1)
dev.off()

#### Replica 2 ####
rep2a <- ggplot(distances[2003:4004,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("red3", "purple")) + a + 
  ggtitle("Replic 2") + ylim(-0.25, 0.25)

rep2b <- ggplot(distances[8009:10010,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("red3", "purple")) + a + 
  ggtitle("Distance between cGAMP and Serine 268")

rep2c <- ggplot(distances[14015:15015,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "red3") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 A")

rep2d <- ggplot(distances[15016:16016,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "purple") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 B")

rep2e <- ggplot(rep2, aes(x = chainA_dist, y = ser268A_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Time (ns)")

rep2f <- ggplot(rep2, aes(x = chainA_dist, y = ser268A_dist, color = as.factor(hbondA))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

rep2g <- ggplot(rep2, aes(x = chainB_dist, y = ser268B_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Time (ns)")

rep2h <- ggplot(rep2, aes(x = chainB_dist, y = ser268B_dist, color = as.factor(hbondB))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))
png("rep2-dist1.png", height = 500, width = 300, units = "mm", res = 300)
ggarrange(rep2a, rep2b, rep2c, rep2d, 
          common.legend = TRUE, legend = "bottom",
          ncol = 1, nrow = 4)
dev.off()
png("rep2-dist2.png", height = 500, width = 500, units = "mm", res = 300)
ggarrange(ggarrange(rep2e, rep2g, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ggarrange(rep2f, rep2h, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ncol = 2, nrow = 1)
dev.off()
#### Replica 3 ####
rep3a <- ggplot(distances[4005:6006,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("royalblue", "orangered")) + a + 
  ggtitle("Replic 3") + ylim(-0.25, 0.25) 

rep3b <- ggplot(distances[10011:12012,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("royalblue", "orangered")) + a + 
  ggtitle("Distance between cGAMP and Serine 268")

rep3c <- ggplot(distances[16017:17017,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "royalblue") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 A")

rep3d <- ggplot(distances[17018:18018,], aes(x = Time, y = Distance)) + a + 
  geom_line(size = .8, color = "orangered") + xlab("Time (ns)") + ylab("Occurrence") + 
  ggtitle("Number of hydrogen bonds between cGAMP and Serine 268 B")

rep3e <- ggplot(rep3, aes(x = chainA_dist, y = ser268A_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Time (ns)")

rep3f <- ggplot(rep3, aes(x = chainA_dist, y = ser268A_dist, color = as.factor(hbondA))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain A and STING com") + 
  ylab("Distance between Serine 268 A and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))

rep3g <- ggplot(rep3, aes(x = chainB_dist, y = ser268B_dist, color = time)) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Time (ns)")

rep3h <- ggplot(rep3, aes(x = chainB_dist, y = ser268B_dist, color = as.factor(hbondB))) + 
  geom_point(size = 0.8) + xlab("Distance between polymerization site of chain B and STING com") + 
  ylab("Distance between Serine 268 B and cGAMP") + a + labs(color = "Hydrogen\nbonds") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A"))
png("rep3-dist1.png", height = 500, width = 300, units = "mm", res = 300)
ggarrange(rep3a, rep3b, rep3c, rep3d, 
          common.legend = TRUE, legend = "bottom",
          ncol = 1, nrow = 4)
dev.off()
png("rep3-dist2.png", height = 500, width = 500, units = "mm", res = 300)
ggarrange(ggarrange(rep3e, rep3g, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ggarrange(rep3f, rep3h, ncol = 1, nrow = 2,
                    common.legend = TRUE, legend = "bottom"),
          ncol = 2, nrow = 1)
dev.off()
  