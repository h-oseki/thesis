abierta_chainA <- read.delim("6nt6_NO_BORRAR/distance/chainA-abierta.xvg")
abierta_chainB <- read.delim("6nt6_NO_BORRAR/distance/chainB-abierta.xvg")
abierta_protein <- read.delim("6nt6_NO_BORRAR/distance/protein-abierta.xvg")
rms_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt6)[names(rms_6nt6)=="TIME.ps."] <- "Time (ns)"
names(rms_6nt6)[names(rms_6nt6)=="RMSD.nm."] <- "RMSD (nm)"
rms_6nt6$`Time (ns)` = rms_6nt6$`Time (ns)`/1000
gr_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_GR.xvg", sep="")
gr_6nt6$Time = gr_6nt6$Time/1000

abierta_rmsd <- ggplot(rms_6nt6, aes(x=rms_6nt6$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt6$`RMSD (nm)`), color = "mediumseagreen", size = 0.8) + 
  ggtitle("Root mean square deviation") + ylim(0,0.6) +
  xlab("Time (ns)") + 
  ylab("RMSD (nm)") + 
  a

abierta_gr <- ggplot(gr_6nt6, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt6$RG), color="mediumseagreen", size=0.8) + 
  ggtitle("Radius of gyration") + xlab("Time (ns)") +   ylab("RG (nm)") + a + ylim(2.15,2.35) +
  annotate("rect", xmin = 60, xmax=105, ymin=2.15, ymax= 2.35, alpha=.2, fill="mediumseagreen") +
  annotate("rect", xmin = 290, xmax=410, ymin=2.15, ymax= 2.35, alpha=.2, fill="mediumseagreen")

abierta <- data.frame("time" = abierta_chainA$time/1000, 
                       "chainA_x" = abierta_chainA$x, "chainA_y" = abierta_chainA$y, "chainA_z" = abierta_chainA$z,
                       "chainB_x" = abierta_chainB$x, "chainB_y" = abierta_chainB$y, "chainB_z" = abierta_chainB$z, 
                       "protein_x" = abierta_protein$x,"protein_y" = abierta_protein$y, "protein_z" = abierta_protein$z)
abierta$chainA_dist <- sqrt((abierta$chainA_x - abierta$protein_x)^2 + 
                               (abierta$chainA_y - abierta$protein_y)^2 + 
                               (abierta$chainA_z - abierta$protein_z)^2)
abierta$chainB_dist <- sqrt((abierta$chainB_x - abierta$protein_x)^2 + 
                               (abierta$chainB_y - abierta$protein_y)^2 + 
                               (abierta$chainB_z - abierta$protein_z)^2)  
abierta$dist_normA <- abierta$chainA_dist - abierta$chainA_dist[1]
abierta$dist_normB <- abierta$chainB_dist - abierta$chainB_dist[1]

distances_abierta <- data.frame("Time" = c(abierta_chainA$time, abierta_chainA$time),
                                 "Distance" = c(abierta$dist_normA, abierta$dist_normB),
                                 "Label" = "label")
rm(abierta_chainA, abierta_chainB, abierta_protein, abierta)

distances_abierta$Label[1:2501] <- "Abeirta - Chain A" 
distances_abierta$Label[2502:5002] <- "Abierta - Chain B" 

abierta_dist <- ggplot(distances_abierta[1:5002,], aes(x = Time/1000, y = Distance, color = Label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("darkmagenta", "chartreuse4")) + a + ylim(-0.25, 0.25) +
  ggtitle("APO STING") 


bfac_6nt6 <- read.pdb("6nt6_NO_BORRAR/400ns_NpT-WT_6nt6_nojump_rot_bfac.pdb")
bfac_6nt7_holo <- read.pdb("6nt7_HOLO_NO_BORRAR/400ns_6nt7_holo_NpT_nojump_rot_bfac.pdb")
lig_rot_bfac <- read.pdb("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-bfac.pdb")
bfac_6nt6_db <- data.frame("NUM" = bfac_6nt6$atom$resno, 
                           "RES" = bfac_6nt6$atom$resid, 
                           "BFAC" = bfac_6nt6$atom$b)
prom_6nt6 <- mean(bfac_6nt6_db$BFAC)
sd_6nt6 <- sd(bfac_6nt6_db$BFAC)
bfac_6nt6_db$NBFAC <- (bfac_6nt6_db$BFAC - prom_6nt6)/sd_6nt6
bfac_6nt6_db <- bfac_6nt6_db[-c(198:201, 399:402),]
j = 146
for (i in 1:394){
  bfac_6nt6_db$ATOM[i] = j
  j = j+1
}
bfac_6nt7_holo_db <- data.frame("NUM" = bfac_6nt7_holo$atom$resno, 
                                "RES" = bfac_6nt7_holo$atom$resid, 
                                "BFAC" = bfac_6nt7_holo$atom$b)
prom_6nt7_holo <- mean(bfac_6nt7_holo_db$BFAC)
sd_6nt7_holo <- sd(bfac_6nt7_holo_db$BFAC)
bfac_6nt7_holo_db$NBFAC <- (bfac_6nt7_holo_db$BFAC - prom_6nt7_holo)/sd_6nt7_holo
j = 146
for (i in 1:394){
  bfac_6nt7_holo_db$ATOM[i] = j
  j = j+1
}

bfac_lig_rot_db <- data.frame("NUM" = lig_rot_bfac$atom$resno, 
                              "RES" = lig_rot_bfac$atom$resid, 
                              "BFAC" = lig_rot_bfac$atom$b)
prom2 <- mean(bfac_lig_rot_db$BFAC)
sd2 <- sd(bfac_lig_rot_db$BFAC)
bfac_lig_rot_db$NBFAC <- (bfac_lig_rot_db$BFAC - prom2)/sd2
j = 146
for (i in 1:394){
  bfac_lig_rot_db$ATOM[i] = j
  j = j+1
}
bfac <- data.frame("residue" = c(bfac_6nt6_db$ATOM, bfac_6nt7_holo_db$ATOM, bfac_lig_rot_db$ATOM), 
                   "nbfac" = c(bfac_6nt6_db$NBFAC, bfac_6nt7_holo_db$NBFAC, bfac_lig_rot_db$NBFAC),
                   "Structure" = "label")
bfac$Structure <- as.character(bfac$Structure)
bfac$Structure[c(1:197)] = as.character("APO-STING")
bfac$Structure[c(198:394)] = as.character("APO-STING")
bfac$Structure[c(395:591)] = as.character("HOLO-STING")
bfac$Structure[c(592:788)] = as.character("HOLO-STING")
bfac$Structure[c(789:985)] = as.character("HOLO-STING rot")
bfac$Structure[c(986:1182)] = as.character("HOLO-STING rot")
rm(bfac_6nt6, bfac_6nt6_db, bfac_6nt7_holo, bfac_6nt7_holo_db, bfac_lig_rot_db, lig_rot_bfac)
rm(i, j, prom_6nt6, prom_6nt7_holo, prom2, sd_6nt6, sd_6nt7_holo, sd2)


