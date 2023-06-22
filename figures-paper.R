#################################################################################################################
#################################################################################################################
#####                                                                                                       #####
#####                                       SCRIPT PAPER - FIGURES                                          #####
#####                                                                                                       #####
#####   1- FILES                                                                                            #####
#####   2- FIGURA 1                                                                                         #####
#####   3- FIGURA 2                                                                                         #####
#####   4- FIGURA 3                                                                                         #####
#####   5- FIGURA 4                                                                                         #####
#####   6- FIGURA 5                                                                                         #####
#####   7- FIGURA 6                                                                                         #####
#####   8- FIGURA S1                                                                                        #####
#####   9- FIGURA S2                                                                                        #####
#####   10- FIGURA S3                                                                                       #####
#####   11- FIGURA S5                                                                                       #####
#####                                                                                             -hoseki   #####
#################################################################################################################
#################################################################################################################

setwd("~/Desktop/Sysbio/cSTING")

#THEME FOR PLOTS
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)
#ggarrange(A, B,  labels = c("A", "B"),
#          common.legend = TRUE, legend = "bottom",
#          widths = c(2,3)) 


################################################# 1- FILES  ######################################################
### RMSD
rms_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt6)[names(rms_6nt6)=="TIME.ps."] <- "Time"
names(rms_6nt6)[names(rms_6nt6)=="RMSD.nm."] <- "RMSD"
rms_6nt6$Time = rms_6nt6$Time/1000
rms_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="Time"] <- "Time"
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="RMSDns"] <- "RMSD"
lig_rot_rmsd <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-rmsd.xvg", sep="")
rms <- data.frame("time" = c(rms_6nt6$Time, rms_6nt7_holo$Time, lig_rot_rmsd$time),
                  "rmsd" = c(rms_6nt6$RMSD, rms_6nt7_holo$RMSD, lig_rot_rmsd$rmsd),
                  "Structure" = as.character("label"))
rms$Structure[0:2501] <- "APO-STING"
rms$Structure[2502:5002] <- "HOLO-STING"
rms$Structure[5003:7503] <- "HOLO-STING rot-lig"
rm(lig_rot_rmsd, rms_6nt6, rms_6nt7_holo)

### RADIUS OF GYRATION
gr_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_GR.xvg", sep="")
gr_6nt6$Time = gr_6nt6$Time/1000
gr_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_gr.xvg", sep="")
gr_6nt7_holo$Time = gr_6nt7_holo$Time/1000
lig_rot_gr <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-gr.xvg", sep="")
lig_rot_gr$time <- lig_rot_gr$time/1000
gr <- data.frame("time" = c(gr_6nt6$Time, gr_6nt7_holo$Time, lig_rot_gr$time),
                 "gr" = c(gr_6nt6$RG, gr_6nt7_holo$GR, lig_rot_gr$all), 
                 "Structure" = as.character("label"))
gr$Structure[0:2501] <- "APO-STING"
gr$Structure[2502:5002] <- "HOLO-STING"
gr$Structure[5003:7503] <- "HOLO-STING rot-lig"
rm(gr_6nt6, gr_6nt7_holo, lig_rot_gr)

### B-FACTORS
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

### PRINCIPAL COMPONENT ANALYSIS
hbonds <- read.csv("6nt7_HOLO_NO_BORRAR/holo_num_hbond.xvg", sep="")
hbonds[2502,] = hbonds[2501,]
hbonds <- hbonds[-c(1:500),]

pdb_6nt6 <- read.pdb("6nt6_NO_BORRAR/6nt6_ns100-500.pdb")
dcd_6n6 <- read.dcd("6nt6_NO_BORRAR/6nt6_ns100-500.dcd")
ca_ndx_6nt6 <- atom.select(pdb_6nt6, elety = "CA") 
pca_6nt6 <- pca.xyz(dcd_6n6[,ca_ndx_6nt6$xyz])
pca_6nt6_db <- data.frame("PC" = pca_6nt6$z)
pca_6nt6_db$time <- c(1:2002)
pca_6nt6_db$time[(1:2002)] = hbonds$TIME
eigen_bd_6nt6 <- data.frame("eigen" = pca_6nt6$L)
sum_eig_6nt6 = sum(eigen_bd_6nt6$eigen)
eigen_bd_6nt6$percent_6nt6 <- (eigen_bd_6nt6$eigen/sum_eig_6nt6)*100
eigen_bd_6nt6$num <- c(1:1206)
eigen_bd_6nt6$sum_percent = 0
eigen_bd_6nt6$sum_percent[1] = eigen_bd_6nt6$percent_6nt6[1]
for (i in 2:1206){
  eigen_bd_6nt6$sum_percent[i] = eigen_bd_6nt6$percent_6nt6[i] + eigen_bd_6nt6$sum_percent[i - 1]
}
contribution_6nt6 <- data.frame("num" = pdb_6nt6$atom$resno[ca_ndx_6nt6$atom],
                                "resnum" = pdb_6nt6$atom$resno[ca_ndx_6nt6$atom],
                                "res" = pdb_6nt6$atom$resid[ca_ndx_6nt6$atom], 
                                "pc1" = pca_6nt6$au[,1], 
                                "pc2" = pca_6nt6$au[,2], 
                                "pc3" = pca_6nt6$au[,3]) 
i = 146
for (j in 1:402){
  contribution_6nt6$num[j] = i
  i = i+1
}
threshold_6nt6_pc1 = (max(contribution_6nt6$pc1) + min(contribution_6nt6$pc1))/2
threshold_6nt6_pc2 = (max(contribution_6nt6$pc2) + min(contribution_6nt6$pc2))/2
rm(ca_ndx_6nt6, dcd_6n6, pca_6nt6, pdb_6nt6, i, j, sum_eig_6nt6)
contribution_6nt6 <- contribution_6nt6[-c(198:201, 399:402),]
j = 146
for(i in 1:394){
  contribution_6nt6$num[i] = j
  j = j+1
}
rm(i, j)

pdb_6nt7_holo <- read.pdb("6nt7_HOLO_NO_BORRAR/6nt7_holo_ns100-500.pdb")
dcd_6nt7_holo <- read.dcd("6nt7_HOLO_NO_BORRAR/6nt7_holo_ns100-500.dcd")
ca_ndx_6nt7_holo <- atom.select(pdb_6nt7_holo, elety = "CA")
pca_6nt7_holo <- pca.xyz(dcd_6nt7_holo[,ca_ndx_6nt7_holo$xyz])
pca_6nt7_holo_db <- data.frame("PC" = pca_6nt7_holo$z)
pca_6nt7_holo_db$time <- c(1:2002)
pca_6nt7_holo_db$time[(1:2002)] = hbonds$TIME
eigen_bd_6nt7_holo <- data.frame("eigen" = pca_6nt7_holo$L)
sum_eig_6nt7_holo = sum(eigen_bd_6nt7_holo$eigen)
eigen_bd_6nt7_holo$percent <- (eigen_bd_6nt7_holo$eigen/sum_eig_6nt7_holo)*100
eigen_bd_6nt7_holo$num <- c(1:1182)
eigen_bd_6nt7_holo$sum_percent = 0
eigen_bd_6nt7_holo$sum_percent[1] = eigen_bd_6nt7_holo$percent[1]
for (i in 2:1182){
  eigen_bd_6nt7_holo$sum_percent[i] = eigen_bd_6nt7_holo$percent[i] + eigen_bd_6nt7_holo$sum_percent[i - 1]
}
contribution_6nt7_holo <- data.frame("num" = pdb_6nt7_holo$atom$resno[ca_ndx_6nt7_holo$atom],
                                     "resnum" = pdb_6nt7_holo$atom$resno[ca_ndx_6nt7_holo$atom],
                                     "res" = pdb_6nt7_holo$atom$resid[ca_ndx_6nt7_holo$atom], 
                                     "pc1" = pca_6nt7_holo$au[,1], 
                                     "pc2" = pca_6nt7_holo$au[,2], 
                                     "pc3" = pca_6nt7_holo$au[,3]) 
i = 146
for (j in 1:394){
  contribution_6nt7_holo$num[j] = i
  i = i+1
}
threshold_6nt7_holo_pc1 = (max(contribution_6nt7_holo$pc1) + min(contribution_6nt7_holo$pc1))/2
threshold_6nt7_holo_pc2 = (max(contribution_6nt7_holo$pc2) + min(contribution_6nt7_holo$pc2))/2
rm(ca_ndx_6nt7_holo, dcd_6nt7_holo, pca_6nt7_holo, pdb_6nt7_holo, i, j, sum_eig_6nt7_holo)

lig_rot_pdb <- read.pdb("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-400ns.pdb")
lig_rot_dcd <- read.dcd("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-400ns.dcd")
ca_ndx_6nt7_holo_rot <- atom.select(lig_rot_pdb, elety = "CA")
pca_6nt7_holo_rot <- pca.xyz(lig_rot_dcd[,ca_ndx_6nt7_holo_rot$xyz])
pca_6nt7_holo_rot_db <- data.frame("PC" = pca_6nt7_holo_rot$z)
pca_6nt7_holo_rot_db$time <- c(1:2002)
pca_6nt7_holo_rot_db$time[(1:2002)] = hbonds$TIME
eigen_bd_6nt7_holo_rot <- data.frame("eigen" = pca_6nt7_holo_rot$L)
sum_eig_6nt7_holo_rot = sum(eigen_bd_6nt7_holo_rot$eigen)
eigen_bd_6nt7_holo_rot$percent <- (eigen_bd_6nt7_holo_rot$eigen/sum_eig_6nt7_holo_rot)*100
eigen_bd_6nt7_holo_rot$num <- c(1:1182)
eigen_bd_6nt7_holo_rot$sum_percent = 0
eigen_bd_6nt7_holo_rot$sum_percent[1] = eigen_bd_6nt7_holo_rot$percent[1]
for (i in 2:1182){
  eigen_bd_6nt7_holo_rot$sum_percent[i] = eigen_bd_6nt7_holo_rot$percent[i] + eigen_bd_6nt7_holo_rot$sum_percent[i - 1]}
contribution_6nt7_holo_rot <- data.frame("num" = lig_rot_pdb$atom$resno[ca_ndx_6nt7_holo_rot$atom],
                                         "resnum" = lig_rot_pdb$atom$resno[ca_ndx_6nt7_holo_rot$atom],
                                         "res" = lig_rot_pdb$atom$resid[ca_ndx_6nt7_holo_rot$atom], 
                                         "pc1" = pca_6nt7_holo_rot$au[,1], 
                                         "pc2" = pca_6nt7_holo_rot$au[,2], 
                                         "pc3" = pca_6nt7_holo_rot$au[,3]) 
i = 146
for (j in 1:394){
  contribution_6nt7_holo_rot$num[j] = i
  i = i+1
}

threshold_6nt7_holo_rot_pc1 = (max(contribution_6nt7_holo_rot$pc1) + min(contribution_6nt7_holo_rot$pc1))/2
threshold_6nt7_holo_rot_pc2 = (max(contribution_6nt7_holo_rot$pc2) + min(contribution_6nt7_holo_rot$pc2))/2
rm(ca_ndx_6nt7_holo_rot, lig_rot_dcd, lig_rot_pdb, pca_6nt7_holo_rot, i, j, sum_eig_6nt7_holo_rot, hbonds)  

### PICS 
F1A <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/fig1/abierta-tiempo0.png")
f1a <- rasterGrob(F1A, interpolate=TRUE)
F1B <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/fig1/abierta-tiempo400.png")
f1b <- rasterGrob(F1B, interpolate=TRUE)
F1C <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/fig1/cerrada-tiempo0.png")
f1c <- rasterGrob(F1C, interpolate=TRUE)
F1D <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/fig1/cerrada-tiempo400.png")
f1d <- rasterGrob(F1D, interpolate=TRUE)

F2C <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/abierta-bfac.png")
F2D <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/cerrada-bfac.png")
F2E <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/abierta-bfac-zoom4.png")
F2F <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/figures/cerrada-bfac-zoom4.png")
f2c <- rasterGrob(F2C, interpolate=TRUE)
f2d <- rasterGrob(F2D, interpolate=TRUE)
f2e <- rasterGrob(F2E, interpolate=TRUE)
f2f <- rasterGrob(F2F, interpolate=TRUE)

F3C <- readPNG("figures/holo-2pc.png")
f3c <- rasterGrob(F3C, interpolate=TRUE)
F3D <- readPNG("figures/holo-2pc3.png")
f3d <- rasterGrob(F3D, interpolate=TRUE)

F4A <- readPNG("figures/hbonds-net4.png")
f4a <- rasterGrob(F4A, interpolate = TRUE)
F4B <- readPNG("figures/Hydrogen Bonds-F4B.png")
f4b <- rasterGrob(F4B, interpolate = TRUE)
F4C <- readPNG("figures/hbonds/chainA-time420-l.png")
f4c <- rasterGrob(F4C, interpolate = TRUE)
F4D <- readPNG("figures/hbonds/chainB-time420-l.png")
f4d <- rasterGrob(F4D, interpolate = TRUE)

F5A <- readPNG("figures/cgamps-derecho.png")
f5a <- rasterGrob(F5A, interpolate = TRUE)
F5Abis <- readPNG("figures/cgamps-rotado.png")
f5abis <- rasterGrob(F5Abis, interpolate = TRUE)
F5B <- readPNG("figures/cgamps2.png")
f5b <- rasterGrob(F5B, interpolate = TRUE)
F5E <- readPNG("figures/Hydrogen Bonds-F5E.png")
f5e <- rasterGrob(F5E, interpolate = TRUE)

F6A <- readPNG("figures/dna/cerrada-chaina3.png")
f6a <- rasterGrob(F6A, interpolate = TRUE)
F6B <- readPNG("figures/dna/rotada-chaina3.png")
f6b <- rasterGrob(F6B, interpolate = TRUE)
F6C <- readPNG("figures/dna/cerrada-chainb3.png")
f6c <- rasterGrob(F6C, interpolate = TRUE)
F6D <- readPNG("figures/dna/rotada-chainb3.png")
f6d <- rasterGrob(F6D, interpolate = TRUE)

FS1A <- readPNG("figures/fig1/pdb-6nt6-6.png")
fs1a <- rasterGrob(FS1A, interpolate = TRUE)
FS1C <- readPNG("figures/fig1/pdb-6nt7-6.png")
fs1c <- rasterGrob(FS1C, interpolate = TRUE)

FS2B <- readPNG("figures/cerrada_con_pc1_restime.png")
fs2b <- rasterGrob(FS2B, interpolate = TRUE)
FS2F <- readPNG("figures/abierta_pc1_restime.png")
fs2f <- rasterGrob(FS2F, interpolate = TRUE)
FS2K <- readPNG("figures/abierta_pc2_restime.png")
fs2k <- rasterGrob(FS2K, interpolate = TRUE)

FS3A <- readPNG("figures/fig1/rotada-tiempo0.png")
fs3a <- rasterGrob(FS3A, interpolate = TRUE)
FS3B <- readPNG("figures/fig1/rotada-tiempo400.png")
fs3b <- rasterGrob(FS3B, interpolate = TRUE)

FS5D <- readPNG("figures/fig 5 parte 2C.png")
fs5d <- rasterGrob(FS5D, interpolate = TRUE)
FS5G <- readPNG("figures/rot-pc2.png")
fs5g <- rasterGrob(FS5G, interpolate = TRUE)

FS6A <- readPNG("figures/dna/abierta-chaina3.png")
fs6a <- rasterGrob(FS6A, interpolate = TRUE)
FS6B <- readPNG("figures/dna/abierta-chainb3.png")
fs6b <- rasterGrob(FS6B, interpolate = TRUE)
###
################################################ 2- FIGURA 1 ####################################################
## 4 PANELES
# PANEL A --> APO-STING tiempo = 100 ns
# PANEL B --> APO-STING tiempo > 100 ns
# PANEL C --> HOLO-STING tiempo = 100 ns
# PANEL D --> HOLO-STING tiempo > 100 ns

A1 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f1a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) 
B1 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f1b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) 
C1 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f1c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
      axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
      plot.title = element_text(colour = 'black', size =15)) 
D1 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f1d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
      axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
      plot.title = element_text(colour = 'black', size =15)) 

png("figures/new/F1.png", height = 122.7, width = 184, units = "mm", res = 300)
ggarrange(A1, C1, B1, D1,  labels = c("A", "B", "C", "D"),
          common.legend = TRUE, legend = "bottom",
          font.label = list(size = 7/.35), heights = c(1,1),
          widths = c(1,1)) 
dev.off()
##
################################################ 3- FIGURA 2 ####################################################
## 6 PANELES
# PANEL A --> B-factors APO&HOLO-STING chain A
# PANEL B --> B-factors APO&HOLO-STING chain B
# PANEL C --> B-factors APO-STING structure
# PANEL D --> B-factors HOLO-STING structure
# PANEL E --> B-factors APO-STING structure zoom polymerization domain
# PANEL F --> B-factors HOLO-STING structure zoom polymerization domain

A2 <- ggplot(bfac[c(1:197, 395:591),], aes(x = residue, y = nbfac)) + ylim(-0.6, 10.15) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  scale_linetype_manual(values=c("twodash", "solid")) + a +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax=290, ymin= -0.6, ymax= 3, fill=NA, colour="black", size = .4) +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
inset_a2 <- ggplot(bfac[c(131:145, 525:539),], aes(x = residue, y = nbfac)) + ylim(-0.6, 2.5) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  scale_linetype_manual(values=c("twodash", "solid"))  +
  scale_x_continuous(breaks = c(276, 278, 280, 282, 284, 286, 288, 290),
                     labels = paste(c("276", "278", "280", "282", "284", "286", "288", "290"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none") 
A2 <- A2 + annotation_custom(ggplotGrob(inset_a2), xmin=245, xmax=315, ymin=4, ymax=10) + 
  annotate("rect", xmin = 245, xmax=315, ymin= 4, ymax= 10, fill=NA, colour="black", size = .4)
B2 <- ggplot(bfac[c(198:394, 592:788),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  scale_linetype_manual(values=c("twodash", "solid")) + a + ylim(-0.6, 10.15) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.6, ymax=10.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax=487, ymin=-0.6, ymax= 3, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
inset_b2 <- ggplot(bfac[c(328:342, 722:736),], aes(x = residue, y = nbfac)) +  ylim(-0.6, 2.5) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  scale_linetype_manual(values=c("twodash", "solid")) + ylab("Normalized B-factor") +
  scale_x_continuous(breaks = c(473, 475, 477, 479, 481, 483, 485, 487),
                     labels = paste(c("276", "278", "280", "282", "284", "286", "288", "290"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )
B2 <- B2 + annotation_custom(ggplotGrob(inset_b2), xmin=442, xmax=512, ymin=4, ymax=10) +
  annotate("rect", xmin = 442, xmax=512, ymin= 4, ymax= 10, fill=NA, colour="black", size = .4)  
C2 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f2c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
D2 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f2d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
E2 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f2e, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
F2 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f2f  , xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
dev.off()
png("figures/new/F2.png", height = 184, width = 170, units = "mm", res = 300)
ggarrange(ggarrange(A2, B2, labels = c("A", "B"), label.x = -0.015, 
                    ncol = 2, nrow = 1, common.legend = TRUE,
                    legend = "bottom", font.label = list(size = 7/.35)),
          ggarrange(C2, D2, E2, F2, labels = c("C", "D", "E", "F"),  label.x = -0.015,
                    ncol = 2, nrow = 2, heights = c(2,1), font.label = list(size = 7/.35)),
          ncol = 1, nrow = 2, heights = c(1, 1.5))
dev.off()
##
################################################ 4- FIGURA 3 ####################################################
## 4 PANELES
# PANEL A --> contribucion de residuos 2CP HOLO-STING chain A
# PANEL B --> contribucion de residuos 2CP HOLO-STING chain B
# PANEL C --> trayectoria 2CP HOLO-STING
# PANEL D --> trayectoria 2CP HOLO-STING zoom polymerization domain
A3 <- ggplot(contribution_6nt7_holo[1:197,], aes(x = contribution_6nt7_holo$num[1:197])) + 
  geom_line(aes(y = contribution_6nt7_holo$pc2[1:197]), color = "mediumorchid2", size = 1) +
  ylab("Contribution 2PC") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax = 288,ymin=0, ymax=.25, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt7_holo_pc2, ymax=0.25, alpha=.1, fill="red")
B3 <- ggplot(contribution_6nt7_holo[198:394,], aes(x = contribution_6nt7_holo$num[198:394])) + 
  geom_line(aes(y = contribution_6nt7_holo$pc2[198:394]), color = "mediumorchid2", size = .8) +
  ylab("Contribution") + xlab("Residue") +ylab(" ") +
  geom_line(aes(y = threshold_6nt7_holo_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.25, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax = 485,ymin=0, ymax=.25, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt7_holo_pc2, ymax=0.25, alpha=.1, fill="red") + a 

C3 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f3c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
D3 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f3d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

png("figures/new/F3.png", height = 122.7, width = 184, units = "mm", res = 300)
ggarrange(A3, B3, C3, D3, labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2, label.x = -0.015,
          font.label = list(size = 7/0.35))
dev.off()
##
################################################ 5- FIGURA 4 ####################################################
## 3 PANELES
# PANEL A --> esquema de puentes de hidrogeno
# PANEL B --> estructura con el esquema de puentes de hidrogeno
# PANEL C --> puente de hidrogeno entre cGAMP y SER268A

A4 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f4a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
B4 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f4b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
C4 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f4c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) +
  annotate("text", x=-Inf, y=2, label="Chain A", hjust=-.2, vjust=2, fontface="bold")
D4 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f4d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) + 
  annotate("text", x=-Inf, y=2, label="Chain B", hjust=-.2, vjust=2, fontface="bold")
dev.off()
png("figures/new/F4.png", height = 184, width = 184, units = "mm", res = 300)
ggarrange(A4, B4, ggarrange(C4, D4, labels = c("C", "D"), font.label = list(size = 7/.35), nrow = 1, ncol = 2),
          labels = c("A", "B", ""), font.label = list(size = 7/.35), nrow = 3, ncol = 1)
dev.off()
       ##
################################################ 6- FIGURA 5 ####################################################
## 6 PANELES
# PANEL A --> Alineamiento de los dos ligandos. 
# PANEL B --> trayectoria 1CP rot-STING
# PANEL C --> B-factors HOLO&rot-STING chain A
# PANEL D --> B-factors HOLO&rot-STING chain B
# PANEL E --> Contribuciones 1CP rot-STING chain A
# PANEL F --> Contribuciones 1CP rot-STING chain B

A5 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f5a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
A5bis <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f5abis, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
B5 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f5b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

C5 <-  ggplot(bfac[c(395:591, 789:985),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("solid", "twodash")) + 
  xlab("Residue") + ylab("Normalized B-factor        ") + ylim(-0.6, 10) +  a +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  annotate("rect", xmin = 189, xmax=200, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax=290, ymin= -0.6, ymax= 4.5, fill=NA, colour="black", size = .4) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.title.x = element_blank()
  ) 

inset_c5 <- ggplot(bfac[c(525:539, 919:933),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("solid", "twodash")) + 
  ylab("Normalized B-factor  ") +  ylim(-0.6, 4.3) + 
  scale_x_continuous(breaks = c(276, 278, 280, 282, 284, 286, 288, 290),
                     labels = paste(c("276", "278", "280", "282", "284", "286", "288", "290"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none") 
C5 <- C5 + annotation_custom(ggplotGrob(inset_c5), xmin=245, xmax=315, ymin=5, ymax=10) +
  annotate("rect", xmin = 245, xmax=315, ymin= 5, ymax= 10, fill=NA, colour="black", size = .4)

D5 <- ggplot(bfac[c(592:788, 986:1182),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("solid", "twodash")) + 
  xlab("Residue") + ylab(" ") + ylim(-0.6, 10) + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("rect", xmin = 473, xmax=487, ymin= -0.6, ymax= 4.5, fill=NA, colour="black", size = .4) +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.6, ymax=10, alpha=.2, fill="grey60") +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black"),
    axis.title.x = element_blank()
  ) 
inset_d5 <- ggplot(bfac[c(722:736, 1116:1130),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("solid", "twodash")) + 
  ylab("Normalized B-factor") +  ylim(-0.6, 4.3) + 
  scale_x_continuous(breaks = c(473, 475, 477, 479, 481, 483, 485, 487),
                     labels = paste(c("276", "278", "280", "282", "284", "286", "288", "290"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none") 
D5 <- D5 + annotation_custom(ggplotGrob(inset_d5), xmin=442, xmax=512, ymin=5, ymax=10) +
  annotate("rect", xmin = 442, xmax=512, ymin= 5, ymax= 10, fill=NA, colour="black", size = .4)
E5 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f5e, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
dev.off()
png("figures/new/F5.png", height = 184, width = 184, units = "mm", res = 300)

ggarrange(ggarrange(A5, A5bis, B5, labels = c("A", "", " "), 
                    font.label = list(size = 7/.35), nrow = 1, ncol = 3, 
                    label.x = -0.015, common.legend = TRUE, legend = "bottom"),
          ggarrange(C5, D5, labels = c("C", "D"),
                    font.label = list(size = 7/.35), nrow = 1, ncol = 2, 
                    label.x = -0.015, common.legend = TRUE, legend = "bottom"),
          E5, labels = c("", "", "E"), nrow = 3, ncol = 1, label.x = -0.0075,
          heights = c(1,1.3,1), font.label = list(size = 7/.35))
dev.off()

################################################ 7- FIGURA 6 ####################################################
# Definir bien
A6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f6a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
B6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f6b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
C6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f6c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
D6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(f6d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
png("figures/new/F6.png", height = 122.7, width = 184, units = "mm", res = 300)

ggarrange(A6, B6, C6, D6, labels = c("A", "B", "C", "D"),
          font.label = list(size = 7/.35), nrow = 2, ncol = 2)
dev.off()
############################################### 8- FIGURA S1  ####################################################
## 6 PANELES 
# PANEL A --> Estructura PDB 6nt6
# PANEL B --> APO-STING a tiempo = 0
# PANEL C --> Estructura PDN 6nt7
# PANEL D --> HOLO-STING a tiempo = 0
# PANEL E --> RMSD APO&HOLO-STING
# PANEL F --> Radio de giro APO&HOLO-STING 

A1S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs1a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15),
        plot.margin = unit(c(1, 4, 1, 4), "pt"))
B1S <- ggplot(rms[0:5002,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a
C1S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs1c, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15),
        plot.margin = unit(c(1, 4, 1, 4), "pt")) 
D1S <-  ggplot(gr[0:5002,], aes(x = time, y = gr)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  scale_linetype_manual(values=c("dashed", "solid")) + a
png("figures/new/FS1.png", height = 122.7, width = 184, units = "mm", res = 300)

ggarrange(A1S, C1S, ggarrange(B1S, D1S, labels = c("C", "D"),
                              nrow = 2, ncol = 1, font.label = list(size = 7/.35),
                              common.legend = TRUE, legend = "bottom", label.x = -0.015),
          labels = c("A", "B", ""), font.label = list(size = 7/.35),
          nrow = 1, ncol = 3, widths = c(1, 1, 1.5))
dev.off()

##
############################################### 9- FIGURA 22 ####################################################
# Analisis de componentes principales, definir bien

A2S <- ggplot(eigen_bd_6nt6[c(1:7),], aes(x = eigen_bd_6nt6$num[c(1:7)], y=eigen_bd_6nt6$percent_6nt6[c(1:7)])) + 
  geom_col(fill="mediumseagreen") + 
  geom_text(aes(y= eigen_bd_6nt6$percent_6nt6[c(1:7)] + 2, 
                label = round(eigen_bd_6nt6$sum_percent[c(1:7)], 1)), 
            color = "black", size = 3, hjust=0.5) +
  xlab("Principal component") +
  ylab("Acumulated\npercentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,10,20,30,40,50)) + 
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA))
B2S <- ggplot(eigen_bd_6nt7_holo[c(1:7),], aes(x = eigen_bd_6nt7_holo$num[c(1:7)], 
                                        y=eigen_bd_6nt7_holo$percent[c(1:7)])) + 
  geom_col(fill="mediumorchid2") + 
  geom_text(aes(y= eigen_bd_6nt7_holo$percent[c(1:7)] + 2, 
                label = round(eigen_bd_6nt7_holo$sum_percent[c(1:7)], 1)), 
            color = "black", size = 3, hjust=0.5) +
  xlab("Principal Component") + 
  a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,10,20,30,40,50)) +
  ylab(" ") + theme(panel.grid.major = element_line(color =NA), 
                    panel.grid.minor = element_line(color =NA)) 
  
C2S <- ggplot(contribution_6nt6[1:197,], aes(x = contribution_6nt6$num[1:197])) + 
  geom_line(aes(y = contribution_6nt6$pc1[1:197]), color = "mediumseagreen", size = 0.8) +
  ylab("Contribution 1PC   ") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_6nt6_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt6_pc1, ymax=0.35, alpha=.1, fill="red")
D2S <- ggplot(contribution_6nt6[198:394,], aes(x = contribution_6nt6$num[198:394])) + 
  geom_line(aes(y = contribution_6nt6$pc1[198:394]), color = "mediumseagreen", size = 0.8) +
  ylab(" ") + xlab(" ") + 
  geom_line(aes(y = threshold_6nt6_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 363, 383, 403, 423, 443, 463, 483, 503, 523),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt6_pc1, ymax=0.35, alpha=.1, fill="red") + a
E2S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs2f, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

F2S <- ggplot(contribution_6nt6[1:197,], aes(x = contribution_6nt6$num[1:197])) + 
  geom_line(aes(y = contribution_6nt6$pc2[1:197]), color = "mediumseagreen", size = 0.8) +
  ylab("Contribution 2PC   ") + xlab(" ") + a + 
  geom_line(aes(y = threshold_6nt6_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt6_pc2, ymax=0.35, alpha=.1, fill="red")
G2S <- ggplot(contribution_6nt6[198:394,], aes(x = contribution_6nt6$num[198:394])) + 
  geom_line(aes(y = contribution_6nt6$pc2[198:394]), color = "mediumseagreen", size = 0.8) +
  ylab(" ") + xlab(" ") + 
  geom_line(aes(y = threshold_6nt6_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 363, 383, 403, 423, 443, 463, 483, 503, 523),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt6_pc2, ymax=0.35, alpha=.1, fill="red") + a
H2S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs2k, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

I2S <- ggplot(contribution_6nt7_holo[1:197,], aes(x = contribution_6nt7_holo$num[1:197])) + 
  geom_line(aes(y = contribution_6nt7_holo$pc1[1:197]), color = "mediumorchid2", size = 0.8) +
  ylab("Contribution 1CP   ") + xlab(" ") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt7_holo_pc1, ymax=0.35, alpha=.1, fill="red")
J2S <- ggplot(contribution_6nt7_holo[198:394,], aes(x = contribution_6nt7_holo$num[198:394])) + 
  geom_line(aes(y = contribution_6nt7_holo$pc1[198:394]), color = "mediumorchid2", size = 0.8) +
  ylab(" ") + xlab(" ") + 
  geom_line(aes(y = threshold_6nt7_holo_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 363, 383, 403, 423, 443, 463, 483, 503, 523),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt7_holo_pc1, ymax=0.35, alpha=.1, fill="red") + a
K2S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs2b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

png("figures/new/FS2.png", height = 184, width = 184, units = "mm", res = 300)

ggarrange(ggarrange(A2S, B2S, nrow = 1, ncol = 2),
          ggarrange(C2S, D2S, E2S, F2S, G2S, H2S, I2S, J2S, K2S,
                    nrow = 3, ncol = 3),
          nrow = 2, ncol = 1, heights = c(1,3)) + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm"))
dev.off()

##
############################################### 10- FIGURA S3 ###################################################
## 4 PANELES
# PANEL A --> rot-STING a tiempo = 0
# PANEL B --> rot-STING a tiempo > 0
# PANEL C --> RMSD rot-STING
# PANEL D --> Radio de giro rot-STING 

A3S <-  qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs3a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) 
B3S <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs3b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15)) 
C3S <- ggplot(rms[2502:7503,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("solid", "twodash")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a
D3S <- ggplot(gr[2502:7503,], aes(x = time, y = gr)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) +
  scale_linetype_manual(values=c("solid", "twodash")) + a

png("figures/new/FS3.png", height = 122.7, width = 184, units = "mm", res = 300)
ggarrange(A3S,C3S, B3S, D3S, labels = c("A", "B", "C", "D"),
          common.legend = TRUE, legend = "bottom", label.x = -0.015,
          font.label = list(size = 7/.35)) 
dev.off()
##
############################################### 11- FIGURA S5 ###################################################
## 4 PANELES
# PANEL A --> Contribucion PCA rot-STING
# PANEL B --> Trayectorias 2CP rot-STING
# PANEL C --> Contribucion 2CP rot-STING chain A
# PANEL D --> Contribucion 2CP rot-STING chain B

AS5 <-  ggplot(eigen_bd_6nt7_holo_rot[c(1:7),], aes(x = eigen_bd_6nt7_holo_rot$num[c(1:7)], 
                                                y=eigen_bd_6nt7_holo_rot$percent[c(1:7)])) + 
  geom_col(fill="dodgerblue3") + 
  geom_text(aes(y= eigen_bd_6nt7_holo_rot$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_holo_rot$sum_percent[c(1:7)], 2)), 
            color = "black", size = 3, hjust=0.5) +
  xlab("Principal Component") + 
  ylab("Acumulated\npercentaje of variance") +
  a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,10,20,30)) + 
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA)) + ylim(0, 25)

BS5 <- ggplot(contribution_6nt7_holo_rot[1:197,], aes(x = contribution_6nt7_holo_rot$num[1:197])) + 
  geom_line(aes(y = contribution_6nt7_holo_rot$pc1[1:197]), color = "dodgerblue3", size = 1) +
  ylab("Contribution 1PC") + xlab("Residues") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7)) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 273, xmax = 290,ymin=0, ymax=.2, fill = NA, color = "black", size=.4 ) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt7_holo_rot_pc1, ymax=0.3, alpha=.1, fill="red") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"))
  
CS5 <- ggplot(contribution_6nt7_holo_rot[198:394,], aes(x = contribution_6nt7_holo_rot$num[198:394])) + 
  geom_line(aes(y = contribution_6nt7_holo_rot$pc1[198:394]), color = "dodgerblue3", size = 1) +
  ylab(" ") + xlab("Residues") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(342, 534) + 
  scale_x_continuous(breaks = c(343, 363, 383, 403, 423, 443, 463, 483, 503, 523),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7)) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 470, xmax = 483, ymin=0, ymax=.2, fill = NA, color = "black", size=.4 ) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt7_holo_rot_pc2, ymax=0.3, alpha=.1, fill="red") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) 
  
DS5 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs5d, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

ES5 <- ggplot(contribution_6nt7_holo_rot[1:197,], aes(x = contribution_6nt7_holo_rot$num[1:197])) + 
  geom_line(aes(y = contribution_6nt7_holo_rot$pc2[1:197]), color = "dodgerblue3", size = 1) +
  ylab("Contribution 2PC") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7)) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 273, xmax = 286,ymin=0, ymax=.25, fill = NA, color = "black", size=.4 ) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_6nt7_holo_rot_pc2, ymax=0.3, alpha=.1, fill="red") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) 
  
FS5 <-  ggplot(contribution_6nt7_holo_rot[198:394,], aes(x = contribution_6nt7_holo_rot$num[198:394])) + 
  geom_line(aes(y = contribution_6nt7_holo_rot$pc2[198:394]), color = "dodgerblue3", size = 1) +
  ylab(" ") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(342, 534) + 
  scale_x_continuous(breaks = c(343, 363, 383, 403, 423, 443, 463, 483, 503, 523),
                     labels = paste(c("ALA146", "TRP166", "GLU186", "LEU206", "TYR226", 
                                      "LYS246", "PHE266", "ARG286", "SER306", "GLU326"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7)) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 470, xmax = 483, ymin=0, ymax=.25, fill = NA, color = "black", size=.4 ) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_6nt7_holo_rot_pc2, ymax=0.3, alpha=.1, fill="red") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) 
GS5 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs5g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))


png("figures/new/FS5.png", height = 122.7, width = 184, units = "mm", res = 300)
ggarrange(AS5, ggarrange(BS5, CS5, DS5, ES5, FS5, GS5,nrow = 2, ncol = 3),
          nrow = 2, ncol = 1, heights = c(1,2)) + theme(plot.margin = margin(0.4,0.1,0.1,0.1, "cm"))
dev.off()
############################################### 12- FIGURA S6 ###################################################
## 2 PANELES
# PANEL A --> DNA APO-STING chain A
# PANEL B --> DNA APO-STING chain B

AS6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs6a, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

BS6 <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(fs6b, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

png("figures/new/FS6.png", height = 92, width = 184, units = "mm", res = 300)
ggarrange(AS6, BS6, labels = c("A", "B"), font.label = list(size = 7/.35))
dev.off()
  ##
#################################################################################################################