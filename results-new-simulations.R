setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic', size = 15),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)
##
############################################# RMSD  #############################################
open1 <- read.csv("abierta-rep1/crop-200/abierta-rep1-rmsd.xvg", sep="")
open2 <- read.csv("abierta-rep2/abierta-rep2-rmsd.xvg", sep="")
#open3
close1 <- read.csv("cerrada-rep1/prod/control_rmsd.xvg", sep="")
close2 <- read.csv("cerrada-rep2/cerrada-rep2-rmsd.xvg", sep="")
close3 <- read.csv("cerrada-rep3/cerrada-rep3-rmsd.xvg", sep="")
rot1 <- read.csv("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-rmsd.xvg", sep="")
rot2 <- read.csv("cerrada-ligrot-rep2/cerrada-ligrot-rep2-rmsd.xvg", sep="")
#rot3
dgmp1 <- read.csv("cerrada+cdiGMP-rep1/prod/cerrada+cdigmp-rmsd.xvg", sep="")
dgmp2 <- read.csv("cerrada+cdiGMP-rep2/cerrada+cdigmp-rep2-rmsd.xvg", sep="")
#dgmp3
s33cgamp1 <- read.csv("cerrada+33cGAMP-rep1/prod_200/cerrada+33cgamp-200ns-rmsd.xvg", sep="")
s33cgamp2 <- read.csv("cerrada+33cGAMP-rep2/cerrada+33cgamp-rep2-rmsd.xvg", sep="")
#s33cgamp3
v160m1 <- read.csv("v160m-rep1/prod_200ns/abierta_v160m_200ns_rmsd.xvg", sep="")
v160m2 <- read.csv("v160m-rep2/abierta-v160m-rep2-rmsd.xvg", sep="")
#v160m3 
sin1 <- read.csv("cerrada-sinlig-rep1/prod_200/cerrada-sinlig-rep1-rmsd.xvg", sep="")
sin2 <- read.csv("cerrada-sinlig-rep2/cerrada-sinlig-rep2-rmsd.xvg", sep="")

rmsd <- data.frame("time" = c(open1$time, open2$time, 
                              close1$time, close2$time, close3$time,
                              rot1$time, rot2$time,
                              s33cgamp1$time, s33cgamp2$time,
                              dgmp1$time, dgmp2$time,
                              v160m1$time, v160m2$time,
                              sin1$time, sin2$time),
                   "rmsd" = c(open1$rmsd, open2$rmsd,
                              close1$rmsd, close2$rmsd, close3$rmsd,
                              rot1$rmsd, rot2$rmsd,
                              s33cgamp1$rmsd, s33cgamp2$rmsd,
                              dgmp1$rmsd, dgmp2$rmsd,
                              v160m1$rmsd, v160m2$rmsd,
                              sin1$rmsd, sin2$rmsd),
                   "label" = "label")

rmsd$label[1:1001] <- "Open STING rep 1" #palegreen
rmsd$label[1002:2002] <- "Open STING rep 2" #limegreen
rmsd$label[2003:3003] <- "Close STING + 2'3'cGAMP rep 1" #blue
rmsd$label[3004:4004] <- "Close STING + 2'3'cGAMP rep 2" #dodgerblue
rmsd$label[4005:5005] <- "Close STING + 2'3'cGAMP rep 3" #deepskyblue
rmsd$label[5006:6006] <- "Close STING + rot 2'3'cGAMP rep 1" #black
rmsd$label[6007:7007] <- "Close STING + rot 2'3'cGAMP rep 2" #slategray4
rmsd$label[7008:8008] <- "Close STING + 3'3' cGAMP rep1" #red
rmsd$label[8009:9009] <- "Close STING + 3'3' cGAMP rep2" #firebrick1
rmsd$label[9010:10010] <- "Close STING + cdiGMP rep 1" #darkorchid2
rmsd$label[10011:11011] <- "Close STING + cdiGMP rep 2" #magenta4
rmsd$label[11012:12012] <-  "Open STING V160M rep 1" #deeppink
rmsd$label[12013:13013] <-  "Open STING V160M rep 2" #hotpink1
rmsd$label[13014:14014] <-  "Close STING w/o ligand rep 1" #chocolate1
rmsd$label[14015:15015] <-  "Close STING w/o ligand rep 2" #gold

  corridas <- c("Open STING rep 1", "Open STING rep 2",
                "Close STING + 2'3'cGAMP rep 1", "Close STING + 2'3'cGAMP rep 2", "Close STING + 2'3'cGAMP rep 3",
                "Close STING + rot 2'3'cGAMP rep 1", "Close STING + rot 2'3'cGAMP rep 2",
                "Close STING + 3'3' cGAMP rep1", "Close STING + 3'3' cGAMP rep2",
                "Close STING + cdiGMP rep 1", "Close STING + cdiGMP rep 2",
                "Open STING V160M rep 1", "Open STING V160M rep 2",
                "Close STING w/o ligand rep 1", "Close STING w/o ligand rep 2")
  
    rmsd$label <- as.factor(rmsd$label)
    rmsd$label <- ordered(rmsd$label, levels=corridas)
rm(open1, open2, open3, close1, close2, close3, rot1, rot2, rot3, dgmp1, dgmp2, dgmp3, 
   s33cgamp1, s33cgamp2, s33cgamp3, v160m1, v160m2, v160m3, sin1, sin2)

  rmsds <- ggplot(rmsd, aes(x = time, y = rmsd)) + 
    geom_line(aes(color = label, linetype = label), size = 1) + 
    scale_color_manual(values = c("palegreen", "limegreen",
                                    "blue", "dodgerblue", "deepskyblue", 
                                  "black", "slategray4", 
                                  "red", "firebrick1", 
                                  "darkorchid2", "magenta4",
                                  "deeppink", "hotpink1",
                                  "chocolate1", "gold")) +
    scale_linetype_manual(values=c("solid", "dashed",
                                   "solid", "dashed", "twodash", 
                                   "solid", "dashed",
                                   "solid", "dashed",
                                   "solid", "dashed",
                                   "solid", "dashed",
                                   "solid", "dashed")) +
    xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("pics-19Sep/rmsds.png", height = 300, width = 750, units = "mm", res = 300)
rmsds
dev.off()

rmsds_open <- ggplot(rmsd[1:2002,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("palegreen", "limegreen")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("pics-19Sep/rmsds_open.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_open
dev.off()

rmsds_close <- ggplot(rmsd[2003:5005,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("blue", "dodgerblue", "deepskyblue")) +
  scale_linetype_manual(values=c("solid", "dashed", "twodash")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("pics-19Sep/rmsds_close.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_close
dev.off()

rmsds_rot <- ggplot(rmsd[5006:7007,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("black", "slategray4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")

png("pics-19Sep/rmsds_rot.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_rot
dev.off()

rmsds_33 <- ggplot(rmsd[7008:9009,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("red", "firebrick1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("pics-19Sep/rmsds_33cgamps.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_33
dev.off()

rmsds_digmp <- ggplot(rmsd[9010:11011,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("darkorchid2", "magenta4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")
png("pics-19Sep/rmsds_digmp.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_digmp
dev.off()

rmsds_v160m <- ggplot(rmsd[11012:13013,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("deeppink", "hotpink1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")

png("pics-19Sep/rmsds_v160m.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_v160m
dev.off()

rmsds_sin<- ggplot(rmsd[13014:15015,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("chocolate1", "gold")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation")

png("pics-19Sep/rmsds_sing.png", height = 300, width = 750, units = "mm", res = 300)
rmsds_sin
dev.off()

########################################## B-FACTORS ###########################################
open1 <- read.pdb("abierta-rep1/crop-200/abierta-rep1-bfac.pdb")
open1_db <- data.frame("NUM" = open1$atom$resno, 
                       "RES" = open1$atom$resid, 
                       "BFAC" = open1$atom$b)
prom1 <- mean(open1_db$BFAC)
sd1 <- sd(open1_db$BFAC)
open1_db$NBFAC <- (open1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  open1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

open2 <- read.pdb("abierta-rep2/abierta-rep2-bfac.pdb")
open2_db <- data.frame("NUM" = open2$atom$resno, 
                       "RES" = open2$atom$resid, 
                       "BFAC" = open2$atom$b)
prom1 <- mean(open2_db$BFAC)
sd1 <- sd(open2_db$BFAC)
open2_db$NBFAC <- (open2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  open2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

#open3
close1 <- read.pdb("cerrada-rep1/prod/control_bfac.pdb")
close1_db <- data.frame("NUM" = close1$atom$resno, 
                        "RES" = close1$atom$resid, 
                        "BFAC" = close1$atom$b)
prom1 <- mean(close1_db$BFAC)
sd1 <- sd(close1_db$BFAC)
close1_db$NBFAC <- (close1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  close1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

close2 <- read.pdb("cerrada-rep2/cerrada-rep2-bfac.pdb")
close2_db <- data.frame("NUM" = close2$atom$resno, 
                        "RES" = close2$atom$resid, 
                        "BFAC" = close2$atom$b)
prom1 <- mean(close2_db$BFAC)
sd1 <- sd(close2_db$BFAC)
close2_db$NBFAC <- (close2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  close2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

close3 <- read.pdb("cerrada-rep3/cerrada-rep3-bfac.pdb")
close3_db <- data.frame("NUM" = close3$atom$resno, 
                        "RES" = close3$atom$resid, 
                        "BFAC" = close3$atom$b)
prom1 <- mean(close3_db$BFAC)
sd1 <- sd(close3_db$BFAC)
close3_db$NBFAC <- (close3_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  close3_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

rot1 <- read.pdb("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-bfac.pdb")
rot1_db <- data.frame("NUM" = rot1$atom$resno, 
                      "RES" = rot1$atom$resid, 
                      "BFAC" = rot1$atom$b)
prom1 <- mean(rot1_db$BFAC)
sd1 <- sd(rot1_db$BFAC)
rot1_db$NBFAC <- (rot1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  rot1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

rot2 <- read.pdb("cerrada-ligrot-rep2/cerrada-ligrot-rep2-bfac.pdb")
rot2_db <- data.frame("NUM" = rot2$atom$resno, 
                      "RES" = rot2$atom$resid, 
                      "BFAC" = rot2$atom$b)
prom1 <- mean(rot2_db$BFAC)
sd1 <- sd(rot2_db$BFAC)
rot2_db$NBFAC <- (rot2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  rot2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

#rot3
dgmp1 <- read.pdb("cerrada+cdiGMP-rep1/prod/cerrada+cdigmp-bfac.pdb")
dgmp1_db <- data.frame("NUM" = dgmp1$atom$resno, 
                       "RES" = dgmp1$atom$resid, 
                       "BFAC" = dgmp1$atom$b)
prom1 <- mean(dgmp1_db$BFAC)
sd1 <- sd(dgmp1_db$BFAC)
dgmp1_db$NBFAC <- (dgmp1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  dgmp1_db$ATOM[i] = j
  j = j+1
}

rm(prom1, sd1, i, j)

#dgmp2
dgmp2 <- read.pdb("cerrada+cdiGMP-rep2/cerrada+cdigmp-rep2-bfac.pdb")
dgmp2_db <- data.frame("NUM" = dgmp2$atom$resno, 
                       "RES" = dgmp2$atom$resid, 
                       "BFAC" = dgmp2$atom$b)
prom1 <- mean(dgmp2_db$BFAC)
sd1 <- sd(dgmp2_db$BFAC)
dgmp2_db$NBFAC <- (dgmp2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  dgmp2_db$ATOM[i] = j
  j = j+1
}

rm(prom1, sd1, i, j)

#dgmp3
s33cgamp1 <- read.pdb("cerrada+33cGAMP-rep1/prod_200/cerrada+33cgamp-200ns-bfac.pdb")
s33cgamp1_db <- data.frame("NUM" = s33cgamp1$atom$resno, 
                           "RES" = s33cgamp1$atom$resid, 
                           "BFAC" = s33cgamp1$atom$b)
prom1 <- mean(s33cgamp1_db$BFAC)
sd1 <- sd(s33cgamp1_db$BFAC)
s33cgamp1_db$NBFAC <- (s33cgamp1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  s33cgamp1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

s33cgamp2 <- read.pdb("cerrada+33cGAMP-rep2/cerrada+33cgamp-rep2-bfac.pdb")
s33cgamp2_db <- data.frame("NUM" = s33cgamp2$atom$resno, 
                           "RES" = s33cgamp2$atom$resid, 
                           "BFAC" = s33cgamp2$atom$b)
prom1 <- mean(s33cgamp2_db$BFAC)
sd1 <- sd(s33cgamp2_db$BFAC)
s33cgamp2_db$NBFAC <- (s33cgamp2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  s33cgamp2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)
#s33cgamp3

v160m1 <- read.pdb("v160m-rep1/prod_200ns/abierta_v160m_200ns_bfac.pdb")
v160m1_db <- data.frame("NUM" = v160m1$atom$resno, 
                           "RES" = v160m1$atom$resid, 
                           "BFAC" = v160m1$atom$b)

which(v160m1_db$ATOM == 346)
v160m1_db <- v160m1_db[-c(198:201, 399:402),]
prom1 <- mean(v160m1_db$BFAC)
sd1 <- sd(v160m1_db$BFAC)
v160m1_db$NBFAC <- (v160m1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  v160m1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

#v160m2
v160m2 <- read.pdb("v160m-rep2/abierta-v160m-rep2-bfac.pdb")
v160m2_db <- data.frame("NUM" = v160m2$atom$resno, 
                        "RES" = v160m2$atom$resid, 
                        "BFAC" = v160m2$atom$b)

which(v160m2_db$ATOM == 346)
v160m2_db <- v160m2_db[-c(198:201, 399:402),]
prom1 <- mean(v160m2_db$BFAC)
sd1 <- sd(v160m2_db$BFAC)
v160m2_db$NBFAC <- (v160m2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  v160m2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

#v160m3

sin1 <- read.pdb("cerrada-sinlig-rep1/prod_200/cerrada-sinlig-rep1-bfac.pdb")
sin1_db <- data.frame("NUM" = sin1$atom$resno, 
                           "RES" = sin1$atom$resid, 
                           "BFAC" = sin1$atom$b)
prom1 <- mean(sin1_db$BFAC)
sd1 <- sd(sin1_db$BFAC)
sin1_db$NBFAC <- (sin1_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  sin1_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

sin2 <- read.pdb("cerrada-sinlig-rep2/cerrada-sinlig-rep2-bfac.pdb")
sin2_db <- data.frame("NUM" = sin2$atom$resno, 
                           "RES" = sin2$atom$resid, 
                           "BFAC" = sin2$atom$b)
prom1 <- mean(sin2_db$BFAC)
sd1 <- sd(sin2_db$BFAC)
sin2_db$NBFAC <- (sin2_db$BFAC - prom1)/sd1
j = 146
for (i in 1:394){
  sin2_db$ATOM[i] = j
  j = j+1
}
rm(prom1, sd1, i, j)

bfacs <- data.frame("res" = c(open1_db$ATOM, open2_db$ATOM,
                              close1_db$ATOM, close2_db$ATOM, close3_db$ATOM,
                              rot1_db$ATOM, rot2_db$ATOM, 
                              s33cgamp1_db$ATOM, s33cgamp2_db$ATOM,
                              dgmp1_db$ATOM, dgmp2_db$ATOM, 
                              v160m1_db$ATOM, v160m2_db$ATOM,
                              sin1_db$ATOM, sin2_db$ATOM),
                    "nbfac" = c(open1_db$NBFAC, open2_db$NBFAC,
                                close1_db$NBFAC, close2_db$NBFAC, close3_db$NBFAC,
                                rot1_db$NBFAC, rot2_db$NBFAC, 
                                s33cgamp1_db$NBFAC, s33cgamp2_db$NBFAC,
                                dgmp1_db$NBFAC, dgmp2_db$NBFAC,
                                v160m1_db$NBFAC, v160m2_db$NBFAC,
                                sin1_db$NBFAC, sin2_db$NBFAC),
                    "label" = "label")

bfacs$label[1:394] <- "Open STING rep 1" #palegreen
bfacs$label[395:788] <- "Open STING rep 2" #limegreen
bfacs$label[789:1182] <- "Close STING + 2'3'cGAMP rep 1" #blue
bfacs$label[1183:1576] <- "Close STING + 2'3'cGAMP rep 2" #dodgerblue
bfacs$label[1577:1970] <- "Close STING + 2'3'cGAMP rep 3" #deepskyblue
bfacs$label[1971:2364] <- "Close STING + rot 2'3'cGAMP rep 1" #black
bfacs$label[2365:2758] <- "Close STING + rot 2'3'cGAMP rep 2" #slategray4
bfacs$label[2759:3152] <- "Close STING + 3'3' cGAMP rep1" #red
bfacs$label[3153:3546] <- "Close STING + 3'3' cGAMP rep2" #firebrick1 
bfacs$label[3547:3940] <- "Close STING + cdiGMP rep 1" #darkorchid2
bfacs$label[3941:4334] <- "Close STING + cdiGMP rep 2" #magenta4
bfacs$label[4335:4728] <- "Open STING V160M rep 1" # deeppink
bfacs$label[4729:5122] <- "Open STING V160M rep 2" # hotpink1
bfacs$label[5123:5516] <- "Close STING w/o ligand rep 1" # deeppink
bfacs$label[5517:5910] <- "Close STING w/o ligand rep 2" # hotpink1

bfacs$label <- as.factor(bfacs$label)
bfacs$label <- ordered(bfacs$label, levels=corridas)
rm(open1, open2, open3, close1, close2, close3, rot1, rot2, rot3, dgmp1, dgmp2, dgmp3, 
   s33cgamp1, s33cgamp2, s33cgamp3, v160m1, v160m2, sin1, sin2)
rm(open1_db, open2_db, open3_db, close1_db, close2_db, close3_db, rot1_db, rot2_db, rot3_db, 
   dgmp1_db, dgmp2_db, dgmp3_db, s33cgamp1_db, s33cgamp2_db, s33cgamp3_db, v160m1_db, v160m2_db, sin1_db, sin2_db)

#### all
bfacs_A <- ggplot(bfacs[c(1:197, 395:591, 789:985, 1183:1379),]
                  , aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("palegreen", "limegreen",
                                "blue", "dodgerblue")) +
  scale_linetype_manual(values=c("solid", "dashed",
                                 "solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_B <- ggplot(bfacs[c(198:394, 592:788, 986:1182, 1380:1576),], 
                  aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("palegreen", "limegreen",
                                "blue", "dodgerblue")) +
  scale_linetype_manual(values=c("solid", "dashed",
                                 "solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_all <- ggarrange(bfacs_A, bfacs_B,  labels = c("A", "B"),
                       common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_all
dev.off()
#### open
bfacs_opens_A <- ggplot(bfacs[c(1:197, 395:591),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("palegreen", "limegreen")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_opens_B <- ggplot(bfacs[c(198:394, 592:788),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("palegreen", "limegreen")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_open <- ggarrange(bfacs_opens_A, bfacs_opens_B,  labels = c("A", "B"),
                        common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_open.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_open
dev.off()
#### closeds
bfacs_close_A <- ggplot(bfacs[c(789:985, 1183:1379),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("blue", "dodgerblue")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_close_B <- ggplot(bfacs[c(986:1182, 1380:1576),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("blue", "dodgerblue")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_close <- ggarrange(bfacs_close_A, bfacs_close_B,  labels = c("A", "B"),
                         common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_close.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_close
dev.off()
#### rots
bfacs_rots_A <- ggplot(bfacs[c(1971:2167, 2365:2561),]
                       , aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("black", "slategray4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_rots_B <- ggplot(bfacs[c(2168:2364, 2562:2758),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("black", "slategray4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_rot <- ggarrange(bfacs_rots_A, bfacs_rots_B,  labels = c("A", "B"),
                       common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_rot.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_rot
dev.off()
#### 33cgamps
bfacs_33_A <- ggplot(bfacs[c(2759:2955, 3153:3349),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("red", "firebrick1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_33_B <- ggplot(bfacs[c(2956:3152, 3350:3546),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("red", "firebrick1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_33 <- ggarrange(bfacs_33_A, bfacs_33_B,  labels = c("A", "B"),
                      common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_33cgamps.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_33
dev.off()

#### digmp
bfacs_digmp_A <- ggplot(bfacs[c(3547:3743, 3941:4137),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("darkorchid2", "magenta4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_digmp_B <- ggplot(bfacs[c(3744:3940, 4138:4334),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("darkorchid2", "magenta4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_digmp <- ggarrange(bfacs_digmp_A, bfacs_digmp_B,  labels = c("A", "B"),
                         common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_digmp.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_digmp
dev.off()

#v160m
bfacs_v160m_A <- ggplot(bfacs[c(3941:4138, 4729:4925),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("deeppink", "hotpink1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_v160m_B <- ggplot(bfacs[c(4138:4333,4926:5122),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("deeppink", "hotpink1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_v160m <- ggarrange(bfacs_v160m_A, bfacs_v160m_B,  labels = c("A", "B"),
                         common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_v160m.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_v160m
dev.off()

## sin ligando
bfacs_sin_A <- ggplot(bfacs[c(5123:5319, 5517:5713),], aes(x = res, y = nbfac)) + ylim(-0.8, 12.4) + a +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("chocolate1", "gold")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 278, xmax=287, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_sin_B <- ggplot(bfacs[c(5320:5516, 5714:5910),], aes(x = res, y = nbfac)) + a + ylim(-0.8, 12.4) +
  geom_line(aes(color = label, linetype = label), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("chocolate1", "gold")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) + 
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 475, xmax=484, ymin= -0.8, ymax=12.4, alpha=.2, fill="yellowgreen") +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

bfacs_sin <- ggarrange(bfacs_sin_A, bfacs_sin_B,  labels = c("A", "B"),
                         common.legend = TRUE, legend = "bottom") 
png("pics-19Sep/bfacs_sin.png", height = 300, width = 750, units = "mm", res = 300)
bfacs_sin
dev.off()

###################################### RADIUS OF GYRATION ######################################
open1 <- read.csv("abierta-rep1/crop-200/abierta-rep1-gr.xvg", sep="")
open2 <- read.csv("abierta-rep2/abierta-rep2-gr.xvg", sep="")
#open3
close1 <- read.csv("cerrada-rep1/prod/control_rg.xvg", sep="")
close2 <- read.csv("cerrada-rep2/cerrada-rep2-rg.xvg", sep="")
close3 <- read.csv("cerrada-rep3/cerrada-rep3-rg.xvg", sep="")
rot1 <- read.csv("cerrada-ligrot-rep1/prod/cerrada+ligrot-rep1-gr.xvg", sep="")
rot2 <- read.csv("cerrada-ligrot-rep2/cerrada-ligrot-rep2-gr.xvg", sep="")
#rot3
dgmp1 <- read.csv("cerrada+cdiGMP-rep1/prod/cerrada+cdigmp-gr.xvg", sep="")
dgmp2 <- read.csv("cerrada+cdiGMP-rep2/cerrada+cdigmp-rep2-gr.xvg", sep="")
#dgmp3
s33cgamp1 <-read.csv("cerrada+33cGAMP-rep1/prod_200/cerrada+33cgamp-200ns-rg.xvg", sep="")
s33cgamp2 <- read.csv("cerrada+33cGAMP-rep2/cerrada+33cgamp-rep2-gr.xvg", sep="")
#s33cgamp3
v160m1 <- read.csv("v160m-rep1/prod_200ns/abierta_v160m_200ns_rg.xvg", sep="")
v160m2 <- read.csv("v160m-rep2/abierta-v160m-rep2-gr.xvg", sep="")
#v160m3
sin1 <- read.csv("cerrada-sinlig-rep1/prod_200/cerrada-sinlig-rep1-gr.xvg", sep="")
sin2 <- read.csv("cerrada-sinlig-rep2/cerrada-sinlig-rep2-gr.xvg", sep="")

gr <- data.frame("time" = c(open1$time, open2$time,
                            close1$time, close2$time, close3$time,
                            rot1$time, rot2$time,
                            s33cgamp1$time, s33cgamp2$time,
                            dgmp1$time, dgmp2$time,
                            v160m1$time, v160m2$time,
                            sin1$time, sin2$time),
                 "gr" = c(open1$all, open2$all,
                          close1$all, close2$all, close3$all,
                          rot1$all, rot2$all,
                          s33cgamp1$all, s33cgamp2$all,
                          dgmp1$all, dgmp2$all,
                          v160m1$all, v160m2$all,
                          sin1$all, sin2$all),
                 "label" = "label")

gr$label[1:1001] <- "Open STING rep 1" #palegreen
gr$label[1002:2002] <- "Open STING rep 2" #limegreen
gr$label[2003:3003] <- "Close STING + 2'3'cGAMP rep 1" #blue
gr$label[3004:4004] <- "Close STING + 2'3'cGAMP rep 2" #dodgerblue
gr$label[4005:5005] <- "Close STING + 2'3'cGAMP rep 3" #deepskyblue
gr$label[5006:6006] <- "Close STING + rot 2'3'cGAMP rep 1" #black
gr$label[6007:7007] <- "Close STING + rot 2'3'cGAMP rep 2" #slategray4
gr$label[7008:8008] <- "Close STING + 3'3' cGAMP rep1" #red
gr$label[8009:9009] <- "Close STING + 3'3' cGAMP rep2" #firebrick1
gr$label[9010:10010] <- "Close STING + cdiGMP rep 1" #darkorchid2
gr$label[10011:11011] <-  "Close STING + cdiGMP rep 2" #magenta4
gr$label[11012:12012] <-  "Open STING V160M rep 1" #deeppink
gr$label[12013:13013] <-  "Open STING V160M rep 2" #hotpink1
gr$label[13014:14014] <- "Close STING w/o ligand rep 1" # deeppink
gr$label[14015:15015] <- "Close STING w/o ligand rep 2" # hotpink1
rm(open1, open2, open3, close1, close2, close3, rot1, rot2, rot3, dgmp1, dgmp2, dgmp3, 
   s33cgamp1, s33cgamp2, s33cgamp3, v160m1, v160m2, v160m3)

gr$label <- as.factor(gr$label)
gr$label <- ordered(gr$label, levels=corridas)

grs <- ggplot(gr, aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("palegreen", "limegreen",
                                "blue", "dodgerblue", "deepskyblue", 
                                "black", "slategray4", 
                                "red", "firebrick1", 
                                "darkorchid2", "magenta4",
                                "deeppink", "hotpink1",
                                "chocolate1", "gold")) +
  scale_linetype_manual(values=c("solid", "dashed",
                                 "solid", "dashed", "twodash", 
                                 "solid", "dashed",
                                 "solid", "dashed",
                                 "solid", "dashed",
                                 "solid", "dashed",
                                 "solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs.png", height = 300, width = 750, units = "mm", res = 300)
grs
dev.off()

grs_open <- ggplot(gr[1:2002,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("palegreen", "limegreen")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_open.png", height = 300, width = 750, units = "mm", res = 300)
grs_open
dev.off()

grs_close <- ggplot(gr[2003:5005,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("blue", "dodgerblue", "deepskyblue")) +
  scale_linetype_manual(values=c("solid", "dashed", "twodash")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_close.png", height = 300, width = 750, units = "mm", res = 300)
grs_close
dev.off()

grs_rot <- ggplot(gr[5006:7007,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("black", "slategray4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_rot.png", height = 300, width = 750, units = "mm", res = 300)
grs_rot
dev.off()

grs_33cgamp <- ggplot(gr[7008:9009,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("red", "firebrick1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_33cgamps.png", height = 300, width = 750, units = "mm", res = 300)
grs_33cgamp
dev.off()

grs_digmp <- ggplot(gr[9010:11011,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("darkorchid2", "magenta4")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_digmp.png", height = 300, width = 750, units = "mm", res = 300)
grs_digmp
dev.off()

grs_v160m <- ggplot(gr[11012:13013,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("deeppink", "hotpink1")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_v160m.png", height = 300, width = 750, units = "mm", res = 300)
grs_v160m
dev.off()

grs_sin <- ggplot(gr[13014:15015,], aes(x = time, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = label, linetype = label), size = 1) + 
  scale_color_manual(values = c("chocolate1", "gold")) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + ggtitle("Radius of gyration")
png("pics-19Sep/grs_sin", height = 300, width = 750, units = "mm", res = 300)
grs_sin
dev.off()
########################################## EXPOSITION ##########################################
open1_chainA <- read.delim("abierta-rep1/crop-200/distance/abierta1_chainA.xvg")
open1_chainB <- read.delim("abierta-rep1/crop-200/distance/abierta1_chainB.xvg")
open1_protein <- read.delim("abierta-rep1/crop-200/distance/abierta1_protein.xvg")
open1 <- data.frame("time" = open1_chainA$time/1000, 
                    "chainA_x" = open1_chainA$x,"chainA_y" = open1_chainA$y, "chainA_z" = open1_chainA$z,
                    "chainB_x" = open1_chainB$x, "chainB_y" = open1_chainB$y,"chainB_z" = open1_chainB$z, 
                    "protein_x" = open1_protein$x,"protein_y" = open1_protein$y, "protein_z" = open1_protein$z)
open1$chainA_dist <- sqrt((open1$chainA_x - open1$protein_x)^2 + 
                            (open1$chainA_y - open1$protein_y)^2 + 
                            (open1$chainA_z - open1$protein_z)^2)
open1$chainB_dist <- sqrt((open1$chainB_x - open1$protein_x)^2 + 
                            (open1$chainB_y - open1$protein_y)^2 + 
                            (open1$chainB_z - open1$protein_z)^2)
open1$chainA_norm <- open1$chainA_dist - open1$chainA_dist[1]
open1$chainB_norm <- open1$chainB_dist - open1$chainB_dist[1]
rm(open1_chainA, open1_chainB, open1_protein)

open2_chainA <-read.delim("abierta-rep2/distance/abierta2_chainA.xvg")
open2_chainB <- read.delim("abierta-rep2/distance/abierta2_chainB.xvg")
open2_protein <- read.delim("abierta-rep2/distance/abierta2_protein.xvg")
open2 <- data.frame("time" = open2_chainA$time/1000, 
                    "chainA_x" = open2_chainA$x,"chainA_y" = open2_chainA$y, "chainA_z" = open2_chainA$z,
                    "chainB_x" = open2_chainB$x, "chainB_y" = open2_chainB$y,"chainB_z" = open2_chainB$z, 
                    "protein_x" = open2_protein$x,"protein_y" = open2_protein$y, "protein_z" = open2_protein$z)
open2$chainA_dist <- sqrt((open2$chainA_x - open2$protein_x)^2 + 
                            (open2$chainA_y - open2$protein_y)^2 + 
                            (open2$chainA_z - open2$protein_z)^2)
open2$chainB_dist <- sqrt((open2$chainB_x - open2$protein_x)^2 + 
                            (open2$chainB_y - open2$protein_y)^2 + 
                            (open2$chainB_z - open2$protein_z)^2)
open2$chainA_norm <- open2$chainA_dist - open2$chainA_dist[1]
open2$chainB_norm <- open2$chainB_dist - open2$chainB_dist[1]
rm(open2_chainA, open2_chainB, open2_protein)


dist_opens <- data.frame("time" = c(open1$time, open1$time, open2$time, open2$time),
                         "distance" = c(open1$chainA_norm, open1$chainB_norm, open2$chainA_norm, open2$chainB_norm),
                         "label" = "label")
dist_opens$label[1:1001] <- "Open STING rep 1 - Chain A"
dist_opens$label[1002:2002] <- "Open STING rep 1 - Chain B"
dist_opens$label[2003:3003] <- "Open STING rep 2 - Chain A"
dist_opens$label[3004:4004] <- "Open STING rep 2 - Chain B"

dist_open1 <- ggplot(dist_opens[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("palegreen", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING rep 1") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")

dist_open2 <- ggplot(dist_opens[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("limegreen", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING rep 2") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")

png("pics-19Sep/dist_opens.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_open1, dist_open2, ncol = 1, nrow = 2)
dev.off()

#open3
#open3_chainA
#open3_chainB
#open3_protein
#open3 <- data.frame("time" = open3_chainA$time/1000, 
#                    "chainA_x" = open3_chainA$x,"chainA_y" = open3_chainA$y, "chainA_z" = open3_chainA$z,
#                    "chainB_x" = open3_chainB$x, "chainB_y" = open3_chainB$y,"chainB_z" = open3_chainB$z, 
#                    "protein_x" = open3_protein$x,"protein_y" = open3_protein$y, "protein_z" = open3_protein$z)
#open3$chainA_dist <- sqrt((open3$chainA_x - open3$protein_x)^2 + 
#                            (open3$chainA_y - open3$protein_y)^2 + 
#                            (open3$chainA_z - open3$protein_z)^2)
#open3$chainB_dist <- sqrt((open3$chainB_x - open3$protein_x)^2 + 
#                            (open3$chainB_y - open3$protein_y)^2 + 
#                            (open3$chainB_z - open3$protein_z)^2)
#open3$chainA_norm <- open3$chainA_dist - open3$chainA_dist[1]
#open3$chainB_norm <- open3$chainB_dist - open3$chainB_dist[1]


close1_chainA <- read.delim("cerrada-rep1/prod/distances/chainA-rep1.xvg")
close1_chainB <- read.delim("cerrada-rep1/prod/distances/chainB-rep1.xvg")
close1_protein <- read.delim("cerrada-rep1/prod/distances/protein-rep1.xvg")
close1 <- data.frame("time" = close1_chainA$time/1000, 
                     "chainA_x" = close1_chainA$x,"chainA_y" = close1_chainA$y, "chainA_z" = close1_chainA$z,
                     "chainB_x" = close1_chainB$x, "chainB_y" = close1_chainB$y,"chainB_z" = close1_chainB$z, 
                     "protein_x" = close1_protein$x,"protein_y" = close1_protein$y, "protein_z" = close1_protein$z)
close1$chainA_dist <- sqrt((close1$chainA_x - close1$protein_x)^2 + 
                             (close1$chainA_y - close1$protein_y)^2 + 
                             (close1$chainA_z - close1$protein_z)^2)
close1$chainB_dist <- sqrt((close1$chainB_x - close1$protein_x)^2 + 
                             (close1$chainB_y - close1$protein_y)^2 + 
                             (close1$chainB_z - close1$protein_z)^2)
close1$chainA_norm <- close1$chainA_dist - close1$chainA_dist[1]
close1$chainB_norm <- close1$chainB_dist - close1$chainB_dist[1]
rm(close1_chainA, close1_chainB, close1_protein)

close2_chainA <- read.delim("cerrada-rep2/distances/chainA-rep2.xvg")
close2_chainB <- read.delim("cerrada-rep2/distances/chainB-rep2.xvg")
close2_protein <- read.delim("cerrada-rep2/distances/protein-rep2.xvg")
close2 <- data.frame("time" = close2_chainA$time/1000, 
                     "chainA_x" = close2_chainA$x,"chainA_y" = close2_chainA$y, "chainA_z" = close2_chainA$z,
                     "chainB_x" = close2_chainB$x, "chainB_y" = close2_chainB$y,"chainB_z" = close2_chainB$z, 
                     "protein_x" = close2_protein$x,"protein_y" = close2_protein$y, "protein_z" = close2_protein$z)
close2$chainA_dist <- sqrt((close2$chainA_x - close2$protein_x)^2 + 
                             (close2$chainA_y - close2$protein_y)^2 + 
                             (close2$chainA_z - close2$protein_z)^2)
close2$chainB_dist <- sqrt((close2$chainB_x - close2$protein_x)^2 + 
                             (close2$chainB_y - close2$protein_y)^2 + 
                             (close2$chainB_z - close2$protein_z)^2)
close2$chainA_norm <- close2$chainA_dist - close2$chainA_dist[1]
close2$chainB_norm <- close2$chainB_dist - close2$chainB_dist[1]
rm(close2_chainA, close2_chainB, close2_protein)

close3_chainA <- read.delim("cerrada-rep3/distances/chainA-rep3.xvg")
close3_chainB <- read.delim("cerrada-rep3/distances/chainB-rep3.xvg")
close3_protein <- read.delim("cerrada-rep3/distances/protein-rep3.xvg")

close3 <- data.frame("time" = close3_chainA$time/1000, 
                     "chainA_x" = close3_chainA$x,"chainA_y" = close3_chainA$y, "chainA_z" = close3_chainA$z,
                     "chainB_x" = close3_chainB$x, "chainB_y" = close3_chainB$y,"chainB_z" = close3_chainB$z, 
                     "protein_x" = close3_protein$x,"protein_y" = close3_protein$y, "protein_z" = close3_protein$z)
close3$chainA_dist <- sqrt((close3$chainA_x - close3$protein_x)^2 + 
                             (close3$chainA_y - close3$protein_y)^2 + 
                             (close3$chainA_z - close3$protein_z)^2)
close3$chainB_dist <- sqrt((close3$chainB_x - close3$protein_x)^2 + 
                             (close3$chainB_y - close3$protein_y)^2 + 
                             (close3$chainB_z - close3$protein_z)^2)
close3$chainA_norm <- close3$chainA_dist - close3$chainA_dist[1]
close3$chainB_norm <- close3$chainB_dist - close3$chainB_dist[1]
rm(close3_chainA, close3_chainB, close3_protein)

dist_closed <- data.frame("time" = c(close1$time, close1$time, close2$time, close2$time),
                          "distance" = c(close1$chainA_norm, close1$chainB_norm, 
                                         close2$chainA_norm, close2$chainB_norm),
                          "label" = "label")
dist_closed$label[1:1001] <- "Closed STING 2'3'cGAMP rep 1 - Chain A"
dist_closed$label[1002:2002] <- "Closed STING 2'3'cGAMP rep 1 - Chain B"
dist_closed$label[2003:3003] <- "Closed STING 2'3'cGAMP rep 2 - Chain A"
dist_closed$label[3004:4004] <- "Closed STING 2'3'cGAMP rep 2 - Chain B"
dist_closed$label[4005:5005] <- "Closed STING 2'3'cGAMP rep 3 - Chain A"
dist_closed$label[5006:6006] <- "Closed STING 2'3'cGAMP rep 3 - Chain B"

dist_close1 <- ggplot(dist_closed[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("blue", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 1") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")
dist_close2 <- ggplot(dist_closed[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("dodgerblue", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 2") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")
dist_close3 <- ggplot(dist_closed[4005:6006,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deepskyblue", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 3") 

png("pics-19Sep/dist_closed.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_close1, dist_close2, ncol = 1, nrow = 2)
dev.off()

  #rot1
rot1_chainA <- read.delim("cerrada-ligrot-rep1/prod/distances/chainA-rot-rep1.xvg")
rot1_chainB <- read.delim("cerrada-ligrot-rep1/prod/distances/chainB-rot-rep1.xvg")
rot1_protein <- read.delim("cerrada-ligrot-rep1/prod/distances/protein-rot-rep1.xvg")
rot1 <- data.frame("time" = rot1_chainA$time/1000, 
                   "chainA_x" = rot1_chainA$x,"chainA_y" = rot1_chainA$y, "chainA_z" = rot1_chainA$z,
                   "chainB_x" = rot1_chainB$x, "chainB_y" = rot1_chainB$y,"chainB_z" = rot1_chainB$z, 
                   "protein_x" = rot1_protein$x,"protein_y" = rot1_protein$y, "protein_z" = rot1_protein$z)
rot1$chainA_dist <- sqrt((rot1$chainA_x - rot1$protein_x)^2 + 
                           (rot1$chainA_y - rot1$protein_y)^2 + 
                           (rot1$chainA_z - rot1$protein_z)^2)
rot1$chainB_dist <- sqrt((rot1$chainB_x - rot1$protein_x)^2 + 
                           (rot1$chainB_y - rot1$protein_y)^2 + 
                           (rot1$chainB_z - rot1$protein_z)^2)
rot1$chainA_norm <- rot1$chainA_dist - rot1$chainA_dist[1]
rot1$chainB_norm <- rot1$chainB_dist - rot1$chainB_dist[1]
rm(rot1_chainA, rot1_chainB, rot1_protein)

#rot2
rot2_chainA <- read.delim("cerrada-ligrot-rep2/distance/chainA-rot-rep2.xvg")
rot2_chainB <- read.delim("cerrada-ligrot-rep2/distance/chainB-rot-rep2.xvg")
rot2_protein <- read.delim("cerrada-ligrot-rep2/distance/protein-rot-rep2.xvg")
rot2 <- data.frame("time" = rot2_chainA$time/1000, 
                   "chainA_x" = rot2_chainA$x,"chainA_y" = rot2_chainA$y, "chainA_z" = rot2_chainA$z,
                   "chainB_x" = rot2_chainB$x, "chainB_y" = rot2_chainB$y,"chainB_z" = rot2_chainB$z, 
                   "protein_x" = rot2_protein$x,"protein_y" = rot2_protein$y, "protein_z" = rot2_protein$z)
rot2$chainA_dist <- sqrt((rot2$chainA_x - rot2$protein_x)^2 + 
                           (rot2$chainA_y - rot2$protein_y)^2 + 
                           (rot2$chainA_z - rot2$protein_z)^2)
rot2$chainB_dist <- sqrt((rot2$chainB_x - rot2$protein_x)^2 + 
                           (rot2$chainB_y - rot2$protein_y)^2 + 
                           (rot2$chainB_z - rot2$protein_z)^2)
rot2$chainA_norm <- rot2$chainA_dist - rot2$chainA_dist[1]
rot2$chainB_norm <- rot2$chainB_dist - rot2$chainB_dist[1]
rm(rot2_chainA, rot2_chainB, rot2_protein)
#rot3
#rot3 <- data.frame("time" = rot3_chainA$time/1000, 
#                   "chainA_x" = rot3_chainA$x,"chainA_y" = rot3_chainA$y, "chainA_z" = rot3_chainA$z,
#                   "chainB_x" = rot3_chainB$x, "chainB_y" = rot3_chainB$y,"chainB_z" = rot3_chainB$z, 
#                   "protein_x" = rot3_protein$x,"protein_y" = rot3_protein$y, "protein_z" = rot3_protein$z)
#rot3$chainA_dist <- sqrt((rot3$chainA_x - rot3$protein_x)^2 + 
#                           (rot3$chainA_y - rot3$protein_y)^2 + 
#                           (rot3$chainA_z - rot3$protein_z)^2)
#rot3$chainB_dist <- sqrt((rot3$chainB_x - rot3$protein_x)^2 + 
#                           (rot3$chainB_y - rot3$protein_y)^2 + 
#                           (rot3$chainB_z - rot3$protein_z)^2)
#rot3$chainA_norm <- rot3$chainA_dist - rot3$chainA_dist[1]
#rot3$chainB_norm <- rot3$chainB_dist - rot3$chainB_dist[1]

dist_rots <- data.frame("time" = c(rot1$time, rot1$time, rot2$time, rot2$time),
                        "distance" = c(rot1$chainA_norm, rot1$chainB_norm, 
                                       rot2$chainA_norm, rot2$chainB_norm),
                        "label" = "label")
dist_rots$label[1:1001] <- "Closed STING rot 2'3'cGAMP rep 1 - Chain A"
dist_rots$label[1002:2002] <- "Closed STING rot 2'3'cGAMP rep 1 - Chain B"
dist_rots$label[2003:3003] <- "Closed STING rot 2'3'cGAMP  rep 2 - Chain A"
dist_rots$label[3004:4004] <- "Closed STING rot 2'3'cGAMP  rep 2 - Chain B"
#dist_closed$label[4005:5005] <- "Closed STING rep 3 - Chain A"
#dist_closed$label[5006:6006] <- "Closed STING rep 3 - Chain B"

dist_rot1 <- ggplot(dist_rots[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("black", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 1") 
dist_rot2 <- ggplot(dist_rots[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("slategray4", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 2") 
#dist_close3 <- ggplot(dist_closed[4005:6006,], aes(x = time, y = distance, color = label)) +
#  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
#  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
#  scale_color_manual(values = c("deepskyblue", "darkorange1")) + a + ylim(-0.5, 0.5) +
#  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 3") 

png("pics-19Sep/dist_rot.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_rot1, dist_rot2, ncol = 1, nrow = 2)
dev.off()

#dgmp1
dgmp1_chainA <- read.delim("cerrada+cdiGMP-rep1/prod/distance/cerrada+digmp-rep1-polA.xvg")
dgmp1_chainB <- read.delim("cerrada+cdiGMP-rep1/prod/distance/cerrada+digmp-rep1-polB.xvg")
dgmp1_protein <- read.delim("cerrada+cdiGMP-rep1/prod/distance/cerrada+digmp-rep1-protein.xvg")
dgmp2_chainA <- read.delim("cerrada+cdiGMP-rep2/distance/cerrada+digmp-rep2-polsiteA.xvg")
dgmp2_chainB <- read.delim("cerrada+cdiGMP-rep2/distance/cerrada+digmp-rep2-polsiteB.xvg")
dgmp2_protein <- read.delim("cerrada+cdiGMP-rep2/distance/cerrada+digmp-rep2-protein.xvg")
dgmp1 <- data.frame("time" = dgmp1_chainA$time/1000, 
                    "chainA_x" = dgmp1_chainA$x,"chainA_y" = dgmp1_chainA$y, "chainA_z" = dgmp1_chainA$z,
                    "chainB_x" = dgmp1_chainB$x, "chainB_y" = dgmp1_chainB$y,"chainB_z" = dgmp1_chainB$z, 
                    "protein_x" = dgmp1_protein$x,"protein_y" = dgmp1_protein$y, "protein_z" = dgmp1_protein$z)
dgmp1$chainA_dist <- sqrt((dgmp1$chainA_x - dgmp1$protein_x)^2 + 
                            (dgmp1$chainA_y - dgmp1$protein_y)^2 + 
                            (dgmp1$chainA_z - dgmp1$protein_z)^2)
dgmp1$chainB_dist <- sqrt((dgmp1$chainB_x - dgmp1$protein_x)^2 + 
                            (dgmp1$chainB_y - dgmp1$protein_y)^2 + 
                            (dgmp1$chainB_z - dgmp1$protein_z)^2)
dgmp1$chainA_norm <- dgmp1$chainA_dist - dgmp1$chainA_dist[1]
dgmp1$chainB_norm <- dgmp1$chainB_dist - dgmp1$chainB_dist[1]
rm(dgmp1_chainA, dgmp1_chainB, dgmp1_protein)
#dgmp2
dgmp2 <- data.frame("time" = dgmp2_chainA$time/1000, 
                    "chainA_x" = dgmp2_chainA$x,"chainA_y" = dgmp2_chainA$y, "chainA_z" = dgmp2_chainA$z,
                    "chainB_x" = dgmp2_chainB$x, "chainB_y" = dgmp2_chainB$y,"chainB_z" = dgmp2_chainB$z, 
                    "protein_x" = dgmp2_protein$x,"protein_y" = dgmp2_protein$y, "protein_z" = dgmp2_protein$z)
dgmp2$chainA_dist <- sqrt((dgmp2$chainA_x - dgmp2$protein_x)^2 + 
                            (dgmp2$chainA_y - dgmp2$protein_y)^2 + 
                            (dgmp2$chainA_z - dgmp2$protein_z)^2)
dgmp2$chainB_dist <- sqrt((dgmp2$chainB_x - dgmp2$protein_x)^2 + 
                            (dgmp2$chainB_y - dgmp2$protein_y)^2 + 
                            (dgmp2$chainB_z - dgmp2$protein_z)^2)
dgmp2$chainA_norm <- dgmp2$chainA_dist - dgmp2$chainA_dist[1]
dgmp2$chainB_norm <- dgmp2$chainB_dist - dgmp2$chainB_dist[1]

#dgmp3
#dgmp3 <- data.frame("time" = dgmp3_chainA$time/1000, 
#                    "chainA_x" = dgmp3_chainA$x,"chainA_y" = dgmp3_chainA$y, "chainA_z" = dgmp3_chainA$z,
#                    "chainB_x" = dgmp3_chainB$x, "chainB_y" = dgmp3_chainB$y,"chainB_z" = dgmp3_chainB$z, 
#                    "protein_x" = dgmp3_protein$x,"protein_y" = dgmp3_protein$y, "protein_z" = dgmp3_protein$z)
#dgmp3$chainA_dist <- sqrt((dgmp3$chainA_x - dgmp3$protein_x)^2 + 
#                            (dgmp3$chainA_y - dgmp3$protein_y)^2 + 
#                            (dgmp3$chainA_z - dgmp3$protein_z)^2)
#dgmp3$chainB_dist <- sqrt((dgmp3$chainB_x - dgmp3$protein_x)^2 + 
#                            (dgmp3$chainB_y - dgmp3$protein_y)^2 + 
#                            (dgmp3$chainB_z - dgmp3$protein_z)^2)
#dgmp3$chainA_norm <- dgmp3$chainA_dist - dgmp3$chainA_dist[1]
#dgmp3$chainB_norm <- dgmp3$chainB_dist - dgmp3$chainB_dist[1]

dist_digmps <- data.frame("time" = c(dgmp1$time, dgmp1$time, dgmp2$time, dgmp2$time),
                          "distance" = c(dgmp1$chainA_norm, dgmp1$chainB_norm, 
                                         dgmp2$chainA_norm, dgmp2$chainB_norm),
                          "label" = "label")
dist_digmps$label[1:1001] <- "Closed STING cdiGMP rep 1 - Chain A"
dist_digmps$label[1002:2002] <- "Closed STING cdiGMP rep 1 - Chain B"
dist_digmps$label[2003:3003] <- "Closed STING cdiGMP rep 2 - Chain A"
dist_digmps$label[3004:4004] <- "Closed STING cdiGMP rep 2 - Chain B"

dist_digmp1 <- ggplot(dist_digmps[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("darkorchid", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING + cdiGMP rep 1")

dist_digmp2 <- ggplot(dist_digmps[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("magenta4", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING + cdiGMP rep 1")

png("pics-19Sep/dist_cdiGMP.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_digmp1, dist_digmp2, ncol = 1, nrow = 2)
dev.off()

#s33cgamp1
s33cgamp1_chainA <- read.delim("cerrada+33cGAMP-rep1/prod_200/distances/chainA-33cgamp-rep1.xvg")
s33cgamp1_chainB <- read.delim("cerrada+33cGAMP-rep1/prod_200/distances/chainB-33cgamp-rep1.xvg")
s33cgamp1_protein <- read.delim("cerrada+33cGAMP-rep1/prod_200/distances/protein-33cgamp-rep1.xvg")
s33cgamp1 <- data.frame("time" = s33cgamp1_chainA$time/1000, 
                        "chainA_x" = s33cgamp1_chainA$x,"chainA_y" = s33cgamp1_chainA$y, "chainA_z" = s33cgamp1_chainA$z,
                        "chainB_x" = s33cgamp1_chainB$x, "chainB_y" = s33cgamp1_chainB$y,"chainB_z" = s33cgamp1_chainB$z, 
                        "protein_x" = s33cgamp1_protein$x,"protein_y" = s33cgamp1_protein$y, "protein_z" = s33cgamp1_protein$z)
s33cgamp1$chainA_dist <- sqrt((s33cgamp1$chainA_x - s33cgamp1$protein_x)^2 + 
                                (s33cgamp1$chainA_y - s33cgamp1$protein_y)^2 + 
                                (s33cgamp1$chainA_z - s33cgamp1$protein_z)^2)
s33cgamp1$chainB_dist <- sqrt((s33cgamp1$chainB_x - s33cgamp1$protein_x)^2 + 
                                (s33cgamp1$chainB_y - s33cgamp1$protein_y)^2 + 
                                (s33cgamp1$chainB_z - s33cgamp1$protein_z)^2)
s33cgamp1$chainA_norm <- s33cgamp1$chainA_dist - s33cgamp1$chainA_dist[1]
s33cgamp1$chainB_norm <- s33cgamp1$chainB_dist - s33cgamp1$chainB_dist[1]
rm(s33cgamp1_chainA, s33cgamp1_chainB, s33cgamp1_protein)

#s33cgamp2
s33cgamp2_chainA <- read.delim("cerrada+33cGAMP-rep2/distances/cerrada+33cgamp2-chainA.xvg")
s33cgamp2_chainB <- read.delim("cerrada+33cGAMP-rep2/distances/cerrada+33cgamp2-chainB.xvg")
s33cgamp2_protein <- read.delim("cerrada+33cGAMP-rep2/distances/cerrada+33cgamp2-protein.xvg")
s33cgamp2 <- data.frame("time" = s33cgamp2_chainA$time/1000, 
                        "chainA_x" = s33cgamp2_chainA$x,"chainA_y" = s33cgamp2_chainA$y, "chainA_z" = s33cgamp2_chainA$z,
                        "chainB_x" = s33cgamp2_chainB$x, "chainB_y" = s33cgamp2_chainB$y,"chainB_z" = s33cgamp2_chainB$z, 
                        "protein_x" = s33cgamp2_protein$x,"protein_y" = s33cgamp2_protein$y, "protein_z" = s33cgamp2_protein$z)
s33cgamp2$chainA_dist <- sqrt((s33cgamp2$chainA_x - s33cgamp2$protein_x)^2 + 
                                (s33cgamp2$chainA_y - s33cgamp2$protein_y)^2 + 
                                (s33cgamp2$chainA_z - s33cgamp2$protein_z)^2)
s33cgamp2$chainB_dist <- sqrt((s33cgamp2$chainB_x - s33cgamp2$protein_x)^2 + 
                                (s33cgamp2$chainB_y - s33cgamp2$protein_y)^2 + 
                                (s33cgamp2$chainB_z - s33cgamp2$protein_z)^2)
s33cgamp2$chainA_norm <- s33cgamp2$chainA_dist - s33cgamp2$chainA_dist[1]
s33cgamp2$chainB_norm <- s33cgamp2$chainB_dist - s33cgamp2$chainB_dist[1]
rm(s33cgamp2_chainA, s33cgamp2_chainB, s33cgamp2_protein)

#s33cgamp3
#s33cgamp3 <- data.frame("time" = s33cgamp3_chainA$time/1000, 
#                        "chainA_x" = s33cgamp3_chainA$x,"chainA_y" = s33cgamp3_chainA$y, "chainA_z" = s33cgamp3_chainA$z,
#                        "chainB_x" = s33cgamp3_chainB$x,"chainB_y" = s33cgamp3_chainB$y, "chainB_z" = s33cgamp3_chainB$z, 
#                        "protein_x" = s33cgamp3_protein$x,"protein_y" = s33cgamp3_protein$y, "protein_z" = s33cgamp3_protein$z)
#s33cgamp3$chainA_dist <- sqrt((s33cgamp3$chainA_x - s33cgamp3$protein_x)^2 + 
#                                (s33cgamp3$chainA_y - s33cgamp3$protein_y)^2 + 
#                                (s33cgamp3$chainA_z - s33cgamp3$protein_z)^2)
#s33cgamp3$chainB_dist <- sqrt((s33cgamp3$chainB_x - s33cgamp3$protein_x)^2 + 
#                                (s33cgamp3$chainB_y - s33cgamp3$protein_y)^2 + 
#                                (s33cgamp3$chainB_z - s33cgamp3$protein_z)^2)
#s33cgamp3$chainA_norm <- s33cgamp3$chainA_dist - s33cgamp3$chainA_dist[1]
#s33cgamp3$chainB_norm <- s33cgamp3$chainB_dist - s33cgamp3$chainB_dist[1]

dist_33cgamp <- data.frame("time" = c(s33cgamp1$time, s33cgamp1$time, s33cgamp2$time, s33cgamp2$time),
                           "distance" = c(s33cgamp1$chainA_norm, s33cgamp1$chainB_norm, 
                                          s33cgamp2$chainA_norm, s33cgamp2$chainB_norm),
                           "label" = "label")
dist_33cgamp$label[1:1001] <- "Closed STING 3'3'cGAMP rep 1 - Chain A"
dist_33cgamp$label[1002:2002] <- "Closed STING 3'3'cGAMP rep 1 - Chain B"
dist_33cgamp$label[2003:3003] <- "Closed STING 3'3'cGAMP  rep 2 - Chain A"
dist_33cgamp$label[3004:4004] <- "Closed STING 3'3'cGAMP  rep 2 - Chain B"
#dist_closed$label[4005:5005] <- "Closed STING rep 3 - Chain A"
#dist_closed$label[5006:6006] <- "Closed STING rep 3 - Chain B"

dist_33cgamp1 <- ggplot(dist_33cgamp[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("red", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING 3'3'cGAMP rep 1") 
dist_33cgamp2 <- ggplot(dist_33cgamp[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("firebrick1", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING 3'3'cGAMP rep 2") 
#dist_close3 <- ggplot(dist_closed[4005:6006,], aes(x = time, y = distance, color = label)) +
#  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
#  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
#  scale_color_manual(values = c("deepskyblue", "darkorange1")) + a + ylim(-0.5, 0.5) +
#  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 3") 

png("pics-19Sep/dist_33cgamp.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_33cgamp1, dist_33cgamp2, ncol = 1, nrow = 2)
dev.off()

rm(open1, open2, open3, close1, close2, close3, rot1, rot2, rot3, dgmp1, dgmp2, dgmp3, 
   s33cgamp1, s33cgamp2, s33cgamp3)

v160m1_chainA <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_polsiteA.xvg")
v160m1_chainB <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_polsiteB.xvg")
v160m1_protein <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_protein.xvg")
v160m1_A <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_chainA.xvg")
v160m1_B <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_chainB.xvg")
v160m1 <- data.frame("time" = v160m1_chainA$time/1000, 
                    "chainA_x" = v160m1_chainA$x,"chainA_y" = v160m1_chainA$y, "chainA_z" = v160m1_chainA$z,
                    "chainB_x" = v160m1_chainB$x, "chainB_y" = v160m1_chainB$y,"chainB_z" = v160m1_chainB$z, 
                    "protein_x" = v160m1_protein$x,"protein_y" = v160m1_protein$y, "protein_z" = v160m1_protein$z,
                    "A_x" = v160m1_A$x, "A_y" = v160m1_A$y, "A_z" = v160m1_A$z,
                    "B_x" = v160m1_B$x, "B_y" = v160m1_B$y, "B_z" = v160m1_B$z)
v160m1$chainA_dist <- sqrt((v160m1$chainA_x - v160m1$protein_x)^2 + 
                            (v160m1$chainA_y - v160m1$protein_y)^2 + 
                            (v160m1$chainA_z - v160m1$protein_z)^2)
v160m1$chainB_dist <- sqrt((v160m1$chainB_x - v160m1$protein_x)^2 + 
                            (v160m1$chainB_y - v160m1$protein_y)^2 + 
                            (v160m1$chainB_z - v160m1$protein_z)^2)
v160m1$chainA_norm <- v160m1$chainA_dist - v160m1$chainA_dist[1]
v160m1$chainB_norm <- v160m1$chainB_dist - v160m1$chainB_dist[1]

v160m1$A_dist <- sqrt((v160m1$chainA_x - v160m1$A_x)^2 + 
                             (v160m1$chainA_y - v160m1$A_y)^2 + 
                             (v160m1$chainA_z - v160m1$A_z)^2)
v160m1$B_dist <- sqrt((v160m1$chainB_x - v160m1$B_x)^2 + 
                             (v160m1$chainB_y - v160m1$B_y)^2 + 
                             (v160m1$chainB_z - v160m1$B_z)^2)
v160m1$chainA_norm <- v160m1$chainA_dist - v160m1$chainA_dist[1]
v160m1$chainB_norm <- v160m1$chainB_dist - v160m1$chainB_dist[1]
v160m1$A_norm <- v160m1$A_dist - v160m1$A_dist[1]
v160m1$B_norm <- v160m1$B_dist - v160m1$B_dist[1]
rm(v160m1_chainA, v160m1_chainB, v160m1_protein, v160m1_A, v160m1_B)

v160m2_chainA <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-polA.xvg")
v160m2_chainB <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-polB.xvg")
v160m2_protein <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-protein.xvg")
v160m2_A <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-chainA.xvg")
v160m2_B <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-chainB.xvg")

v160m2 <- data.frame("time" = v160m2_chainA$time/1000, 
                     "chainA_x" = v160m2_chainA$x,"chainA_y" = v160m2_chainA$y, "chainA_z" = v160m2_chainA$z,
                     "chainB_x" = v160m2_chainB$x, "chainB_y" = v160m2_chainB$y,"chainB_z" = v160m2_chainB$z, 
                     "protein_x" = v160m2_protein$x,"protein_y" = v160m2_protein$y, "protein_z" = v160m2_protein$z,
                     "A_x" = v160m2_A$x, "A_y" = v160m2_A$y, "A_z" = v160m2_A$z,
                     "B_x" = v160m2_B$x, "B_y" = v160m2_B$y, "B_z" = v160m2_B$z)
v160m2$chainA_dist <- sqrt((v160m2$chainA_x - v160m2$protein_x)^2 + 
                             (v160m2$chainA_y - v160m2$protein_y)^2 + 
                             (v160m2$chainA_z - v160m2$protein_z)^2)
v160m2$chainB_dist <- sqrt((v160m2$chainB_x - v160m2$protein_x)^2 + 
                             (v160m2$chainB_y - v160m2$protein_y)^2 + 
                             (v160m2$chainB_z - v160m2$protein_z)^2)
v160m2$A_dist <- sqrt((v160m2$chainA_x - v160m2$A_x)^2 + 
                             (v160m2$chainA_y - v160m2$A_y)^2 + 
                             (v160m2$chainA_z - v160m2$A_z)^2)
v160m2$B_dist <- sqrt((v160m2$chainB_x - v160m2$B_x)^2 + 
                             (v160m2$chainB_y - v160m2$B_y)^2 + 
                             (v160m2$chainB_z - v160m2$B_z)^2)
v160m2$chainA_norm <- v160m2$chainA_dist - v160m2$chainA_dist[1]
v160m2$chainB_norm <- v160m2$chainB_dist - v160m2$chainB_dist[1]
v160m2$A_norm <- v160m2$A_dist - v160m2$A_dist[1]
v160m2$B_norm <- v160m2$B_dist - v160m2$B_dist[1]
rm(v160m2_chainA, v160m2_chainB, v160m2_protein, v160m2_A, v160m2_B)

dist_v160m <- data.frame("time" = c(v160m1$time, v160m1$time, v160m2$time, v160m2$time),
                          "distance" = c(v160m1$chainA_norm, v160m1$chainB_norm,
                                         v160m2$chainA_norm, v160m2$chainB_norm,
                                         v160m1$A_norm, v160m1$B_norm,
                                         v160m2$A_norm, v160m2$B_norm),
                          "label" = "label")
dist_v160m$label[1:1001] <- "Open STING V160M rep 1 - Chain A"
dist_v160m$label[1002:2002] <- "Open STING V160M rep 1 - Chain B"
dist_v160m$label[2003:3003] <- "Open STING V160M rep 2 - Chain A"
dist_v160m$label[3004:4004] <- "Open STING V160M rep 2 - Chain B"
dist_v160m$label[4005:5005] <- "Open STING V160M rep 1 - Chain A : Polsite A"
dist_v160m$label[5006:6006] <- "Open STING V160M rep 1 - Chain B : Polsite B"
dist_v160m$label[6007:7007] <- "Open STING V160M rep 2 - Chain A : Polsite A"
dist_v160m$label[7008:8008] <- "Open STING V160M rep 2 - Chain B : Polsite B"

dist_v160m1 <- ggplot(dist_v160m[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 1") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")
dist_v160m2 <- ggplot(dist_v160m[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("hotpink1", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 2") +
  geom_ysidedensity(aes(x = after_stat(density)), position = "stack")

dist_v160m1_pol <- ggplot(dist_v160m[4005:6006,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 1") 
dist_v160m2_pol <- ggplot(dist_v160m[6007:8008,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("hotpink1", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 2") 

png("pics-19Sep/dist_v160m.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m1, dist_v160m2, ncol = 1, nrow = 2)
dev.off()

png("pics-19Sep/dist_v160m_pol-chain.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m1_pol, dist_v160m2_pol, ncol = 1, nrow = 2)
dev.off()


##Cerrada sin ligando
sin1_chainA <- read.delim("cerrada-sinlig-rep1/prod_200/distance/cerrada-sinlig-rep1-polA.xvg")
sin1_chainB <- read.delim("cerrada-sinlig-rep1/prod_200/distance/cerrada-sinlig-rep1-polB.xvg")
sin1_protein <- read.delim("cerrada-sinlig-rep1/prod_200/distance/cerrada-sinlig-rep1-protein.xvg")
sin1 <- data.frame("time" = sin1_chainA$time/1000, 
                   "chainA_x" = sin1_chainA$x,"chainA_y" = sin1_chainA$y, "chainA_z" = sin1_chainA$z,
                   "chainB_x" = sin1_chainB$x, "chainB_y" = sin1_chainB$y,"chainB_z" = sin1_chainB$z, 
                   "protein_x" = sin1_protein$x,"protein_y" = sin1_protein$y, "protein_z" = sin1_protein$z)
sin1$chainA_dist <- sqrt((sin1$chainA_x - sin1$protein_x)^2 + 
                           (sin1$chainA_y - sin1$protein_y)^2 + 
                           (sin1$chainA_z - sin1$protein_z)^2)
sin1$chainB_dist <- sqrt((sin1$chainB_x - sin1$protein_x)^2 + 
                           (sin1$chainB_y - sin1$protein_y)^2 + 
                           (sin1$chainB_z - sin1$protein_z)^2)
sin1$chainA_norm <- sin1$chainA_dist - sin1$chainA_dist[1]
sin1$chainB_norm <- sin1$chainB_dist - sin1$chainB_dist[1]
rm(sin1_chainA, sin1_chainB, sin1_protein)

sin2_chainA <-read.delim("cerrada-sinlig-rep2/distance/cerrada-sinlig-rep2-polsiteA.xvg")
sin2_chainB <- read.delim("cerrada-sinlig-rep2/distance/cerrada-sinlig-rep2-polsiteB.xvg")
sin2_protein <- read.delim("cerrada-sinlig-rep2/distance/cerrada-sinlig-rep2-protein.xvg")
sin2 <- data.frame("time" = sin2_chainA$time/1000, 
                   "chainA_x" = sin2_chainA$x,"chainA_y" = sin2_chainA$y, "chainA_z" = sin2_chainA$z,
                   "chainB_x" = sin2_chainB$x, "chainB_y" = sin2_chainB$y,"chainB_z" = sin2_chainB$z, 
                   "protein_x" = sin2_protein$x,"protein_y" = sin2_protein$y, "protein_z" = sin2_protein$z)
sin2$chainA_dist <- sqrt((sin2$chainA_x - sin2$protein_x)^2 + 
                           (sin2$chainA_y - sin2$protein_y)^2 + 
                           (sin2$chainA_z - sin2$protein_z)^2)
sin2$chainB_dist <- sqrt((sin2$chainB_x - sin2$protein_x)^2 + 
                           (sin2$chainB_y - sin2$protein_y)^2 + 
                           (sin2$chainB_z - sin2$protein_z)^2)
sin2$chainA_norm <- sin2$chainA_dist - sin2$chainA_dist[1]
sin2$chainB_norm <- sin2$chainB_dist - sin2$chainB_dist[1]
rm(sin2_chainA, sin2_chainB, sin2_protein)


dist_sins <- data.frame("time" = c(sin1$time, sin1$time, sin2$time, sin2$time),
                        "distance" = c(sin1$chainA_norm, sin1$chainB_norm, sin2$chainA_norm, sin2$chainB_norm),
                        "label" = "label")
dist_sins$label[1:1001] <- "Closed STING w/o ligand rep 1 - Chain A"
dist_sins$label[1002:2002] <- "Closed STING w/o ligand rep 1 - Chain B"
dist_sins$label[2003:3003] <- "Closed STING w/o ligand rep 2 - Chain A"
dist_sins$label[3004:4004] <- "Closed STING w/o ligand rep 2 - Chain B"

dist_sin1 <- ggplot(dist_sins[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("chocolate1", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING w/o ligand rep 1") 

dist_sin2 <- ggplot(dist_sins[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("gold", "black")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING w/o ligand rep 2") 

png("pics-19Sep/dist_sin.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(dist_sin1, dist_sin2, ncol = 1, nrow = 2)
dev.off()


