setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

rg <- read.csv("~/Desktop/Sysbio/cSTING/THR272VAL/Prod_500ns/cerrada-t272v-500_gr.xvg", sep="")
rg$time <- rg$time/1000
gr_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_gr.xvg", sep="")
gr_6nt7_holo$Time = gr_6nt7_holo$Time/1000
gr <- data.frame("time" = c(gr_6nt7_holo$Time, rg$time),
                 "gr" = c(gr_6nt7_holo$GR, rg$all), 
                 "Structure" = as.character("label"))
gr$Structure[0:2501] <- "WT HOLO-STING"
gr$Structure[2502:5002] <- "T272V HOLO-STING"
rm(rg, gr_6nt7_holo)
gr$Structure <- factor(gr$Structure, levels = c("WT HOLO-STING", "T272V HOLO-STING"))

rmsd <- read.csv("~/Desktop/Sysbio/cSTING/THR272VAL/Prod_500ns/cerrada-t272v-500_rmsd.xvg", sep="")
rms_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="Time"] <- "Time"
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="RMSDns"] <- "RMSD"
rms <- data.frame("time" = c(rms_6nt7_holo$Time, rmsd$time),
                  "rmsd" = c(rms_6nt7_holo$RMSD, rmsd$rmsd),
                  "Structure" = as.character("label"))
rms$Structure[0:2501] <- "WT HOLO-STING"
rms$Structure[2502:5002] <- "T272V HOLO-STING"
rm(rmsd, rms_6nt7_holo)
rms$Structure <- factor(rms$Structure, levels = c("WT HOLO-STING", "T272V HOLO-STING"))

bfac_6nt7_holo <- read.pdb("6nt7_HOLO_NO_BORRAR/400ns_6nt7_holo_NpT_nojump_rot_bfac.pdb")
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
bfac <- read.pdb("~/Desktop/Sysbio/cSTING/THR272VAL/Prod_500ns/cerrada-t272v-500_bfac.pdb")
bfac_db <- data.frame("NUM" = bfac$atom$resno, 
                           "RES" = bfac$atom$resid, 
                           "BFAC" = bfac$atom$b)
prom <- mean(bfac_db$BFAC)
sd <- sd(bfac_db$BFAC)
bfac_db$NBFAC <- (bfac_db$BFAC - prom)/sd
j = 146
for (i in 1:394){
  bfac_db$ATOM[i] = j
  j = j+1
}
bfacs <- data.frame("residue" = c(bfac_6nt7_holo_db$ATOM, bfac_db$ATOM), 
                   "nbfac" = c(bfac_6nt7_holo_db$NBFAC, bfac_db$NBFAC),
                   "Structure" = "label")
bfacs$Structure <- as.character(bfacs$Structure)
bfacs$Structure[1:394] <- "WT HOLO-STING"
bfacs$Structure[395:788] <- "T272V HOLO-STING"
rm(bfac, bfac_6nt7_holo, bfac_6nt7_holo_db, bfac_db, i, j, prom, prom_6nt7_holo, sd, sd_6nt7_holo)
bfacs$Structure <- factor(bfacs$Structure, levels = c("WT HOLO-STING", "T272V HOLO-STING"))


ggplot(rms, aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumorchid2", "firebrick1")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a

ggplot(gr, aes(x = time, y = gr)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumorchid2", "firebrick1")) +
  scale_linetype_manual(values=c("dashed", "solid")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + ylim(2,2.5) + a

ggplot(bfacs[c(1:197, 395:591),], aes(x = residue, y = nbfac)) + ylim(-0.6, 10.15) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("mediumorchid2", "chartreuse4")) +
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
ggplot(bfacs[c(198:394, 592:788),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylab(" ") + xlab("Residues") + 
  scale_color_manual(values = c("mediumorchid2", "chartreuse4")) +
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
  annotate("rect", xmin = 473, xmax=487, ymin=-0.6, ymax= 7, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank()) 

ggplot(bfacs[c(131:145, 525:539),], aes(x = residue, y = nbfac)) + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylim(-0.5, 6) + 
  scale_color_manual(values = c("mediumorchid2", "chartreuse4")) +
  scale_linetype_manual(values=c("twodash", "solid"))  +
  scale_x_continuous(breaks = c(276, 277, 278, 279, 280, 281, 282, 283, 
                                284, 285, 286, 287,  288, 289,  290),
                     labels = paste(c("276MET", "277SER", "278GLN", "279ASP", 
                                      "280ASP", "281CYS", "282ALA", "283ALA", 
                                      "284PHE", "285SER", "286ARG", "287GLU", 
                                      "288GLN", "289ARG", "290LEU"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
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


ggplot(bfacs[c(328:342, 722:736),], aes(x = residue, y = nbfac))  + 
  geom_line(aes(color = Structure, linetype = Structure), size = 1) + ylim(-0.5, 6) + 
  scale_color_manual(values = c("mediumorchid2", "chartreuse4")) +
  scale_linetype_manual(values=c("twodash", "solid")) + 
  scale_x_continuous(breaks = c(473, 474, 475, 476, 477, 478, 479,
                                480, 481, 482, 483, 484, 485, 486, 487),
                     labels = paste(c("276MET", "277SER", "278GLN", "279ASP", 
                                      "280ASP", "281CYS", "282ALA", "283ALA", 
                                      "284PHE", "285SER", "286ARG", "287GLU", 
                                      "288GLN", "289ARG", "290LEU"))) +
  geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
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

hbonds <- read.csv("6nt7_HOLO_NO_BORRAR/holo_num_hbond.xvg", sep="")
hbonds[2502,] = hbonds[2501,]
hbonds <- hbonds[-c(1:500),]

pdb_t272v <- read.pdb("THR272VAL/Prod_500ns/cerrada_T272V.pdb")
dcd_t272v <- read.dcd("THR272VAL/Prod_500ns/cerrada_T272V.dcd")
ca_ndx_t272v <- atom.select(pdb_t272v, elety = "CA") 
pca_t272v <- pca.xyz(dcd_t272v[,ca_ndx_t272v$xyz])
plot(pca_t272v, col=topo.colors(nrow(dcd_t272v)))

pca_t272v_db <- data.frame("PC" = pca_t272v$z)
pca_t272v_db$time <- c(1:2002)
pca_t272v_db$time[(1:2002)] = hbonds$TIME
eigen_bd_t272v <- data.frame("eigen" = pca_t272v$L)
sum_eig_t272v = sum(eigen_bd_t272v$eigen)
eigen_bd_t272v$percent_t272v <- (eigen_bd_t272v$eigen/sum_eig_t272v)*100
eigen_bd_t272v$num <- c(1:1182)
eigen_bd_t272v$sum_percent = 0
eigen_bd_t272v$sum_percent[1] = eigen_bd_t272v$percent_t272v[1]
for (i in 2:1182){
  eigen_bd_t272v$sum_percent[i] = eigen_bd_t272v$percent_t272v[i] + eigen_bd_t272v$sum_percent[i - 1]
}
contribution_t272v <- data.frame("num" = pdb_t272v$atom$resno[ca_ndx_t272v$atom],
                                "resnum" = pdb_t272v$atom$resno[ca_ndx_t272v$atom],
                                "res" = pdb_t272v$atom$resid[ca_ndx_t272v$atom], 
                                "pc1" = pca_t272v$au[,1], 
                                "pc2" = pca_t272v$au[,2], 
                                "pc3" = pca_t272v$au[,3]) 
i = 146
for (j in 1:394){
  contribution_t272v$num[j] = i
  i = i+1
}
threshold_t272v_pc1 = (max(contribution_t272v$pc1) + min(contribution_t272v$pc1))/2
threshold_t272v_pc2 = (max(contribution_t272v$pc2) + min(contribution_t272v$pc2))/2
rm(ca_ndx_t272v, dcd_t272v, pca_t272v, pdb_t272v, i, j, sum_eig_t272v)
contribution_t272v <- contribution_t272v[-c(198:201, 399:402),]
rm(i, j)

ggplot(eigen_bd_t272v[c(1:7),], aes(x = eigen_bd_t272v$num[c(1:7)], y=eigen_bd_t272v$percent_t272v[c(1:7)])) + 
  geom_col(fill="chartreuse4") + 
  geom_text(aes(y= eigen_bd_t272v$percent_t272v[c(1:7)] + 1, 
                label = round(eigen_bd_t272v$sum_percent[c(1:7)], 1)), 
            color = "black", size = 3, hjust=0.5) +
  xlab("Principal component") +
  ylab("Acumulated\npercentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,10,20,30,40,50)) + 
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA))


ggplot(contribution_t272v[1:197,], aes(x = contribution_t272v$num[1:197])) + 
  geom_line(aes(y = contribution_t272v$pc1[1:197]), color = "chartreuse4", size = 1) +
  ylab("Contribution 1PC") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_t272v_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax = 288,ymin=0, ymax=.35, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_t272v_pc1, ymax=0.35, alpha=.1, fill="red")
ggplot(contribution_t272v[198:394,], aes(x = contribution_t272v$num[198:394])) + 
  geom_line(aes(y = contribution_t272v$pc1[198:394]), color = "chartreuse4", size = .8) +
  ylab("Contribution") + xlab("Residue") +ylab(" ") +
  geom_line(aes(y = threshold_t272v_pc1), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax = 485,ymin=0, ymax=.35, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_t272v_pc1, ymax=0.35, alpha=.1, fill="red") + a



ggplot(contribution_t272v[1:197,], aes(x = contribution_t272v$num[1:197])) + 
  geom_line(aes(y = contribution_t272v$pc2[1:197]), color = "chartreuse4", size = 1) +
  ylab("Contribution 2PC (5.5%)") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_t272v_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(145, 343) + 
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax = 288,ymin=0, ymax=.35, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 146, xmax=342, ymin= threshold_t272v_pc2, ymax=0.35, alpha=.1, fill="red")
ggplot(contribution_t272v[198:394,], aes(x = contribution_t272v$num[198:394])) + 
  geom_line(aes(y = contribution_t272v$pc2[198:394]), color = "chartreuse4", size = .8) +
  ylab("Contribution") + xlab("Residue") +ylab(" ") +
  geom_line(aes(y = threshold_t272v_pc2), color = "red", alpha = 0.5, linetype = 1) + xlim(343, 539) + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= 0, ymax=0.35, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax = 485,ymin=0, ymax=.35, fill = NA, color = "black", size=.4 )+
  theme(axis.text.x = element_text(angle = 90, face = "bold"),
        panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95"),
        axis.title.x = element_blank()) +
  annotate("rect", xmin = 343, xmax=539, ymin= threshold_t272v_pc2, ymax=0.35, alpha=.1, fill="red") + a 
