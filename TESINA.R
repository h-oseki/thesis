###########################################################################################################################
###########################################################################################################################
#####                                                                                                                 #####
#####                                                SCRIPT TESINA                                                    #####
#####                                                                                                                 #####
##### 0- Archivos y cosas generales                                                                                   #####
##### 1- Analisis STING abierta                                                                                       #####
##### 2- Analisis STING cerrada con ligando                                                                           #####
##### 3- Analisis STING cerrada sin ligando                                                                           #####
##### 4- Analisis STING cerrada con ligando rotado                                                                    #####
##### 5- Comparaciones                                                                                                #####
##### 6- Distancia de glicinas 163                                                                                    #####
##### 7- RMSD de distintas regiones                                                                                   #####
##### 8- B-factors                                                                                                    #####
##### 9- Analisis de componentes principales (PCA)                                                                    #####
##### 10- Puentes de hidrógeno                                                                                        #####
##### 11- Mas cosas del sitio de polimerizacion                                                                       #####
#####                                                                                                       -hoseki   #####
###########################################################################################################################
###########################################################################################################################

############################################# 0- ARCHIVOS Y COSAS GENERALES   ###############################################
setwd("~/Desktop/Sysbio/cSTING")

#THEME FOR PLOTS
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic'),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

#ARCHIVOS 
rms_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt6)[names(rms_6nt6)=="TIME.ps."] <- "Time (ns)"
names(rms_6nt6)[names(rms_6nt6)=="RMSD.nm."] <- "RMSD (nm)"
rms_6nt6$`Time (ns)` = rms_6nt6$`Time (ns)`/1000
rmsf_6nt6 <- read.csv("6nt6_NO_BORRAR/400ns_NpT-WT_6nt6_nojump_rot_rmsf.xvg", sep="")
i=146
for (j in 1:402){
  rmsf_6nt6$ATOM[j]=i
  i=i+1
}
gr_6nt6 <- read.csv("6nt6_NO_BORRAR/500ns_NpT-WT_6nt6_nojump_rot_GR.xvg", sep="")
gr_6nt6$Time = gr_6nt6$Time/1000
pdb_6nt6 <- read.pdb("6nt6_NO_BORRAR/6nt6_ns100-500.pdb")
dcd_6n6 <- read.dcd("6nt6_NO_BORRAR/6nt6_ns100-500.dcd")  

rms_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="Time"] <- "Time (ns)"
names(rms_6nt7_holo)[names(rms_6nt7_holo)=="RMSDns"] <- "RMSD (nm)"
rmsf_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/400ns_6nt7_holo_NpT_nojump_rot_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_holo$ATOM[j]=i
  i=i+1
}
gr_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_gr.xvg", sep="")
gr_6nt7_holo$Time = gr_6nt7_holo$Time/1000
pdb_6nt7_holo <- read.pdb("6nt7_HOLO_NO_BORRAR/6nt7_holo_ns100-500.pdb")
dcd_6nt7_holo <- read.dcd("6nt7_HOLO_NO_BORRAR/6nt7_holo_ns100-500.dcd")

rms_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/500ns_NpT_6nt7_apo_WO_nojump_rot_rmsd.xvg", sep="")
names(rms_6nt7_apo)[names(rms_6nt7_apo)=="Time"] <- "Time (ns)"
names(rms_6nt7_apo)[names(rms_6nt7_apo)=="RMSDns"] <- "RMSD (nm)"
rmsf_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/400ns_NpT_6nt7_apo_wo_nojump_rot_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_apo$ATOM[j]=i
  i=i+1
}
gr_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/500ns_NpT_6nt7_apo_WO_nojump_rot_GR.xvg", sep="")
gr_6nt7_apo$Time = gr_6nt7_apo$Time/1000
pdb_6nt7_apo <- read.pdb("6nt7_APO_NO_BORRAR/6nt7_apo_ns100:500.pdb")
dcd_6nt7_apo <- read.dcd("6nt7_APO_NO_BORRAR/6nt7_apo_ns100:500.dcd") 

pol_6nt6 <- read.csv("6nt6_NO_BORRAR/6nt6_pol_rmsd.xvg", sep="")
pol_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/6nt7_apo_polsite_rmsd.xvg", sep="")
pol_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/6nt7_holo_pol_rmsd.xvg", sep="")

helix_6nt6 <- read.csv("6nt6_NO_BORRAR/6nt6_con_hel_rmsd.xvg", sep="")
helix_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/6nt7_apo_conHelix_rmsd.xvg", sep="")
helix_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/6nt7_holo_conHelix_rmsd.xvg", sep="")

loop_6nt6 <- read.csv("6nt6_NO_BORRAR/6nt6_con_loop_rmsd.xvg", sep="")
loop_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/6nt7_apo_conLoop_rmsd.xvg", sep="")
loop_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/6nt7_holo_conLoop_rmsd.xvg", sep="")

gly_6nt6 <- read.csv("6nt6_NO_BORRAR/6nt6_gly163_dist.xvg", sep="")
gly_6nt7_apo <- read.csv("6nt7_APO_NO_BORRAR/6nt7_apo_gly163_dist.xvg", sep="")
gly_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/6nt7_holo_gly163_dist.xvg", sep="")

bfac_6nt6 <- read.pdb("6nt6_NO_BORRAR/400ns_NpT-WT_6nt6_nojump_rot_bfac.pdb")
bfac_6nt7_holo <- read.pdb("6nt7_HOLO_NO_BORRAR/400ns_6nt7_holo_NpT_nojump_rot_bfac.pdb")
bfac_6nt7_apo <- read.pdb("6nt7_APO_NO_BORRAR/400ns_NpT_6nt7_apo_wo_nojump_rot_bfac.pdb")

nbond <- read.csv("6nt7_HOLO_NO_BORRAR/holo_num_hbond.xvg", sep="")
hbonds <- read.table("6nt7_HOLO_NO_BORRAR/hbonds-lig-ns100:500.dat", quote="\"", comment.char="")
names(hbonds)[names(hbonds)=="V1"] <- "TIME"
names(hbonds)[names(hbonds)=="V2"] <- "HBOND"
nbond[2502,] = nbond[2501,]
nbond <- nbond[-c(1:500),]
hbonds$TIME <- nbond$TIME
NFOP <- read.delim("6nt7_HOLO_NO_BORRAR/hbonds-details-lig-ns100:500.dat")
hbonds_6nt6 <- read.table("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds_all-system.dat", quote="\"", comment.char="")
names(hbonds_6nt6)[] <- c("time", "hbonds")
hbonds_6nt6$time <- c(0, rms_6nt6$`Time (ns)`)
hbonds_6nt7_holo <- read.table("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds_all-system.dat", quote="\"", comment.char="")
names(hbonds_6nt7_holo)[] <- c("time", "hbonds")
hbonds_6nt7_holo$time <- c(0, rms_6nt6$`Time (ns)`)

lig_rot_rmsd <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-rmsd.xvg", sep="")
lig_rot_rmsf <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-rmsf.xvg", sep="")
i=146
for (j in 1:394){
  lig_rot_rmsf$res[j]=i
  i=i+1
}
lig_rot_gr <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-gr.xvg", sep="")
lig_rot_gr$time <- lig_rot_gr$time/1000
lig_rot_bfac <- read.pdb("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-500ns-bfac.pdb")
lig_rot_gly <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-rot-gly163.xvg", sep="")
lig_rot_gly$TIME = lig_rot_gly$TIME/1000
lig_rot_pdb <- read.pdb("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-400ns.pdb")
lig_rot_dcd <- read.dcd("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada-lig-rot-400ns.dcd")

############################################## 1- ANALISIS STING ABIERTA ##################################################
#RMSD
ggplot(rms_6nt6, aes(x=rms_6nt6$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt6$`RMSD (nm)`), color = "mediumseagreen", size = 0.8) + 
  ggtitle("RMSD OPEN STING (Z RESTRICTION)") + ylim(0,0.6) +
  xlab("Time (ns)") + 
  ylab("RMSD (nm)") + 
  a

mean(rms_6nt6$`RMSD (nm)`[501:2501])
sd(rms_6nt6$`RMSD (nm)`[501:2501])
#RMSF
ggplot(rmsf_6nt6, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "mediumseagreen", size=0.8 ) + ggtitle("RMSF OPEN STING (Z RESTRICTION)") +
  xlab("Residue") + ylab("RMSF (nm)") + a +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                446, 456, 466, 476, 486, 496, 506, 516, 526, 536, 546),
                     labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                      "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                      "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                      "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                      "TYR346A", "SER155B", "ALA165B", "VAL175B", "GLU185B", 
                                      "ARG195B", "ILE205B", "ASP215B", "GLN225B", "LEU235B", 
                                      "TYR245B", "LYS255B", "GLU265B", "ALA275B", "SER285B", 
                                      "LEU295B", "GLY205B", "ARG315B", "PRO325B", "LEU335B", "GLU345B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax= 1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax= 1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)


#RADIO DE GIRO
ggplot(gr_6nt6, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt6$RG), color="mediumseagreen", size=0.8) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (Z RESTRICTION)") + xlab("Time (ns)") +   ylab("RG (nm)") + a + ylim(2.15,2.35) +
  annotate("rect", xmin = 60, xmax=105, ymin=2.15, ymax= 2.35, alpha=.2, fill="mediumseagreen") +
  annotate("rect", xmin = 290, xmax=410, ymin=2.15, ymax= 2.35, alpha=.2, fill="mediumseagreen")

ggplot(gr_6nt6, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt6$RG), color="mediumseagreen", size=0.8) + 
  geom_line(aes(y=gr_6nt6$RGX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt6$RGY), color="magenta2", size=0.8) + 
  geom_line(aes(y=gr_6nt6$RGZ), color="orange", size=0.8) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (Z RESTRICTION)") + xlab("Time (ns)") + 
  ylab("GR (nm)") + a +
  annotate("rect", xmin = 60, xmax=105, ymin=1.4, ymax= 2.4, alpha=.2, fill="mediumseagreen") +
  annotate("rect", xmin = 290, xmax=410, ymin=1.4, ymax= 2.4, alpha=.2, fill="mediumseagreen") +
  annotate("text", x=10, y=1.6, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.88, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.15, label="Y-RG", color="magenta2")
  

#70-105
#280-400
#400-500
############################################# 2- ANALISIS STING CERRADA CON LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_holo, aes(x=rms_6nt7_holo$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_holo$`RMSD (nm)`), color = "mediumorchid2", size=0.8) + 
  ggtitle("RMSD CLOSE STING + LIGAND (Z RESTRICTION)") + ylim(0,0.6)+
  xlab("Time (ns)") + ylab("RMSD (nm)") + a 
mean(rms_6nt7_holo$`RMSD (nm)`[501:2501])
sd(rms_6nt7_holo$`RMSD (nm)`[501:2501])
#RMSF
ggplot(rmsf_6nt7_holo, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "mediumorchid2", size=0.8) + ggtitle("RMSF CLOSE STING + LIGAND (Z RESTRICTION)") + 
  xlab("Residue") + ylab("RMSF (nm)") + a +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                     labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                      "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                      "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                      "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                      "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                      "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                      "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                      "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
  theme(axis.text.x = element_text(angle = 90, face="bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

#RADIO DE GIRO
ggplot(gr_6nt7_holo, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo$GR), color="mediumorchid2", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (Z RESTRICTION)") + ylim(2.15,2.35) +
  xlab("Time (ns)") + ylab("RG (nm)") + a

ggplot(gr_6nt7_holo, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo$GR), color="mediumorchid2", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo$GRY), color="brown", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo$GRZ), color="orange", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (Z RESTRICTION)") + 
  xlab("Time (ns)") + ylab("RG (nm)") + a +
  annotate("text", x=10, y=1.6, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.8, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.1, label="Y-RG", color="brown")
############################################# 3- ANALISIS STING CERRADA SIN LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_apo, aes(x=rms_6nt7_apo$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_apo$`RMSD (nm)`), color = "#110BD0", size=0.8) + 
  ggtitle("RMSD CLOSE STING NO LIGAND (Z RESTRICTION)") + a + ylim(0,0.6) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + 
  annotate("rect", xmin=360, xmax=390, ymin=0, ymax= 0.6, alpha=.2, fill="#110BD0") + 
  annotate("rect", xmin = 25, xmax=55, ymin=0, ymax= 0.6, alpha=.2, fill="#110BD0")

#RMSF
ggplot(rmsf_6nt7_apo, aes(x=ATOM, y=nm)) + 
  geom_line(color = "#110BD0", size=0.8) + ggtitle("RMSD CLOSE STING NO LIGAND (Z RESTRICTION)") + 
  xlab("Residues") + ylab("RMSF (nm)") + a + ylim(0,1)  +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                     labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                      "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                      "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                      "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                      "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                      "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                      "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                      "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 1, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

#RADIO DE GIRO
ggplot(gr_6nt7_apo, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo$GR), color="#110BD0", size=0.8) + ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (Z RESTRICTION)") + a +
  xlab("Time (ns)") + ylab("RG (nm)") + ylim(2.15,2.35) +
  annotate("rect", xmin=360, xmax=390, ymin=2.15, ymax= 2.35, alpha=.2, fill="#110BD0") #+ 
  annotate("rect", xmin = 25, xmax=55, ymin=2.15, ymax= 2.35, alpha=.2, fill="#110BD0") 

ggplot(gr_6nt7_apo, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo$GR), color="#110BD0", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo$GRY), color="magenta3", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo$GRZ), color="orange", size=0.8) + 
  xlab("Time (ns)") + ylab("RG (nm)") + ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (Z RESTRICTION)") +
  annotate("rect", xmin=360, xmax=390, ymin=1.5, ymax= 2.3, alpha=.2, fill="#110BD0") + 
  annotate("text", x=10, y=1.6, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.75, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.05, label="Y-RG", color="magenta3") + a
########################################### 4- ANALISIS STING CERRADA CON LIGANDO ROTADO ##################################
## RMSD
ggplot(lig_rot_rmsd, aes(x=time)) + 
  geom_line(aes(y=rmsd), color = "dodgerblue3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a 
mean(lig_rot_rmsd$rmsd[501:2501])
sd(lig_rot_rmsd$rmsd[501:2501])

  ## RMSF
ggplot(lig_rot_rmsf, aes(x=res, y=rmsf)) +
  geom_line(color = "dodgerblue3", size=0.8) + 
  ggtitle("RMSF CLOSE STING + ROTATED LIGAND (Z RESTRICTION)") + 
  xlab("Residue") + ylab("RMSF (nm)") + a +
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                     labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                      "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                      "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                      "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                      "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                      "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                      "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                      "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
  theme(axis.text.x = element_text(angle = 90, face="bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.5, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

## RADIO DE GIRO
ggplot(lig_rot_gr, aes(x=time)) + 
  geom_line(aes(y=all), color="dodgerblue3", size=0.8) + 
  geom_line(aes(y=x), color="blue", size=0.8) + 
  geom_line(aes(y=y), color="brown", size=0.8) + 
  geom_line(aes(y=z), color="orange", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + ROTATED LIGAND (Z RESTRICTION)") + 
  xlab("Time (ns)") + ylab("RG (nm)") + a +
  annotate("text", x=10, y=1.6, label="X-RG", color="blue") +
  annotate("text", x=10, y=1.8, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.1, label="Y-RG", color="brown")
pos <- seq(from = 1, to = 1001, by = 5)
ggplot(lig_rot_gr, aes(x=time)) + ylim(2,2.5) +
  geom_line(aes(y = all), color="dodgerblue3", size=0.8) +
  xlab("Time (ns)") + ylab("RG (nm)") + a 
################################################## 5- COMPARACIONES #######################################################
#RMSD
RMSD <- data.frame("Time(ns)" = rms_6nt6$`Time (ns)`, 
                   "6NT6" = rms_6nt6$`RMSD (nm)`, 
                   "6NT7 APO" = rms_6nt7_apo$`RMSD (nm)`, 
                   "6NT7 HOLO" = rms_6nt7_holo$`RMSD (nm)`)
ggplot(RMSD, aes(x=RMSD$Time.ns.)) + 
  geom_line(aes(y=RMSD$X6NT6), color = "mediumseagreen", size=0.8) + 
  geom_line(aes(y=RMSD$X6NT7.HOLO), color = "mediumorchid2", size=0.8) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + 
  annotate("text", x=500, y=0.6, label="APO-STING", color="mediumseagreen", hjust = 1) +
  annotate("text", x=500, y=0.55, label="HOLO-STING", color="mediumorchid2", hjust = 1)
abierta_rmsd <- data.frame("time" = rms_6nt6$`Time (ns)`, 
                           "rmsd" = rms_6nt6$`RMSD (nm)`, 
                           "label" = as.character("APO-STING"))
cerrada_rms <- data.frame("time" = rms_6nt7_holo$`Time (ns)`,
                          "rmsd" = rms_6nt7_holo$`RMSD (nm)`,
                          "label" = "HOLO-STING")
rotada_rms <- data.frame("time" = lig_rot_rmsd$time,
                         "rms" = lig_rot_rmsd$rmsd, 
                         "label" = "HOLO-STING rot")
rms <- data.frame("time" = c(abierta_rmsd$time, cerrada_rms$time, rotada_rms$time),
                  "rmsd" = c(abierta_rmsd$rmsd, cerrada_rms$rmsd, rotada_rms$rms),
                  "Structure" = c(as.character(abierta_rmsd$label), as.character(cerrada_rms$label), as.character(rotada_rms$label)))
A <- ggplot(rms[0:5002,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  xlab("Time (ns)") + ylab("RMSD (nm)")
  

#RMSF
d <- c(343, 344, 345, 346, 544, 545, 546, 547)
b <- c(1:8)

for (l in 1:8){
  b[l] <- which(rmsf_6nt6$ATOM == d[l])
}

apo <- rmsf_6nt6[c(-b),]
i=146

for (j in 1:394){
  apo$ATOM[j]=i
  i=i+1
}

RMSF <- data.frame("aa" = apo$ATOM, 
                   "6nt6" = apo$nm, 
                   "6nt7-apo" = rmsf_6nt7_apo$nm, 
                   "6nt7-holo"= rmsf_6nt7_holo$nm)
ggplot(RMSF, aes(x=RMSF$aa)) +
  geom_line(aes(y=RMSF$X6nt6), color = 'mediumseagreen', size=0.8) +
  geom_line(aes(y=RMSF$X6nt7.apo), color ='#110BD0', size=0.8) +
  geom_line(aes(y=RMSF$X6nt7.holo), color ='mediumorchid2', size=0.8) +
  ylab('RMSF (nm)') + xlab('aminoacid') + a +
  ggtitle('Comparación RSMF')

#RADIO DE GIRO
GR <- data.frame("Time (ns)"= gr_6nt6$Time, 
                 "6NT6"=gr_6nt6$RG, 
                 "6NT7-APO"= gr_6nt7_apo$GR, 
                 "6NT7-HOLO"= gr_6nt7_holo$GR,
                 "label abierta" = "APO-STING",
                 "label cerrada" = "HOLO-STING")
ggplot(GR, aes(x= GR$Time..ns.)) + 
  geom_line(aes(y=GR$X6NT6), color="mediumseagreen", size=0.8) + 
  geom_line(aes(y=GR$X6NT7.APO), color="#110BD0", size=0.8) + 
  geom_line(aes(y=GR$X6NT7.HOLO), color="mediumorchid2", size=0.8) +
  ggtitle("Comparación Radio de Giro") +
  xlab("Time (ns)") + ylab("GR (nm)") + a
gr <- data.frame("time" = c(GR$Time..ns., GR$Time..ns.),
                 "gr" = c(GR$X6NT6, GR$X6NT7.HOLO), 
                 "Structure" = c(as.character(GR$label.abierta), as.character(GR$label.cerrada)))
B <- ggplot(gr, aes(x = time, y = gr)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2")) +
  xlab("Time (ns)") + ylab("Radius of gyration (nm)") + ylim(2,2.5) 

ggarrange(A, B, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

############################################### 6- DISTANCIA DE GLICINAS 163 ##############################################
pos <- seq(from = 1, to = 2501, by = 5)
ggplot(gly_6nt6[pos,], aes(x=TIME)) + 
  geom_line(aes(y=DIS), color = "mediumseagreen", size=0.8, alpha = .7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a +
  ggtitle("GLYCINE 163 DISTANCE OPEN STING (Z RESTRICTION)") #+
  geom_line(data = GLY_AVG, aes(x = Time/1000, y = avg_6NT6), color="darkgreen")

ggplot(gly_6nt7_apo, aes(x=gly_6nt7_apo$TIME/1000)) + 
  geom_line(aes(y=gly_6nt7_apo$DIS), color = "#110BD0", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING NO LIGAND (Z RESTRICTION)")# +
  geom_line(data = GLY_AVG, aes(x = Time/1000, y = avg_APO), color="blue")

ggplot(gly_6nt7_holo[pos,], aes(x=TIME)) + 
  geom_line(aes(y=DIS), color = "mediumorchid2", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING + LIGAND (Z RESTRICTION)")# +
  geom_line(data = GLY_AVG, aes(x = Time/1000, y = avg_HOLO), color="slateblue3")
  

ggplot(GLY[pos,], aes(x = TIME)) + 
  geom_line(aes(y=abierta), color = "mediumseagreen", size=0.8, alpha = .7) + 
  geom_line(aes(y=cerrada), color = "mediumorchid2", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING + LIGAND (Z RESTRICTION)") +
  annotate("text", x=450, y=1, label="Open-STING", color="mediumseagreen") +
  annotate("text", x=450, y=1.05, label="Close-STING", color="mediumorchid2")

ggplot(lig_rot_gly[pos,], aes(x=TIME)) + 
  geom_line(aes(y=DIS), color = "magenta4", size=0.8, alpha = 1) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING + ROTATED LIGAND (Z RESTRICTION)") 
  

GLY <- data.frame("TIME" = gly_6nt6$TIME, 
                  "abierta" = gly_6nt6$DIS,
                  "cerrada" = gly_6nt7_holo$DIS)


GLY_AVG <- data.frame("Time" = c(0:50), "GLY_6NT6" = 0, "GLY_HOLO"=0, "GLY_APO"=0, "count_6NT6" = 0, "count_APO" = 0, "count_HOLO" = 0)
GLY_AVG$Time = GLY_AVG$Time*10000
GLY_AVG$GLY_6NT6 <- as.integer(GLY_AVG$GLY_6NT6)
GLY_AVG$GLY_HOLO <- as.integer(GLY_AVG$GLY_HOLO)
GLY_AVG$GLY_APO <- as.integer(GLY_AVG$GLY_APO)
j = 1
for(i in 1: 2501){
  if(gly_6nt6$TIME[i] == GLY_AVG$Time[j]){
    GLY_AVG$GLY_6NT6[j] = GLY_AVG$GLY_6NT6[j] + gly_6nt6$DIS[i]
    GLY_AVG$count_6NT6[j] = GLY_AVG$count_6NT6[j] + 1
  }
  else{
    if(gly_6nt6$TIME[i] <GLY_AVG$Time[j]){
      GLY_AVG$GLY_6NT6[j] = GLY_AVG$GLY_6NT6[j] + gly_6nt6$DIS[i]
      GLY_AVG$count_6NT6[j] = GLY_AVG$count_6NT6[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG$GLY_6NT6[j] = GLY_AVG$GLY_6NT6[j] + gly_6nt6$DIS[i]
      GLY_AVG$count_6NT6[j] = GLY_AVG$count_6NT6[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_holo$TIME[i] == GLY_AVG$Time[j]){
    GLY_AVG$GLY_HOLO[j] = GLY_AVG$GLY_HOLO[j] + gly_6nt7_holo$DIS[i]
    GLY_AVG$count_HOLO[j] = GLY_AVG$count_HOLO[j] + 1
  }
  else{
    if(gly_6nt7_holo$TIME[i] <GLY_AVG$Time[j]){
      GLY_AVG$GLY_HOLO[j] = GLY_AVG$GLY_HOLO[j] + gly_6nt7_holo$DIS[i]
      GLY_AVG$count_HOLO[j] = GLY_AVG$count_HOLO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG$GLY_HOLO[j] = GLY_AVG$GLY_HOLO[j] + gly_6nt7_holo$DIS[i]
      GLY_AVG$count_HOLO[j] = GLY_AVG$count_HOLO[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_apo$TIME[i] == GLY_AVG$Time[j]){
    GLY_AVG$GLY_APO[j] = GLY_AVG$GLY_APO[j] + gly_6nt7_apo$DIS[i]
    GLY_AVG$count_APO[j] = GLY_AVG$count_APO[j] + 1
  }
  else{
    if(gly_6nt7_apo$TIME[i] <GLY_AVG$Time[j]){
      GLY_AVG$GLY_APO[j] = GLY_AVG$GLY_APO[j] + gly_6nt7_apo$DIS[i]
      GLY_AVG$count_APO[j] = GLY_AVG$count_APO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG$GLY_APO[j] = GLY_AVG$GLY_APO[j] + gly_6nt7_apo$DIS[i]
      GLY_AVG$count_APO[j] = GLY_AVG$count_APO[j] + 1
    }
  }
}



GLY_AVG$avg_6NT6 <- GLY_AVG$GLY_6NT6/GLY_AVG$count_6NT6
GLY_AVG$avg_HOLO <- GLY_AVG$GLY_HOLO/GLY_AVG$count_HOLO
GLY_AVG$avg_APO <- GLY_AVG$GLY_APO/GLY_AVG$count_APO



ggplot(GLY, aes(x=GLY$TIME)) +
  geom_line(aes(y= GLY$X6NT6), color = "mediumseagreen", size=0.8) + 
  geom_line(aes(y= GLY$X6NT7.HOLO), color = "mediumorchid2", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") + a 

ggplot(GLY, aes(x=GLY$TIME)) +
  geom_line(aes(y= GLY$X6NT6), color = "mediumseagreen", size=0.8, alpha = 0.6) + 
  geom_line(aes(y= GLY$X6NT7.HOLO), color = "mediumorchid2", size=0.8, alpha=0.6) +
  geom_line(aes(y= GLY$X6NT7.APO), color = "#110BD0", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") + a 

############################################### 7- RMSD DE DISTINTOS DOMINIOS #############################################
#SITIO DE POLIMERIZACION
ggplot(pol_6nt6, aes(x=pol_6nt6$TIME)) + 
  geom_line(aes(y=pol_6nt6$RMSD), color = "mediumseagreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.7) + 
  ggtitle("RMSD POLYMERIZATION DOMAIN OPEN STING (Z RESTRICTION)")

ggplot(pol_6nt7_apo, aes(x=pol_6nt7_apo$TIME)) + 
  geom_line(aes(y=pol_6nt7_apo$RMSD), color = "#110BD0", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.7)  + 
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING NO LIGAND (Z RESTRICTION)")

ggplot(pol_6nt7_holo, aes(x=pol_6nt7_holo$TIME)) + 
  geom_line(aes(y=pol_6nt7_holo$RMSD), color = "mediumorchid2", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.7) + 
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING + LIGAND (Z RESTRICTION)")

#HELICE CONTECTORA
ggplot(helix_6nt6, aes(x=helix_6nt6$TIME)) + 
  geom_line(aes(y=helix_6nt6$RMSD), color = "mediumseagreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.65) + 
  ggtitle("RMSD CONECTOR HELIX OPEN STING (Z RESTRICTION)")

ggplot(helix_6nt7_apo, aes(x=helix_6nt7_apo$TIME)) + 
  geom_line(aes(y=helix_6nt7_apo$RMSD), color = "#110BD0", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.65)  + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING NO LIGAND (Z RESTRICTION)")

ggplot(helix_6nt7_holo, aes(x=helix_6nt7_holo$TIME)) + 
  geom_line(aes(y=helix_6nt7_holo$RMSD), color = "mediumorchid2", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.65)  + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING + LIGAND (Z RESTRICTION)")

#LOOP CONECTOR
ggplot(loop_6nt6, aes(x=loop_6nt6$TIME)) + 
  geom_line(aes(y=loop_6nt6$RMSD), color = "mediumseagreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.30) + 
  ggtitle("RMSD CONECTOR LOOP OPEN STING (Z RESTRICTION)")

ggplot(loop_6nt7_apo, aes(x=loop_6nt7_apo$TIME)) + 
  geom_line(aes(y=loop_6nt7_apo$RMSD), color = "#110BD0", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.30) +  
  ggtitle("RMSD CONECTOR LOOP CLOSE STING NO LIGAND (Z RESTRICTION)")

ggplot(loop_6nt7_holo, aes(x=loop_6nt7_holo$TIME)) + 
  geom_line(aes(y=loop_6nt7_holo$RMSD), color = "mediumorchid2", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.30) + 
  ggtitle("RMSD CONECTOR LOOP CLOSE STING + LIGAND (Z RESTRICTION)")

#################################################### 8- B-FACTORS #########################################################
  #6nt6
  bfac_6nt6_db <- data.frame("NUM" = bfac_6nt6$atom$resno, 
                             "RES" = bfac_6nt6$atom$resid, 
                             "BFAC" = bfac_6nt6$atom$b)
  prom_6nt6 <- mean(bfac_6nt6_db$BFAC)
  sd_6nt6 <- sd(bfac_6nt6_db$BFAC)
  bfac_6nt6_db$NBFAC <- (bfac_6nt6_db$BFAC - prom_6nt6)/sd_6nt6
  bfac_6nt6_db$RMSF <- rmsf_6nt6$nm
  bfac_6nt6_db <- bfac_6nt6_db[-c(198:201, 399:402),]
  
  j = 146
  for (i in 1:394){
    bfac_6nt6_db$ATOM[i] = j
    j = j+1
  }
  
  ggplot(bfac_6nt6_db, aes(x = bfac_6nt6_db$ATOM)) +
    geom_line(aes(y = bfac_6nt6_db$NBFAC), color = "mediumseagreen", size = 0.8) +
    a + xlab("Residue") + ylab("Normalized B-factor") + ggtitle("NORMALIZED B-FACTORS OPEN STING (Z RESTRICTION)") +
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5) +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536, 546),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "TYR346A", "SER155B", "ALA165B", "VAL175B", "GLU185B", 
                                        "ARG195B", "ILE205B", "ASP215B", "GLN225B", "LEU235B", 
                                        "TYR245B", "LYS255B", "GLU265B", "ALA275B", "SER285B", 
                                        "LEU295B", "GLY205B", "ARG315B", "PRO325B", "LEU335B", "GLU345B"))) +
    theme(axis.text.x = element_text(angle = 90, size=rel(0.7), face = "bold")) +
    theme(panel.grid.major.x = element_blank(),          
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray95")) +
    annotate("rect", xmin = 189, xmax=200, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 390, xmax=401, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 232, xmax=244, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 433, xmax=445, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 254, xmax=258, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 455, xmax=459, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 322, xmax=328, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 523, xmax=529, ymin=-0.4, ymax= 9.2, alpha=.2, fill="grey60") +
    geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)
  
  bfac_6nt6_75 <- which(bfac_6nt6_db$NBFAC > 0.75)
  formattable(bfac_6nt6_db[bfac_6nt6_75,])
  bfac_6nt6_50 <- which(bfac_6nt6_db$NBFAC > 0.5)
  bfac_6nt6_db[bfac_6nt6_50,]
  
  #6nt7 holo
  bfac_6nt7_holo_db <- data.frame("NUM" = bfac_6nt7_holo$atom$resno, 
                                  "RES" = bfac_6nt7_holo$atom$resid, 
                                  "BFAC" = bfac_6nt7_holo$atom$b)
  prom_6nt7_holo <- mean(bfac_6nt7_holo_db$BFAC)
  sd_6nt7_holo <- sd(bfac_6nt7_holo_db$BFAC)
  bfac_6nt7_holo_db$NBFAC <- (bfac_6nt7_holo_db$BFAC - prom_6nt7_holo)/sd_6nt7_holo
  bfac_6nt7_holo_db$RMSF <- rmsf_6nt7_holo$nm
  
  j = 146
  for (i in 1:394){
    bfac_6nt7_holo_db$ATOM[i] = j
    j = j+1
  }
  
  ggplot(bfac_6nt7_holo_db, aes(x = bfac_6nt7_holo_db$ATOM)) +
    geom_line(aes(y = bfac_6nt7_holo_db$NBFAC), color = "mediumorchid2", size = 0.8) +
    a + xlab("Residue") + ylab("Normalized B-factor") + ggtitle("NORMALIZED B-FACTORS CLOSE STING + LIGAND (Z RESTRICTION)") + 
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5)  +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                        "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                        "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                        "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
    theme(axis.text.x = element_text(angle = 90, size=rel(0.7), face = "bold")) +
    annotate("rect", xmin = 191, xmax=196, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 388, xmax=393, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 255, xmax=258, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 452, xmax=455, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 323, xmax=327, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 520, xmax=524, ymin=-0.48, ymax= 10.15, alpha=.2, fill="grey60") +
    theme(panel.grid.major.x = element_blank(),          
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray95")) + 
    geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)
  
  bfac_6nt7_holo_75 <- which(bfac_6nt7_holo_db$NBFAC > 0.75)
  formattable(bfac_6nt7_holo_db[bfac_6nt7_holo_75,])
  bfac_6nt7_holo_50 <- which(bfac_6nt7_holo_db$NBFAC > 0.5)
  bfac_6nt7_holo_db[bfac_6nt7_holo_50,]
  
  #6nt7 apo
  bfac_6nt7_apo_db <- data.frame("NUM" = bfac_6nt7_apo$atom$resno, 
                                 "RES" = bfac_6nt7_apo$atom$resid, 
                                 "BFAC" = bfac_6nt7_apo$atom$o) 
  bfac_6nt7_apo_db$BFAC <- bfac_6nt7_apo_db$BFAC + bfac_6nt7_apo$atom$b
  prom_6nt7_apo <- mean(bfac_6nt7_apo_db$BFAC)
  sd_6nt7_apo <- sd(bfac_6nt7_apo_db$BFAC)
  bfac_6nt7_apo_db$NBFAC <- (bfac_6nt7_apo_db$BFAC - prom_6nt7_apo)/sd_6nt7_apo
  bfac_6nt7_apo_db$RMSF <- rmsf_6nt7_apo$nm
  
  j = 146
  for (i in 1:394){
    bfac_6nt7_apo_db$ATOM[i] = j
    j = j+1
  }
  
  ggplot(bfac_6nt7_apo_db, aes(x = bfac_6nt7_apo_db$ATOM)) +
    geom_line(aes(y = bfac_6nt7_apo_db$NBFAC), color = "#110BD0", size = 0.8) +
    a + xlab("Residue") + ylab("Normalized B-Factor") + ggtitle("NORMALIZED B-FACTORS CLOSE STING NO LIGAND (Z RESTRICTION)") +
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5)  +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                        "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                        "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                        "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
    theme(axis.text.x = element_text(angle = 90, size=rel(0.7), face = "bold")) +
    annotate("rect", xmin = 191, xmax=196, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 388, xmax=393, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 255, xmax=258, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 452, xmax=455, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 323, xmax=327, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 520, xmax=524, ymin=-0.45, ymax= 9, alpha=.2, fill="grey60") +
    theme(panel.grid.major.x = element_blank(),          
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray95")) + 
    geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)
  
  bfac_6nt7_apo_75 <- which(bfac_6nt7_apo_db$NBFAC > 0.75)
  bfac_6nt7_apo_db[bfac_6nt7_apo_75,]
  bfac_6nt7_apo_50 <- which(bfac_6nt7_apo_db$NBFAC > 0.5)
  bfac_6nt7_apo_db[bfac_6nt7_apo_50,]
  
  ggplot() +
    geom_line(data = bfac_6nt7_apo_db, aes(x = bfac_6nt7_apo_db$ATOM, y = bfac_6nt7_apo_db$NBFAC), color = "#110BD0", size = 0.8) +
    geom_line(data = bfac_6nt7_holo_db, aes(x= bfac_6nt7_holo_db$ATOM, y= bfac_6nt7_holo_db$NBFAC), color = "mediumorchid2", size = 0.8) +
    geom_line(data = bfac_6nt6_db, aes(x = bfac_6nt6_db$ATOM, y = bfac_6nt6_db$NBFAC), color = "mediumseagreen", size = 0.8 ) +
    a + xlab("Residuo") + ylab("B-factor normalizado")  +
    geom_hline(yintercept=0.75, linetype="dashed", color = "black") +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                        "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                        "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                        "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
    theme(axis.text.x = element_text(angle = 90, face = "bold")) 
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
  
  
  ### colors -> "mediumseagreen", "green3","mediumorchid2", "purple4", "dodgerblue3", "blue3"
  ### scale_linetype_manual(values=c("solid", "twodash", "solid", "twodash", "solid", "twodash")) +
ggplot(bfac[c(395:591, 789:985),], aes(x = residue, y = nbfac)) + 
    geom_line(aes(color = Structure, linetype = Structure), size = 0.8) + 
    scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
    xlab("Residue") + ylab("Normalized B-factor") +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                        "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                        "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                        "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
    theme(axis.text.x = element_text(angle = 90, face = "bold")) +
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 1) +
    annotate("rect", xmin = 191, xmax=196, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") + 
    annotate("rect", xmin = 388, xmax=393, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") +
    annotate("rect", xmin = 255, xmax=258, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") +
    annotate("rect", xmin = 452, xmax=455, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") +
    annotate("rect", xmin = 323, xmax=327, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") +
    annotate("rect", xmin = 520, xmax=524, ymin=-0.45, ymax= 10.15, alpha=.3, fill="hotpink4") +
    annotate("rect", xmin = 189, xmax=200, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") + 
    annotate("rect", xmin = 386, xmax=397, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    annotate("rect", xmin = 232, xmax=244, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    annotate("rect", xmin = 429, xmax=441, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    annotate("rect", xmin = 254, xmax=258, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    annotate("rect", xmin = 451, xmax=455, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") + 
    annotate("rect", xmin = 322, xmax=328, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    annotate("rect", xmin = 519, xmax=525, ymin=-0.45, ymax= 10.15, alpha=.3, fill="greenyellow") +
    theme(
      panel.grid.major = element_line(color = NA), 
      panel.grid.minor = element_line(color = NA),
      plot.title = element_text(colour = 'black', size = 15),
      axis.title = element_text(colour = 'black', face = 'italic'),
      panel.background = element_rect(fill = "white", colour = "black")
    ) +
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5) +
    geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)
    
  
  
    #STING CON LIGANDO ROTADO
  bfac_lig_rot_db <- data.frame("NUM" = lig_rot_bfac$atom$resno, 
                                "RES" = lig_rot_bfac$atom$resid, 
                                "BFAC" = lig_rot_bfac$atom$b)
  prom2 <- mean(bfac_lig_rot_db$BFAC)
  sd2 <- sd(bfac_lig_rot_db$BFAC)
  bfac_lig_rot_db$NBFAC <- (bfac_lig_rot_db$BFAC - prom2)/sd2
  bfac_lig_rot_db$RMSF <- lig_rot_rmsf$rmsf
  
  j = 146
  for (i in 1:394){
    bfac_lig_rot_db$ATOM[i] = j
    j = j+1
  }
  
  ggplot(bfac_lig_rot_db, aes(x = ATOM)) +
    geom_line(aes(y = NBFAC), color = "dodgerblue3", size = 0.8) +
    a + xlab("Residue") + ylab("Normalized B-factor") + ggtitle("NORMALIZED B-FACTORS CLOSE STING + ROTATED LIGAND (Z RESTRICTION)") + 
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5)  +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336, 
                                  346, 356, 366, 376, 386, 396, 406, 416, 426, 436, 
                                  446, 456, 466, 476, 486, 496, 506, 516, 526, 536),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A", 
                                        "VAL149B", "ASN159B", "TYR169B", "ARG179B", "ARG189B", 
                                        "ASP199B", "LEU209B", "LYS219B", "ASP229B", "GLY239B", 
                                        "LEU249B", "LEU259B", "PRO269B", "ASP279B", "ARG289B", 
                                        "SER299B", "GLU309B", "TYR319B", "PHE329B", "GLN339B"))) +
    theme(axis.text.x = element_text(angle = 90, size=rel(0.7), face = "bold")) +
    annotate("rect", xmin = 191, xmax=196, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 388, xmax=393, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 255, xmax=258, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 452, xmax=455, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 323, xmax=327, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 520, xmax=524, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    theme(panel.grid.major.x = element_blank(),          
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray95")) + 
    geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5) + 
    geom_segment(aes(x = 310, y = 8, xend = 286, yend = 4.5), arrow = arrow(length = unit(0.7, "cm")))
  
  bfac_6nt7_holo_75 <- which(bfac_lig_rot_db$NBFAC > 0.75)
  bfac_lig_rot_db[bfac_6nt7_holo_75,]
  bfac_6nt7_holo_50 <- which(bfac_lig_rot_db$NBFAC > 0.5)
  bfac_lig_rot_db[bfac_6nt7_holo_50,]
  
  db <- data.frame("atom" = lig_rot_bfac$atom$resno,
                   "bfac" = lig_rot_bfac$atom$b)
  prom <- mean(db$bfac)
  sd <- sd(db$bfac)
  db$nbfac <- (db$bfac - prom)/sd
  atom <- c(146:342)
  chaina <- db$nbfac[1:197]
  chainb <- db$nbfac[198:394]
  
  db2 <- data.frame("atom" = atom,
                    "chaina" = chaina, 
                    "chainb" = chainb)
  ggplot(db2, aes(x = atom)) + 
    geom_line(aes(y = chaina), color = "blue", size = 0.8) + 
    geom_line(aes(y = chainb), color = "darkgreen", size = 0.8) + 
    xlab("Residue") + ylab("Normalized B-factor") + a + 
    annotate("rect", xmin = 191, xmax=196, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") + 
    annotate("rect", xmin = 255, xmax=258, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") +
    annotate("rect", xmin = 323, xmax=327, ymin=-0.72, ymax= 7.7, alpha=.2, fill="grey60") + 
    theme(axis.text.x = element_text(angle = 90, size=rel(0.7), face = "bold")) +
    geom_line(aes(y = 0.75), color = "red", alpha = 0.5, linetype = 5)  +
    scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                  246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                       labels = paste(c("ALA146A", "SER156A", "TRP166A", "VAL176A", "GLU186A", 
                                        "ALA196A", "LEU206A", "ASP216A", "TYR226A", "THR236A", 
                                        "LYS246A", "ASP256A", "PHE266A", "MET276A", "ARG286A", 
                                        "PHE296A", "SER306A", "LEU316A", "GLU326A", "TRP336A"))) +
    theme(panel.grid.major.x = element_blank(),          
          panel.grid.minor.x = element_blank(), 
          panel.grid.major.y = element_line(color = "gray95")) +
    ggtitle("NORMALIZED B-FACTORS CLOSE STING + ROTATED LIGAND (Z RESTRICTION)") + 
    annotate("text", x=320, y=7, label="Chain A", color="blue") +
    annotate("text", x=320, y=6, label="Chain B", color="darkgreen") +
    geom_segment(aes(x = 300, y = 8, xend = 286, yend = 4.5), arrow = arrow(length = unit(0.7, "cm")))
  
################################## 9- ANALISIS DE COMPONENTES PRINCIPALES (PCA) ###########################################
#6NT6
ca_ndx_6nt6 <- atom.select(pdb_6nt6, elety = "CA") 
pca_6nt6 <- pca.xyz(dcd_6n6[,ca_ndx_6nt6$xyz])
plot(pca_6nt6, col=topo.colors(nrow(dcd_6n6)))
  
pca_6nt6_db <- data.frame("PC" = pca_6nt6$z)
pca_6nt6_db$time <- c(1:2002)
pca_6nt6_db$time[(1:2002)] = hbonds$TIME

ggplot(pca_6nt6_db, aes(x=pca_6nt6_db$PC.1, y=pca_6nt6_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.88%)") + xlab("CP1 (50.33%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
ggplot(pca_6nt6_db, aes(x=pca_6nt6_db$PC.3, y=pca_6nt6_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.88%)") + xlab("CP3 (4.31%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt6_db, aes(x=pca_6nt6_db$PC.1, y=pca_6nt6_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (50.33%)") + ylab("CP3 (4.31%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

eigen_bd_6nt6 <- data.frame("eigen" = pca_6nt6$L)
sum_eig_6nt6 = sum(eigen_bd_6nt6$eigen)
eigen_bd_6nt6$percent_6nt6 <- (eigen_bd_6nt6$eigen/sum_eig_6nt6)*100
eigen_bd_6nt6$num <- c(1:1206)
eigen_bd_6nt6$sum_percent = 0
eigen_bd_6nt6$sum_percent[1] = eigen_bd_6nt6$percent_6nt6[1]

for (i in 2:1206){
  eigen_bd_6nt6$sum_percent[i] = eigen_bd_6nt6$percent_6nt6[i] + eigen_bd_6nt6$sum_percent[i - 1]
}

ggplot(eigen_bd_6nt6[c(1:7),], aes(x = eigen_bd_6nt6$num[c(1:7)], y=eigen_bd_6nt6$percent_6nt6[c(1:7)])) + 
  geom_line(color="mediumseagreen", size=0.8) + 
  geom_point(color="mediumseagreen") + 
  geom_text(aes(y= eigen_bd_6nt6$percent_6nt6[c(1:7)] + 1, 
                label = round(eigen_bd_6nt6$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))
  
pc1_6nt6 <- readPNG("pics/abierta_pc1_restime.png")
d <- rasterGrob(pc1_6nt6, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
pc2_6nt6 <- readPNG("pics/abierta_pc2_restime.png")
f <- rasterGrob(pc2_6nt6, interpolate=TRUE)
  
grid.arrange(arrangeGrob(B, ncol = 1, nrow = 1),
            arrangeGrob(A, ncol = 1, nrow = 1), widths = c(2,3))
  
#mktrj.pca(pca_6nt6, pc=1, b=pca_6nt6$au[,1], file="6nt6_NO_BORRAR/6nt6_pc1_restime.pdb")
#mktrj.pca(pca_6nt6, pc=2, b=pca_6nt6$au[,2], file="6nt6_NO_BORRAR/6nt6_pc2_restime.pdb")
#mktrj.pca(pca_6nt6, pc=3, b=pca_6nt6$au[,3], file="6nt6_NO_BORRAR/6nt6_pc3_restime.pdb")
  
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
  
#threshold (max-min)/2
threshold_6nt6_pc1 = (max(contribution_6nt6$pc1) + min(contribution_6nt6$pc1))/2
contribution_6nt6[c(which(contribution_6nt6$pc1>threshold_6nt6_pc1)),]
threshold_6nt6_pc2 = (max(contribution_6nt6$pc2) + min(contribution_6nt6$pc2))/2
contribution_6nt6[c(which(contribution_6nt6$pc2>threshold_6nt6_pc2)),]
threshold_6nt6_pc3 = (max(contribution_6nt6$pc3) + min(contribution_6nt6$pc3))/2
contribution_6nt6[c(which(contribution_6nt6$pc3>threshold_6nt6_pc3)),]
  
ggplot(contribution_6nt6, aes(x = contribution_6nt6$num)) + 
  geom_col(aes(y = contribution_6nt6$pc1), color = "chartreuse3", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526, 546),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A", 
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "TYR346A", "ALA165B", "GLU185B", "ILE205B", "GLN225B", 
                                      "TYR245B", "GLU265B", "SER285B", "GLY305B", "PRO325B","GLU345B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

  
ggplot(contribution_6nt6, aes(x = contribution_6nt6$num)) + 
  geom_col(aes(y = contribution_6nt6$pc2), color = "darkolivegreen3", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526, 546),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A", 
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "TYR346A", "ALA165B", "GLU185B", "ILE205B", "GLN225B", 
                                      "TYR245B", "GLU265B", "SER285B", "GLY305B", "PRO325B","GLU345B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)
  
ggplot(contribution_6nt6, aes(x = contribution_6nt6$num)) + 
  geom_col(aes(y = contribution_6nt6$pc3), color = "springgreen1", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526, 546),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A", 
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "TYR346A", "ALA165B", "GLU185B", "ILE205B", "GLN225B", 
                                      "TYR245B", "GLU265B", "SER285B", "GLY305B", "PRO325B","GLU345B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) +
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.4, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

#6NT7 HOLO
ca_ndx_6nt7_holo <- atom.select(pdb_6nt7_holo, elety = "CA")
pca_6nt7_holo <- pca.xyz(dcd_6nt7_holo[,ca_ndx_6nt7_holo$xyz])
plot(pca_6nt7_holo, col=topo.colors(nrow(dcd_6nt7_holo)))

pca_6nt7_holo_db <- data.frame("PC" = pca_6nt7_holo$z)
pca_6nt7_holo_db$time <- c(1:2002)
pca_6nt7_holo_db$time[(1:2002)] = hbonds$TIME

ggplot(pca_6nt7_holo_db, aes(x=pca_6nt7_holo_db$PC.1, y=pca_6nt7_holo_db$PC.2, colour=time)) + 
    geom_point(size = 2) + 
    ylab("CP2 (9.81%)") + xlab("CP1 (27.9%)") + a +
    scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt7_holo_db, aes(x=pca_6nt7_holo_db$PC.1, y=pca_6nt7_holo_db$PC.3, colour=time)) + 
    geom_point(size = 2) + 
      ylab("CP3 (5%)") + xlab("CP1 (27.9%)") + a +
    scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt7_holo_db, aes(x=pca_6nt7_holo_db$PC.3, y=pca_6nt7_holo_db$PC.2, colour=time)) + 
    geom_point(size = 2) + 
    ylab("CP2 (9.81%)") + xlab("CP3 (5%)") + a +
    scale_color_gradientn(colours = c("orange", "red", "magenta"))
  
pc1_6nt7_holo <- readPNG("pics/cerrada_con_pc1_restime.png")
j <- rasterGrob(pc1_6nt7_holo, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
  
pc2_6nt7_holo <- readPNG("pics/6nt7_holo_pc2_restime.png")
h <- rasterGrob(pc2_6nt7_holo, interpolate=TRUE)
B <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(h, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
  
grid.arrange(arrangeGrob(A, ncol = 1, nrow = 1),
             arrangeGrob(B, ncol = 1, nrow = 1), widths = c(2,3))
  
eigen_bd_6nt7_holo <- data.frame("eigen" = pca_6nt7_holo$L)
sum_eig_6nt7_holo = sum(eigen_bd_6nt7_holo$eigen)
eigen_bd_6nt7_holo$percent <- (eigen_bd_6nt7_holo$eigen/sum_eig_6nt7_holo)*100
eigen_bd_6nt7_holo$num <- c(1:1182)
eigen_bd_6nt7_holo$sum_percent = 0
eigen_bd_6nt7_holo$sum_percent[1] = eigen_bd_6nt7_holo$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_holo$sum_percent[i] = eigen_bd_6nt7_holo$percent[i] + eigen_bd_6nt7_holo$sum_percent[i - 1]
}

ggplot(eigen_bd_6nt7_holo[c(1:7),], aes(x = eigen_bd_6nt7_holo$num[c(1:7)], 
                                              y=eigen_bd_6nt7_holo$percent[c(1:7)])) + 
  geom_line(color="mediumorchid2", size=0.8) + 
  geom_point(color="mediumorchid2") + 
  geom_text(aes(y= eigen_bd_6nt7_holo$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_holo$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente") + 
  ylab("Acumulated percentaje of variance") + 
  a + scale_x_continuous(breaks = c(1:8)) + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

#mktrj.pca(pca_6nt7_holo, pc=1, b=pca_6nt7_holo$au[,1], file="6nt7_HOLO_NO_BORRAR/6nt7_holo_pc1_restime.pdb")
#mktrj.pca(pca_6nt7_holo, pc=2, b=pca_6nt7_holo$au[,2], file="6nt7_HOLO_NO_BORRAR/6nt7_holo_pc2_restime.pdb")
#mktrj.pca(pca_6nt7_holo, pc=3, b=pca_6nt7_holo$au[,3], file="6nt7_HOLO_NO_BORRAR/6nt7_holo_pc3_restime.pdb")
#mktrj.pca(pca_6nt7_holo, pc=4, b=pca_6nt7_holo$au[,4], file="6nt7_HOLO_NO_BORRAR/6nt7_holo_pc4_restime.pdb")
  
contribution_6nt7_holo <- data.frame("num" = pdb_6nt7_holo$atom$resno[ca_ndx_6nt7_holo$atom],
                                       "resnum" = pdb_6nt7_holo$atom$resno[ca_ndx_6nt7_holo$atom],
                                       "res" = pdb_6nt7_holo$atom$resid[ca_ndx_6nt7_holo$atom], 
                                       "pc1" = pca_6nt7_holo$au[,1], 
                                       "pc2" = pca_6nt7_holo$au[,2], 
                                       "pc3" = pca_6nt7_holo$au[,3]) 
  
i = 146
for (j in 1:392){
  contribution_6nt7_holo$num[j] = i
  i = i+1
}
  
threshold_6nt7_holo_pc1 = (max(contribution_6nt7_holo$pc1) + min(contribution_6nt7_holo$pc1))/2
contribution_6nt7_holo[c(which(contribution_6nt7_holo$pc1>threshold_6nt7_holo_pc1)),]
threshold_6nt7_holo_pc2 = (max(contribution_6nt7_holo$pc2) + min(contribution_6nt7_holo$pc2))/2
contribution_6nt7_holo[c(which(contribution_6nt7_holo$pc2>threshold_6nt7_holo_pc2)),]
threshold_6nt7_holo_pc3 = (max(contribution_6nt7_holo$pc3) + min(contribution_6nt7_holo$pc3))/2
contribution_6nt7_holo[c(which(contribution_6nt7_holo$pc3>threshold_6nt7_holo_pc3)),]
  
ggplot(contribution_6nt7_holo, aes(x = contribution_6nt7_holo$num)) + 
  geom_line(aes(y = contribution_6nt7_holo$pc1), color = "mediumorchid", size = 0.4)  + ylim(0, 0.4) +
  ylab("Contribution") + xlab("Residue") + ggtitle("FIRST PRINCIPAL COMPONENT") + a +
  geom_line(aes(y = threshold_6nt7_holo_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_holo, aes(x = contribution_6nt7_holo$num)) + 
  geom_line(aes(y = contribution_6nt7_holo$pc2), color = "darkorchid4", size = 0.4) +
  ylab("Contribution") + xlab("Residue") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.25, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_holo, aes(x = contribution_6nt7_holo$num)) + 
  geom_col(aes(y = contribution_6nt7_holo$pc3), color = "mediumpurple1", size = 0.4)   + ylim(0, 0.4) +
  ylab("Contribution") + xlab("Residue") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

#6NT7 APO
ca_ndx_6nt7_apo <- atom.select(pdb_6nt7_apo, elety = "CA") 
pca_6nt7_apo <- pca.xyz(dcd_6nt7_apo[,ca_ndx_6nt7_apo$xyz])
plot(pca_6nt7_apo, col=topo.colors(nrow(dcd_6nt7_apo)))
  
pca_6nt7_apo_db <- data.frame("PC" = pca_6nt7_apo$z)
pca_6nt7_apo_db$time <- c(1:2003)
pca_6nt7_apo_db$time[1] = 100
pca_6nt7_apo_db$time[(2:2003)] = hbonds$TIME

A <- ggplot(pca_6nt7_apo_db, aes(x=pca_6nt7_apo_db$PC.1, y=pca_6nt7_apo_db$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (13.28%)") + xlab("CP1 (30.54%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
B <- ggplot(pca_6nt7_apo_db, aes(x=pca_6nt7_apo_db$PC.3, y=pca_6nt7_apo_db$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (13.28%)") + xlab("CP3 (6.17%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt7_apo_db, aes(x=pca_6nt7_apo_db$PC.1, y=pca_6nt7_apo_db$PC.3)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP3 (6.17%)") + xlab("CP1 (30.54%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
  
pc1_6nt7_apo <- readPNG("pics/cerrada_sin_PC1_restime.png")
pc2_6nt7_apo  <- readPNG("pics/cerrada_sin_PC2_restime.png")
g <- rasterGrob(pc1_6nt7_apo, interpolate=TRUE)
i <- rasterGrob(pc2_6nt7_apo, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(i, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))
  
grid.arrange(arrangeGrob(D, ncol = 1, nrow = 1),
             arrangeGrob(B, ncol = 1, nrow = 1), widths = c(2,3))
  
eigen_bd_6nt7_apo <- data.frame("eigen" = pca_6nt7_apo$L)
sum_eig_6nt7_apo = sum(eigen_bd_6nt7_apo$eigen)
eigen_bd_6nt7_apo$percent <- (eigen_bd_6nt7_apo$eigen/sum_eig_6nt7_apo)*100
eigen_bd_6nt7_apo$num <- c(1:1182)
eigen_bd_6nt7_apo$sum_percent = 0
eigen_bd_6nt7_apo$sum_percent[1] = eigen_bd_6nt7_apo$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_apo$sum_percent[i] = eigen_bd_6nt7_apo$percent[i] + eigen_bd_6nt7_apo$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt7_apo[c(1:7),], aes(x = eigen_bd_6nt7_apo$num[c(1:7)], 
                                           y=eigen_bd_6nt7_apo$percent[c(1:7)])) + 
  geom_line(color="#110BD0", size=0.8) + 
  geom_point(color="#110BD0") +
  geom_text(aes(y= eigen_bd_6nt7_apo$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_apo$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente Principal") + scale_x_continuous(breaks = c(1:7)) + 
  ylab("Acumulated percentaje of variance") + a + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35))
  
#mktrj.pca(pca_6nt7_apo, pc=1, b=pca_6nt7_apo$au[,1], file="6nt7_APO_NO_BORRAR/6nt7_apo_pc1_restime.pdb")
#mktrj.pca(pca_6nt7_apo, pc=2, b=pca_6nt7_apo$au[,2], file="6nt7_APO_NO_BORRAR/6nt7_apo_pc2_restime.pdb")
#mktrj.pca(pca_6nt7_apo, pc=3, b=pca_6nt7_apo$au[,3], file="6nt7_APO_NO_BORRAR/6nt7_apo_pc3_restime.pdb")
  
contribution_6nt7_apo <- data.frame("num" = pdb_6nt7_apo$atom$resno[ca_ndx_6nt7_apo$atom], 
                                    "res" = pdb_6nt7_apo$atom$resid[ca_ndx_6nt7_apo$atom],
                                    "resnum" = pdb_6nt7_apo$atom$resno[ca_ndx_6nt7_apo$atom],
                                    "pc1" = pca_6nt7_apo$au[,1], 
                                    "pc2" = pca_6nt7_apo$au[,2], 
                                    "pc3" = pca_6nt7_apo$au[,3]) 
  
i = 146
for (j in 1:392){
  contribution_6nt7_apo$num[j] = i
  i = i+1
}
  
threshold_6nt7_apo_pc1 = (max(contribution_6nt7_apo$pc1) + min(contribution_6nt7_apo$pc1))/2
contribution_6nt7_apo[c(which(contribution_6nt7_apo$pc1>threshold_6nt7_apo_pc1)),]
threshold_6nt7_apo_pc2 = (max(contribution_6nt7_apo$pc2) + min(contribution_6nt7_apo$pc2))/2
contribution_6nt7_apo[c(which(contribution_6nt7_apo$pc2>threshold_6nt7_apo_pc2)),]
threshold_6nt7_apo_pc3 = (max(contribution_6nt7_apo$pc3) + min(contribution_6nt7_apo$pc3))/2
contribution_6nt7_apo[c(which(contribution_6nt7_apo$pc3>threshold_6nt7_apo_pc3)),]
  
ggplot(contribution_6nt7_apo, aes(x = contribution_6nt7_apo$num)) + 
  geom_col(aes(y = contribution_6nt7_apo$pc1), color = "dodgerblue3", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_apo_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_apo, aes(x = contribution_6nt7_apo$num)) + 
  geom_col(aes(y = contribution_6nt7_apo$pc2), color = "steelblue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("SECOND PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_apo, aes(x = contribution_6nt7_apo$num)) + 
  geom_col(aes(y = contribution_6nt7_apo$pc3), color = "mediumblue", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("THIRD PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc3), color = "black", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.4, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

## STING CERRADA+LIGANDO ROTADO
ca_ndx_6nt7_holo_rot <- atom.select(lig_rot_pdb, elety = "CA")
pca_6nt7_holo_rot <- pca.xyz(lig_rot_dcd[,ca_ndx_6nt7_holo_rot$xyz])
plot(pca_6nt7_holo_rot, col=topo.colors(nrow(lig_rot_dcd)))

pca_6nt7_holo_rot_db <- data.frame("PC" = pca_6nt7_holo_rot$z)
pca_6nt7_holo_rot_db$time <- c(1:2002)
pca_6nt7_holo_rot_db$time[(1:2002)] = hbonds$TIME

ggplot(pca_6nt7_holo_rot_db, aes(x=pca_6nt7_holo_rot_db$PC.1, y=pca_6nt7_holo_rot_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.81%)") + xlab("CP1 (27.9%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt7_holo_rot_db, aes(x=pca_6nt7_holo_rot_db$PC.1, y=pca_6nt7_holo_rot_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP3 (5%)") + xlab("CP1 (27.9%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt7_holo_rot_db, aes(x=pca_6nt7_holo_rot_db$PC.3, y=pca_6nt7_holo_rot_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.81%)") + xlab("CP3 (5%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))

eigen_bd_6nt7_holo_rot <- data.frame("eigen" = pca_6nt7_holo_rot$L)
sum_eig_6nt7_holo_rot = sum(eigen_bd_6nt7_holo_rot$eigen)
eigen_bd_6nt7_holo_rot$percent <- (eigen_bd_6nt7_holo_rot$eigen/sum_eig_6nt7_holo_rot)*100
eigen_bd_6nt7_holo_rot$num <- c(1:1182)
eigen_bd_6nt7_holo_rot$sum_percent = 0
eigen_bd_6nt7_holo_rot$sum_percent[1] = eigen_bd_6nt7_holo_rot$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_holo_rot$sum_percent[i] = eigen_bd_6nt7_holo_rot$percent[i] + eigen_bd_6nt7_holo_rot$sum_percent[i - 1]
}

ggplot(eigen_bd_6nt7_holo_rot[c(1:7),], aes(x = eigen_bd_6nt7_holo_rot$num[c(1:7)], 
                                                 y=eigen_bd_6nt7_holo_rot$percent[c(1:7)])) + 
  geom_line(color="dodgerblue3", size=0.8) + 
  geom_point(color="dodgerblue3") + 
  geom_text(aes(y= eigen_bd_6nt7_holo_rot$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_holo_rot$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente") + 
  ylab("Acumulated percentaje of variance") + 
  a + scale_x_continuous(breaks = c(1:8)) + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

mktrj.pca(pca_6nt7_holo_rot, pc=1, b=pca_6nt7_holo_rot$au[,1], file="~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada_lig_rot_pc1.pdb")
mktrj.pca(pca_6nt7_holo_rot, pc=2, b=pca_6nt7_holo_rot$au[,2], file="~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada_lig_rot_pc2.pdb")
mktrj.pca(pca_6nt7_holo_rot, pc=3, b=pca_6nt7_holo_rot$au[,3], file="~/Desktop/Sysbio/cSTING/rot_cgamp/production/cerrada_lig_rot_pc3.pdb")

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
contribution_6nt7_holo_rot[c(which(contribution_6nt7_holo_rot$pc1>threshold_6nt7_holo_rot_pc1)),]
threshold_6nt7_holo_rot_pc2 = (max(contribution_6nt7_holo_rot$pc2) + min(contribution_6nt7_holo_rot$pc2))/2
contribution_6nt7_holo_rot[c(which(contribution_6nt7_holo_rot$pc2>threshold_6nt7_holo_rot_pc2)),]
threshold_6nt7_holo_rot_pc3 = (max(contribution_6nt7_holo_rot$pc3) + min(contribution_6nt7_holo_rot$pc3))/2
contribution_6nt7_holo_rot[c(which(contribution_6nt7_holo_rot$pc3>threshold_6nt7_holo_rot_pc3)),]

D <- ggplot(contribution_6nt7_holo_rot, aes(x = num)) + 
  geom_line(aes(y = pc1), color = "dodgerblue3", size = 0.4) +
  ylab("Contribution") + xlab("Residue") + a +
  geom_line(aes(y = threshold_6nt7_holo_rot_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

B <- ggplot(contribution_6nt7_holo_rot, aes(x = num)) + 
  geom_col(aes(y = pc2), color = "blue4", size = 0.4)  +
  ylab("Contribution") + xlab("Residue") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

C <- ggplot(contribution_6nt7_holo_rot, aes(x = num)) + 
  geom_col(aes(y = pc3), color = "royalblue3", size = 0.4)   +
  ylab("Contribution") + xlab("Residue") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_rot_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.3, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

pc1_rot  <- readPNG("/home/hoseki/Desktop/Sysbio/cSTING/rot_cgamp/production/PC1.png")
g <- rasterGrob(pc1_rot, interpolate=TRUE)
C <- qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

############################################### 10- PUENTES DE HIDRÓGENO ###################################################
pos <- seq(from = 1, to = 2502, by = 5)
hbonds <- data.frame("time" = hbonds_6nt6$time,
                     "abierta" = hbonds_6nt6$hbonds,
                     "cerrada" = hbonds_6nt7_holo$hbonds)
ggplot(hbonds[pos,], aes(x = time)) + 
  geom_line(aes(y = abierta), color = "mediumseagreen") + 
  geom_line(aes(y = cerrada), color = "mediumorchid2") +
  xlab("Time (ns)") + ylab("Number of Hydrogen Bonds") +
  ggtitle("Number of hydrogen bonds per time") + a
ggplot(hbonds_6nt6[pos,], aes(time, hbonds)) + geom_line() + a


describe(hbonds)

ggplot(hbonds, aes(x = HBOND)) + geom_histogram(binwidth = 0.5, fill = "black") + a +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  ylab("Ocurrencia") + xlab("Cantidad de puentes de hidrógeno") 
  

ggplot(hbonds, aes(x=TIME)) + 
  geom_point(aes(y = hbonds$HBOND), color = "mediumorchid2", size = 0.6) + 
  geom_line(aes(y = hbonds$HBOND), color = "black", size = 0.2, alpha = 0.5) + 
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  a + xlab("Tiempo (ns)") + ylab("Cantidad de puentes de hidrógeno") +
  geom_line(data = df, aes(x= df$Time, y=df$avg), color = "blue", size = 0.8)

df <- data.frame("Time" = c(0:100), "hbond" = 0, "count" = 0)
df$Time = df$Time*5
df$hbond <- as.integer(df$hbond)
j = 1
df <- df[-(1:20),]
for(i in 1: 2002){
  if(hbonds$TIME[i] == df$Time[j]){
    df$hbond[j] = df$hbond[j] + hbonds$HBOND[i]
    df$count[j] = df$count[j] + 1
  }
  else{
    if(hbonds$TIME[i] < df$Time[j]){
      df$hbond[j] = df$hbond[j] + hbonds$HBOND[i]
      df$count[j] = df$count[j] + 1
    }
    else{
      j = j + 1
      df$hbond[j] = df$hbond[j] + hbonds$HBOND[i]
      df$count[j] = df$count[j] + 1
    }
  }
}
df$avg <- df$hbond/df$count
ggplot(df, aes(Time, avg)) + geom_line() + geom_point()

abierta_hbonds <- read.table("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds-protein.dat", quote="\"", comment.char="")
cerrada_hbonds <- read.table("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-protein.dat", quote="\"", comment.char="")
rotada_hbonds <- read.table("~/Desktop/Sysbio/cSTING/rot_cgamp/production/HBONDS/hbonds-protein.dat", quote="\"", comment.char="")

names(abierta_hbonds)[] <- c("time", "hb")
names(cerrada_hbonds)[] <- c("time", "hb")
names(rotada_hbonds)[] <- c("time", "hb")

abierta_hbonds$time <- c(0, rms_6nt6$`Time (ns)`)
cerrada_hbonds$time <- c(0, rms_6nt6$`Time (ns)`)
rotada_hbonds$time <-  c(0, rms_6nt6$`Time (ns)`)

abierta_hbonds$label <- "APO-STING"
cerrada_hbonds$label <- "HOLO-STING"
rotada_hbonds$label <- "HOLO-STING rot"

hbond_df <- data.frame("time" = c(abierta_hbonds$time, cerrada_hbonds$time, rotada_hbonds$time),
                       "hb" = c(abierta_hbonds$hb, cerrada_hbonds$hb, rotada_hbonds$hb),
                       "Structure" = c(abierta_hbonds$label, cerrada_hbonds$label, rotada_hbonds$label)
)
pos <- seq(from = 1, to = 2502, by = 5)
pos2 <- seq(from = 2503, to = 5004, by = 5)
pos3 <- seq(from = 5005, to = 7506, by = 5)
ggplot(hbond_df[c(pos3),], aes(x = time, y = hb)) + 
  geom_line(aes(color = Structure, linetype = Structure)) + 
  scale_color_manual(values = c("mediumseagreen", "mediumorchid2", "dodgerblue3")) +
  xlab("Time (ns)") + ylab("Hydrogen Bonds (nm)") + a 

mean(abierta_hbonds$hb[502:2502])
sd(abierta_hbonds$hb[502:2502])
mean(cerrada_hbonds$hb[502:2502])
sd(cerrada_hbonds$hb[502:2502])
mean(rotada_hbonds$hb[502:2502])
sd(rotada_hbonds$hb[502:2502])
############################################# 11- MAS COSAS DEL SITIO DE POLIMERIZACION ###################################
`pol-bfac-ca` <- read.csv("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-polsite-negative/pol-bfac-ca-rmsd.xvg", sep="")
`pol-bfac` <- read.csv("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-polsite-negative/pol-bfac-rmsd.xvg", sep="")
`pol-ca` <- read.csv("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-polsite-negative/pol-ca-rmsd.xvg", sep="")
pol <- read.csv("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/hbonds-polsite-negative/pol-rmsd.xvg", sep="")
rms_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/500ns_6nt7_holo_NpT_nojump_rot_rmsd.xvg", sep="")
gly_6nt7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/6nt7_holo_gly163_dist.xvg", sep="")

db <- data.frame("time" = pol$TIME,
                 "pol" = pol$RMSD,
                 "pol-ca" = `pol-ca`$RMSD,
                 "pol-bfac" = `pol-bfac`$RMSD,
                 "pol-bfac-ca" = `pol-bfac-ca`$RMSD, 
                 "all" = rms_6nt7_holo$`RMSDns`,
                 "gly" = gly_6nt7_holo$DIS)
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic'),
  panel.background = element_rect(fill = "white", colour = "black")
)

ggplot(db, aes(x = time)) + 
  geom_line(aes(y = pol), color = "blue") + 
  geom_line(aes(y = pol.ca), color = "green") + 
  geom_line(aes(y = pol.bfac), color = "red") + 
  geom_line(aes(y = pol.bfac.ca), color = "orange") + 
  geom_line(aes(y = gly), color = "purple") + 
  a + xlab("Time (ns)") + ylab("nm") 

  dist_ci <- data.frame("time" = com_cerrada$time)
dist_ci$res279 <- sqrt((com_cerrada$x[1] - pol.traj$X279)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y279)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z279)^2)
dist_ci$res280 <- sqrt((com_cerrada$x[1] - pol.traj$X280)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y280)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z280)^2)
dist_ci$res281 <- sqrt((com_cerrada$x[1] - pol.traj$X281)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y281)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z281)^2)
dist_ci$res282 <- sqrt((com_cerrada$x[1] - pol.traj$X282)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y282)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z282)^2)
dist_ci$res283 <- sqrt((com_cerrada$x[1] - pol.traj$X283)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y283)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z283)^2)
dist_ci$res284 <- sqrt((com_cerrada$x[1] - pol.traj$X284)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y284)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z284)^2)
dist_ci$res285 <- sqrt((com_cerrada$x[1] - pol.traj$X285)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y285)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z285)^2)
dist_ci$res286 <- sqrt((com_cerrada$x[1] - pol.traj$X286)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y286)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z286)^2)
dist_ci$res287 <- sqrt((com_cerrada$x[1] - pol.traj$X287)^2 + 
                         (com_cerrada$y[1] - pol.traj$Y287)^2 + 
                         (com_cerrada$z[1] - pol.traj$Z287)^2)

A <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res279), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES279")
B <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res280), color = "darkgreen") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES280")
C <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res281), color = "darkgreen") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES281")
D <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res282), color = "darkgreen") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES282")
E <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res283), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES283")
F <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res284), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES284")
G <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res285), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES285")
H <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res286), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES286")
I <- ggplot(dist_ci, aes(x = time)) + 
  geom_line(aes(y = res287), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES287")
grid.arrange(A,B,C,D,E,F,G,H,I, ncol=3, nrow=3)



dist_sim <- data.frame("time" = com_cerrada$time)
dist_sim$res279 <- sqrt((com_cerrada$x - pol.traj$X279)^2 + 
                          (com_cerrada$y - pol.traj$Y279)^2 + 
                          (com_cerrada$z - pol.traj$Z279)^2)
dist_sim$res280 <- sqrt((com_cerrada$x - pol.traj$X280)^2 + 
                          (com_cerrada$y - pol.traj$Y280)^2 + 
                          (com_cerrada$z - pol.traj$Z280)^2)
dist_sim$res281 <- sqrt((com_cerrada$x - pol.traj$X281)^2 + 
                          (com_cerrada$y - pol.traj$Y281)^2 + 
                          (com_cerrada$z - pol.traj$Z281)^2)
dist_sim$res282 <- sqrt((com_cerrada$x - pol.traj$X282)^2 + 
                          (com_cerrada$y - pol.traj$Y282)^2 + 
                          (com_cerrada$z - pol.traj$Z282)^2)
dist_sim$res283 <- sqrt((com_cerrada$x - pol.traj$X283)^2 + 
                          (com_cerrada$y - pol.traj$Y283)^2 + 
                          (com_cerrada$z - pol.traj$Z283)^2)
dist_sim$res284 <- sqrt((com_cerrada$x - pol.traj$X284)^2 + 
                          (com_cerrada$y - pol.traj$Y284)^2 + 
                          (com_cerrada$z - pol.traj$Z284)^2)
dist_sim$res285 <- sqrt((com_cerrada$x - pol.traj$X285)^2 + 
                          (com_cerrada$y - pol.traj$Y285)^2 + 
                          (com_cerrada$z - pol.traj$Z285)^2)
dist_sim$res286 <- sqrt((com_cerrada$x - pol.traj$X286)^2 + 
                          (com_cerrada$y - pol.traj$Y286)^2 + 
                          (com_cerrada$z - pol.traj$Z286)^2)
dist_sim$res287 <- sqrt((com_cerrada$x - pol.traj$X287)^2 + 
                          (com_cerrada$y - pol.traj$Y287)^2 + 
                          (com_cerrada$z - pol.traj$Z287)^2)

A2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res279), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES279")
B2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res280), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES280")
C2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res281), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES281")
D2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res282), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES282")
E2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res283), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES283")
F2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res284), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES284")
G2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res285), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES285")
H2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res286), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES286")
I2 <- ggplot(dist_sim, aes(x = time)) + 
  geom_line(aes(y = res287), color = "black") + 
  a + ylab("distance") + xlab("time") + ggtitle("RES287")
grid.arrange(A2,B2,C2,D2,E2,F2,G2,H2,I2, ncol=3, nrow=3)

###########################################################################################################################
  
plots  <- list(A,B,C)
do.call(grid.arrange, plots)

grid.arrange(A,B,C,D, ncol=2, nrow=2)
