###########################################################################################################################
###########################################################################################################################
#####                                                                                                                 #####
#####                                                SCRIPT TESINA                                                    #####
#####                                                                                                                 #####
##### 0- Archivos y cosas generales                                                                                   #####
##### 1- Analisis STING ABIERTA/NO_FIT                                                                                #####
##### 2- Analisis STING cerrada con ligando                                                                           #####
##### 3- Analisis STING cerrada sin ligando                                                                           #####
##### 4- Comparaciones                                                                                                #####
##### 5- Distancia de glicinas 163                                                                                    #####
##### 6- RMSD de distintas regiones                                                                                   #####
##### 7- B-factors                                                                                                    #####
##### 8- Analisis de componentes principales (PCA)                                                                    #####
##### 9- Puentes de hidrógeno                                                                                         #####
#####                                                                                                       -hoseki   #####
###########################################################################################################################
###########################################################################################################################

############################################# 0- ARCHIVOS Y COSAS GENERALES   ###############################################
setwd("~/Desktop/Sysbio/cSTING/corridas restringidas")

#THEME FOR PLOTS
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic'),
  panel.background = element_rect(fill = "white", colour = "black")
)

#ARCHIVOS 
rms_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_rmsd.xvg", sep="")
names(rms_6nt6_rest_nofit)[names(rms_6nt6_rest_nofit)=="Time"] <- "Time (ns)"
names(rms_6nt6_rest_nofit)[names(rms_6nt6_rest_nofit)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_rmsf.xvg", sep="")
rmsf_6nt6_rest_nofit <- rmsf_6nt6_rest_nofit[-(c(1:6, 208:213)),]
i=146
for (j in 1:402){
  rmsf_6nt6_rest_nofit$ATOM[j]=i
  i=i+1
}
gr_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_gr.xvg", sep="")
gr_6nt6_rest_nofit$Time = gr_6nt6_rest_nofit$Time/1000
pdb_6nt6_rest_nofit <- read.pdb("ABIERTA/NO_FIT/6nt6_rest_nofit_400ns.pdb")
dcd_6n6_rest_nofit <- read.dcd("ABIERTA/NO_FIT/6nt6_rest_nofit_400ns.dcd")  

rms_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_rmsd.xvg", sep="")
names(rms_6nt7_holo_rest_nofit)[names(rms_6nt7_holo_rest_nofit)=="Time"] <- "Time (ns)"
names(rms_6nt7_holo_rest_nofit)[names(rms_6nt7_holo_rest_nofit)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_holo_rest_nofit$ATOM[j]=i
  i=i+1
}
gr_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_gr.xvg", sep="")
gr_6nt7_holo_rest_nofit$Time = gr_6nt7_holo_rest_nofit$Time/1000
    pdb_6nt7_holo_rest_nofit <- read.pdb("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit.pdb")
dcd_6nt7_holo_rest_nofit <- read.dcd("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit.dcd")

rms_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_rmsd.xvg", sep="")
names(rms_6nt7_apo_rest_nofit)[names(rms_6nt7_apo_rest_nofit)=="Time"] <- "Time (ns)"
names(rms_6nt7_apo_rest_nofit)[names(rms_6nt7_apo_rest_nofit)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_apo_rest_nofit$ATOM[j]=i
  i=i+1
}
gr_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_gr.xvg", sep="")
gr_6nt7_apo_rest_nofit$Time = gr_6nt7_apo_rest_nofit$Time/1000
pdb_6nt7_apo_rest_nofit <- read.pdb("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_400ns.pdb")
dcd_6nt7_apo_rest_nofit <- read.dcd("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_400ns.dcd") 

pol_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_pol_rmsd.xvg", sep="")
pol_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_pol_rmsd.xvg", sep="")
pol_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_pol_rmsd.xvg", sep="")

helix_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_helix_rmsd.xvg", sep="")
helix_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_helix_rmsd.xvg", sep="")
helix_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_helix_rmsd.xvg", sep="")

loop_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_loop_rmsd.xvg", sep="")
loop_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_loop_rmsd.xvg", sep="")
loop_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_loop_rmsd.xvg", sep="")

gly_6nt6_rest_nofit <- read.csv("ABIERTA/NO_FIT/6nt6_rest_nofit_gly.xvg", sep="")
gly_6nt7_apo_rest_nofit <- read.csv("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_gly.xvg", sep="")
gly_6nt7_holo_rest_nofit <- read.csv("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_distgly163.xvg", sep="")

bfac_6nt6_rest_nofit <- read.pdb("ABIERTA/NO_FIT/6nt6_rest_nofit_bfac.pdb")
bfac_6nt7_holo_rest_nofit <- read.pdb("CERRADA_LIGANDO/NO_FIT/6nt7_holo_rest_nofit_bfac.pdb")
bfac_6nt7_apo_rest_nofit <- read.pdb("CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_rest_nofit_bfac.pdb")

hbonds_rest_nofit <- read.table("CERRADA_LIGANDO/NO_FIT/hbonds.dat", quote="\"", comment.char="")
names(hbonds_rest_nofit)[names(hbonds_rest_nofit)=="V1"] <- "TIME"
names(hbonds_rest_nofit)[names(hbonds_rest_nofit)=="V2"] <- "HBOND"
hbonds_rest_nofit$TIME <- rms_6nt7_holo_rest_nofit$`Time (ns)`[c(500:2501)]

############################################## 1- ANALISIS STING ABIERTA/NO_FIT ##################################################
#RMSD
ggplot(rms_6nt6_rest_nofit, aes(x=rms_6nt6_rest_nofit$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt6_rest_nofit$`RMSD (nm)`), color = "yellowgreen", size = 0.8) + 
  ggtitle("RMSD OPEN STING (XYZ RESTRICTION) NO FIT") + 
  xlab("Time (ns)") + 
  ylab("RMSD (nm)") + 
  a + ylim(0, 0.6)

#RMSF
ggplot(rmsf_6nt6_rest_nofit, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "yellowgreen", size=0.8 ) + 
  ggtitle("RMSF OPEN STING (XYZ RESTRICTION) NO FIT") +
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
ggplot(gr_6nt6_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt6_rest_nofit$GR), color="yellowgreen", size=0.8) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (XYZ RESTRICTION) NO FIT") + 
  xlab("Time (ns)") +   ylab("RG (nm)") + a + ylim(2.1, 2.4)

ggplot(gr_6nt6_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt6_rest_nofit$GR), color="yellowgreen", size=0.8) + 
  geom_line(aes(y=gr_6nt6_rest_nofit$GRX), color="darkgreen", size=0.7) + 
  geom_line(aes(y=gr_6nt6_rest_nofit$GRY), color="magenta2", size=0.7) + 
  geom_line(aes(y=gr_6nt6_rest_nofit$GRZ), color="orange", size=0.7) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (XYZ RESTRICTION) NO FIT") + xlab("Time (ns)") + 
  ylab("GR (nm)") + a +
  annotate("text", x=28, y=1.54, label="X-RG", color="darkgreen") +
  annotate("text", x=490, y=1.6, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.15, label="Y-RG", color="magenta2")


#70-105
#280-400
#400-500
############################################# 2- ANALISIS STING CERRADA CON LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_holo_rest_nofit, aes(x=rms_6nt7_holo_rest_nofit$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_holo_rest_nofit$`RMSD (nm)`), color = "slateblue4", size=0.8) + 
  ggtitle("RMSD CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") + ylim(0,0.6)+
  xlab("Time (ns)") + ylab("RMSD (nm)") + a 

#RMSF
ggplot(rmsf_6nt7_holo_rest_nofit, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "slateblue4", size=0.8) + 
  ggtitle("RMSF CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") + 
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
ggplot(gr_6nt7_holo_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo_rest_nofit$GR), color="slateblue4", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") + ylim(2.1,2.4) +
  xlab("Time (ns)") + ylab("RG (nm)") + a

ggplot(gr_6nt7_holo_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo_rest_nofit$GR), color="slateblue4", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest_nofit$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest_nofit$GRY), color="blue", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest_nofit$GRZ), color="darkgreen", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") + 
  xlab("Time (ns)") + ylab("RG (nm)") + a +
  annotate("text", x=25, y=1.55, label="X-RG", color="#DC143C") +
  annotate("text", x=480, y=1.9, label="Z-RG", color="darkgreen") +
  annotate("text", x=22, y=2.1, label="Y-RG", color="blue")
############################################# 3- ANALISIS STING CERRADA SIN LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_apo_rest_nofit, aes(x=rms_6nt7_apo_rest_nofit$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_apo_rest_nofit$`RMSD (nm)`), color = "royalblue", size=0.8) + 
  ggtitle("RMSD CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") + a + ylim(0,0.6) +
  xlab("Time (ns)") + ylab("RMSD (nm)") 

#RMSF
ggplot(rmsf_6nt7_apo_rest_nofit, aes(x=ATOM, y=nm)) + 
  geom_line(color = "royalblue", size=0.8) + 
  ggtitle("RMSD CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") + 
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
ggplot(gr_6nt7_apo_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo_rest_nofit$GR), color="royalblue", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") + a +
  xlab("Time (ns)") + ylab("RG (nm)") + ylim(2.1,2.4) 

ggplot(gr_6nt7_apo_rest_nofit, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo_rest_nofit$GR), color="royalblue", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest_nofit$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest_nofit$GRY), color="magenta3", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest_nofit$GRZ), color="orange", size=0.8) + 
  xlab("Time (ns)") + ylab("RG (nm)") + 
  ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") +
  annotate("text", x=1, y=1.5, label="X-RG", color="#DC143C") +
  annotate("text", x=480, y=1.5, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.06, label="Y-RG", color="magenta3") + a
################################################## 4- COMPARACIONES #######################################################
#RMSD
RMSD_rest_nofit <- data.frame("Time(ns)" = rms_6nt6_rest_nofit$`Time (ns)`, 
                        "6NT6" = rms_6nt6_rest_nofit$`RMSD (nm)`, 
                        "6NT7 APO" = rms_6nt7_apo_rest_nofit$`RMSD (nm)`, 
                        "6NT7 HOLO" = rms_6nt7_holo_rest_nofit$`RMSD (nm)`)
ggplot(RMSD_rest_nofit, aes(x=RMSD_rest_nofit$Time.ns.)) + 
  geom_line(aes(y=RMSD_rest_nofit$X6NT6), color = "yellowgreen", size=0.8) + 
  geom_line(aes(y=RMSD_rest_nofit$X6NT7.APO), color = "royalblue", size=0.8) + 
  geom_line(aes(y=RMSD_rest_nofit$X6NT7.HOLO), color = "slateblue4", size=0.8) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + ggtitle("Comparación de RMSD")

#RMSF
d <- c(343, 344, 345, 346, 544, 545, 546, 547)
b <- c(1:8)

for (l in 1:8){
  b[l] <- which(rmsf_6nt6_rest_nofit$ATOM == d[l])
}

apo_rest_nofit <- rmsf_6nt6_rest_nofit[c(-b),]
i=146

for (j in 1:394){
  apo_rest_nofit$ATOM[j]=i
  i=i+1
}

RMSF_rest_nofit <- data.frame("aa" = apo_rest_nofit$ATOM, 
                        "6nt6" = apo_rest_nofit$nm, 
                        "6nt7-apo" = rmsf_6nt7_apo_rest_nofit$nm, 
                        "6nt7-holo"= rmsf_6nt7_holo_rest_nofit$nm)
ggplot(RMSF_rest_nofit, aes(x=RMSF_rest_nofit$aa)) +
  geom_line(aes(y=RMSF_rest_nofit$X6nt6), color = 'yellowgreen', size=0.8) +
  geom_line(aes(y=RMSF_rest_nofit$X6nt7.apo), color ='royalblue', size=0.8) +
  geom_line(aes(y=RMSF_rest_nofit$X6nt7.holo), color ='slateblue4', size=0.8) +
  ylab('RMSF (nm)') + xlab('aminoacid') + a +
  ggtitle('Comparación RSMF')

#RADIO DE GIRO
GR_rest_nofit <- data.frame("Time (ns)"= gr_6nt6_rest_nofit$Time, 
                      "6NT6"=gr_6nt6_rest_nofit$GR, 
                      "6NT7-APO"= gr_6nt7_apo_rest_nofit$GR, 
                      "6NT7-HOLO"= gr_6nt7_holo_rest_nofit$GR)
ggplot(GR_rest_nofit, aes(x= GR_rest_nofit$Time..ns.)) + 
  geom_line(aes(y=GR_rest_nofit$X6NT6), color="yellowgreen", size=0.8) + 
  geom_line(aes(y=GR_rest_nofit$X6NT7.APO), color="royalblue", size=0.8) + 
  geom_line(aes(y=GR_rest_nofit$X6NT7.HOLO), color="slateblue4", size=0.8) +
  ggtitle("Comparación Radio de Giro") +
  xlab("Time (ns)") + ylab("GR (nm)") + a

############################################### 5- DISTANCIA DE GLICINAS 163 ##############################################
ggplot(gly_6nt6_rest_nofit, aes(x=gly_6nt6_rest_nofit$TIME)) + 
  geom_line(aes(y=gly_6nt6_rest_nofit$DIS), color = "yellowgreen", size=0.8, alpha = .7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a + ylim(0.35, 1.1) +
  ggtitle("GLYCINE 163 DISTANCE OPEN STING (XYZ RESTRICTION) NO FIT") +
  geom_line(data = GLY_AVG_rest_nofit, aes(x = Time, y = avg_6NT6), color="darkgreen")

ggplot(gly_6nt7_apo_rest_nofit, aes(x=gly_6nt7_apo_rest_nofit$TIME)) + 
  geom_line(aes(y=gly_6nt7_apo_rest_nofit$DIS), color = "royalblue", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a + ylim(0.35, 1.1) +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") +
  geom_line(data = GLY_AVG_rest_nofit, aes(x = Time, y = avg_APO), color="blue")

ggplot(gly_6nt7_holo_rest_nofit, aes(x=gly_6nt7_holo_rest_nofit$TIME)) + 
  geom_line(aes(y=gly_6nt7_holo_rest_nofit$DIS), color = "slateblue4", size=0.8, alpha = 0.65) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a + ylim(0.35, 1.1) +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") +
  geom_line(data = GLY_AVG_rest_nofit, aes(x = Time, y = avg_HOLO), color="slateblue3")

GLY <- data.frame("TIME" = gly_6nt6_rest_nofit$TIME, 
                  "6NT6" = gly_6nt6_rest_nofit$DIS, 
                  "6NT7 APO" = gly_6nt7_apo_rest_nofit$DIS, 
                  "6NT7 HOLO" = gly_6nt7_holo_rest_nofit$DIS)

GLY_AVG_rest_nofit <- data.frame("Time" = c(0:50), "GLY_6NT6" = 0, "GLY_HOLO"=0, "GLY_APO"=0, "count_6NT6" = 0, "count_APO" = 0, "count_HOLO" = 0)
GLY_AVG_rest_nofit$Time = GLY_AVG_rest_nofit$Time*10
GLY_AVG_rest_nofit$GLY_6NT6 <- as.integer(GLY_AVG_rest_nofit$GLY_6NT6)
GLY_AVG_rest_nofit$GLY_HOLO <- as.integer(GLY_AVG_rest_nofit$GLY_HOLO)
GLY_AVG_rest_nofit$GLY_APO <- as.integer(GLY_AVG_rest_nofit$GLY_APO)
j = 1
for(i in 1: 2501){
  if(gly_6nt6_rest_nofit$TIME[i] == GLY_AVG_rest_nofit$Time[j]){
    GLY_AVG_rest_nofit$GLY_6NT6[j] = GLY_AVG_rest_nofit$GLY_6NT6[j] + gly_6nt6_rest_nofit$DIS[i]
    GLY_AVG_rest_nofit$count_6NT6[j] = GLY_AVG_rest_nofit$count_6NT6[j] + 1
  }
  else{
    if(gly_6nt6_rest_nofit$TIME[i] <GLY_AVG_rest_nofit$Time[j]){
      GLY_AVG_rest_nofit$GLY_6NT6[j] = GLY_AVG_rest_nofit$GLY_6NT6[j] + gly_6nt6_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_6NT6[j] = GLY_AVG_rest_nofit$count_6NT6[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest_nofit$GLY_6NT6[j] = GLY_AVG_rest_nofit$GLY_6NT6[j] + gly_6nt6_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_6NT6[j] = GLY_AVG_rest_nofit$count_6NT6[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_holo_rest_nofit$TIME[i] == GLY_AVG_rest_nofit$Time[j]){
    GLY_AVG_rest_nofit$GLY_HOLO[j] = GLY_AVG_rest_nofit$GLY_HOLO[j] + gly_6nt7_holo_rest_nofit$DIS[i]
    GLY_AVG_rest_nofit$count_HOLO[j] = GLY_AVG_rest_nofit$count_HOLO[j] + 1
  }
  else{
    if(gly_6nt7_holo_rest_nofit$TIME[i] <GLY_AVG_rest_nofit$Time[j]){
      GLY_AVG_rest_nofit$GLY_HOLO[j] = GLY_AVG_rest_nofit$GLY_HOLO[j] + gly_6nt7_holo_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_HOLO[j] = GLY_AVG_rest_nofit$count_HOLO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest_nofit$GLY_HOLO[j] = GLY_AVG_rest_nofit$GLY_HOLO[j] + gly_6nt7_holo_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_HOLO[j] = GLY_AVG_rest_nofit$count_HOLO[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_apo_rest_nofit$TIME[i] == GLY_AVG_rest_nofit$Time[j]){
    GLY_AVG_rest_nofit$GLY_APO[j] = GLY_AVG_rest_nofit$GLY_APO[j] + gly_6nt7_apo_rest_nofit$DIS[i]
    GLY_AVG_rest_nofit$count_APO[j] = GLY_AVG_rest_nofit$count_APO[j] + 1
  }
  else{
    if(gly_6nt7_apo_rest_nofit$TIME[i] <GLY_AVG_rest_nofit$Time[j]){
      GLY_AVG_rest_nofit$GLY_APO[j] = GLY_AVG_rest_nofit$GLY_APO[j] + gly_6nt7_apo_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_APO[j] = GLY_AVG_rest_nofit$count_APO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest_nofit$GLY_APO[j] = GLY_AVG_rest_nofit$GLY_APO[j] + gly_6nt7_apo_rest_nofit$DIS[i]
      GLY_AVG_rest_nofit$count_APO[j] = GLY_AVG_rest_nofit$count_APO[j] + 1
    }
  }
}



GLY_AVG_rest_nofit$avg_6NT6 <- GLY_AVG_rest_nofit$GLY_6NT6/GLY_AVG_rest_nofit$count_6NT6
GLY_AVG_rest_nofit$avg_HOLO <- GLY_AVG_rest_nofit$GLY_HOLO/GLY_AVG_rest_nofit$count_HOLO
GLY_AVG_rest_nofit$avg_APO <- GLY_AVG_rest_nofit$GLY_APO/GLY_AVG_rest_nofit$count_APO


ggplot(GLY_rest_nofit, aes(x=GLY_rest_nofit$TIME)) +
  geom_line(aes(y= GLY_rest_nofit$X6NT6), color = "yellowgreen", size=0.8) + 
  geom_line(aes(y= GLY_rest_nofit$X6NT7.HOLO), color = "slateblue4", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") + a 

ggplot(GLY_rest_nofit, aes(x=GLY_rest_nofit$TIME)) +
  geom_line(aes(y= GLY_rest_nofit$X6NT6), color = "yellowgreen", size=0.8, alpha = 0.6) + 
  geom_line(aes(y= GLY_rest_nofit$X6NT7.HOLO), color = "slateblue4", size=0.8, alpha=0.6) +
  geom_line(aes(y= GLY_rest_nofit$X6NT7.APO), color = "royalblue", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") + a 

############################################### 6- RMSD DE DISTINTOS DOMINIOS #############################################
#SITIO DE POLIMERIZACION
ggplot(pol_6nt6_rest_nofit, aes(x=pol_6nt6_rest_nofit$Time)) + 
  geom_line(aes(y=pol_6nt6_rest_nofit$RMSD), color = "yellowgreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.5) +
  ggtitle("RMSD POLYMERIZATION DOMAIN OPEN STING (XYZ RESTRICTION) NO FIT")

ggplot(pol_6nt7_apo_rest_nofit, aes(x = Time)) + 
  geom_line(aes(y = RMSD), color = "royalblue", size = 0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.5) +
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT")

ggplot(pol_6nt7_holo_rest_nofit, aes(x=pol_6nt7_holo_rest_nofit$Time)) + 
  geom_line(aes(y=pol_6nt7_holo_rest_nofit$RMSD), color = "slateblue4", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  +  ylim(0, 0.5) +
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT")

#HELICE CONTECTORA
ggplot(helix_6nt6_rest_nofit, aes(x=helix_6nt6_rest_nofit$Time)) + 
  geom_line(aes(y=helix_6nt6_rest_nofit$RMSD), color = "yellowgreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.75) + 
  ggtitle("RMSD CONECTOR HELIX OPEN STING (XYZ RESTRICTION) NO FIT")

ggplot(helix_6nt7_apo_rest_nofit, aes(x=helix_6nt7_apo_rest_nofit$Time)) + 
  geom_line(aes(y=helix_6nt7_apo_rest_nofit$RMSD), color = "royalblue", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.75)  + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT")

ggplot(helix_6nt7_holo_rest_nofit, aes(x=helix_6nt7_holo_rest_nofit$Time)) + 
  geom_line(aes(y=helix_6nt7_holo_rest_nofit$RMSD), color = "slateblue4", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.75) + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT")

#LOOP CONECTOR
ggplot(loop_6nt6_rest_nofit, aes(x=loop_6nt6_rest_nofit$Time)) + 
  geom_line(aes(y=loop_6nt6_rest_nofit$RMSD), color = "yellowgreen", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.45) +
  ggtitle("RMSD CONECTOR LOOP OPEN STING (XYZ RESTRICTION) NO FIT")

ggplot(loop_6nt7_apo_rest_nofit, aes(x=loop_6nt7_apo_rest_nofit$Time)) + 
  geom_line(aes(y=loop_6nt7_apo_rest_nofit$RMSD), color = "royalblue", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.45) +  
  ggtitle("RMSD CONECTOR LOOP CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT")

ggplot(loop_6nt7_holo_rest_nofit, aes(x=loop_6nt7_holo_rest_nofit$Time)) + 
  geom_line(aes(y=loop_6nt7_holo_rest_nofit$RMSD), color = "slateblue4", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.45) + 
  ggtitle("RMSD CONECTOR LOOP CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT")

#################################################### 7- B-FACTORS #########################################################
#6nt6
bfac_6nt6_db_rest_nofit <- data.frame("NUM" = bfac_6nt6_rest_nofit$atom$resno, 
                                "RES" = bfac_6nt6_rest_nofit$atom$resid, 
                                "BFAC" = bfac_6nt6_rest_nofit$atom$b)
bfac_6nt6_db_rest_nofit <- bfac_6nt6_db_rest_nofit[-c(198:201, 399:402),]
prom_6nt6_rest_nofit <- mean(bfac_6nt6_db_rest_nofit$BFAC)
sd_6nt6_rest_nofit <- sd(bfac_6nt6_db_rest_nofit$BFAC)
bfac_6nt6_db_rest_nofit$NBFAC <- (bfac_6nt6_db_rest_nofit$BFAC - prom_6nt6_rest_nofit)/sd_6nt6_rest_nofit
bfac_6nt6_db_rest_nofit$RMSF <- rmsf_6nt6_rest_nofit$nm


j = 146
for (i in 1:394){
  bfac_6nt6_db_rest_nofit$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt6_db_rest_nofit, aes(x = bfac_6nt6_db_rest_nofit$ATOM)) +
  geom_line(aes(y = bfac_6nt6_db_rest_nofit$NBFAC), color = "yellowgreen", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-factor") + 
  ggtitle("NORMALIZED B-FACTORS OPEN STING (XYZ RESTRICTION) NO FIT") +
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
  annotate("rect", xmin = 189, xmax=200, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

bfac_6nt6_75_rest_nofit <- which(bfac_6nt6_db_rest_nofit$NBFAC > 0.75)
bfac_6nt6_db_rest_nofit[bfac_6nt6_75_rest_nofit,]
bfac_6nt6_50_rest_nofit <- which(bfac_6nt6_db_rest_nofit$NBFAC > 0.5)
bfac_6nt6_db_rest_nofit[bfac_6nt6_50_rest_nofit,]

#6nt7 holo
bfac_6nt7_holo_db_rest_nofit <- data.frame("NUM" = bfac_6nt7_holo_rest_nofit$atom$resno, 
                                     "RES" = bfac_6nt7_holo_rest_nofit$atom$resid, 
                                     "BFAC" = bfac_6nt7_holo_rest_nofit$atom$b)
prom_6nt7_holo_rest_nofit <- mean(bfac_6nt7_holo_db_rest_nofit$BFAC)
sd_6nt7_holo_rest_nofit <- sd(bfac_6nt7_holo_db_rest_nofit$BFAC)
bfac_6nt7_holo_db_rest_nofit$NBFAC <- (bfac_6nt7_holo_db_rest_nofit$BFAC - prom_6nt7_holo_rest_nofit)/sd_6nt7_holo_rest_nofit
bfac_6nt7_holo_db_rest_nofit$RMSF <- rmsf_6nt7_holo_rest_nofit$nm

j = 146
for (i in 1:394){
  bfac_6nt7_holo_db_rest_nofit$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt7_holo_db_rest_nofit, aes(x = bfac_6nt7_holo_db_rest_nofit$ATOM)) +
  geom_line(aes(y = bfac_6nt7_holo_db_rest_nofit$NBFAC), color = "slateblue4", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-factor") + 
  ggtitle("NORMALIZED B-FACTORS CLOSE STING + LIGAND (XYZ RESTRICTION) NO FIT") + 
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
  annotate("rect", xmin = 191, xmax=196, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

bfac_6nt7_holo_75_rest_nofit <- which(bfac_6nt7_holo_db_rest_nofit$NBFAC > 0.75)
bfac_6nt7_holo_db_rest_nofit[bfac_6nt7_holo_75_rest_nofit,]
bfac_6nt7_holo_50_rest_nofit <- which(bfac_6nt7_holo_db_rest_nofit$NBFAC > 0.5)
bfac_6nt7_holo_db_rest_nofit[bfac_6nt7_holo_50_rest_nofit,]

#6nt7 apo
bfac_6nt7_apo_db_rest_nofit <- data.frame("NUM" = bfac_6nt7_apo_rest_nofit$atom$resno, 
                                    "RES" = bfac_6nt7_apo_rest_nofit$atom$resid, 
                                    "BFAC" = bfac_6nt7_apo_rest_nofit$atom$o) 
bfac_6nt7_apo_db_rest_nofit$BFAC <- bfac_6nt7_apo_db_rest_nofit$BFAC + bfac_6nt7_apo_rest_nofit$atom$b
prom_6nt7_apo_rest_nofit <- mean(bfac_6nt7_apo_db_rest_nofit$BFAC)
sd_6nt7_apo_rest_nofit <- sd(bfac_6nt7_apo_db_rest_nofit$BFAC)
bfac_6nt7_apo_db_rest_nofit$NBFAC <- (bfac_6nt7_apo_db_rest_nofit$BFAC - prom_6nt7_apo_rest_nofit)/sd_6nt7_apo_rest_nofit
bfac_6nt7_apo_db_rest_nofit$RMSF <- rmsf_6nt7_apo_rest_nofit$nm

j = 146
for (i in 1:394){
  bfac_6nt7_apo_db_rest_nofit$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt7_apo_db_rest_nofit, aes(x = bfac_6nt7_apo_db_rest_nofit$ATOM)) +
  geom_line(aes(y = bfac_6nt7_apo_db_rest_nofit$NBFAC), color = "royalblue", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-Factor") + 
  ggtitle("NORMALIZED B-FACTORS CLOSE STING NO LIGAND (XYZ RESTRICTION) NO FIT") +
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
  annotate("rect", xmin = 191, xmax=196, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=-0.6, ymax= 12.6, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

bfac_6nt7_apo_75_rest_nofit <- which(bfac_6nt7_apo_db_rest_nofit$NBFAC > 0.75)
bfac_6nt7_apo_db_rest_nofit[bfac_6nt7_apo_75_rest_nofit,]
bfac_6nt7_apo_50_rest_nofit <- which(bfac_6nt7_apo_db_rest_nofit$NBFAC > 0.5)
bfac_6nt7_apo_db_rest_nofit[bfac_6nt7_apo_50_rest_nofit,]

ggplot() +
  geom_line(data = bfac_6nt7_apo_db_rest_nofit, aes(x = bfac_6nt7_apo_db_rest_nofit$ATOM, y = bfac_6nt7_apo_db_rest_nofit$NBFAC), color = "royalblue", size = 0.8) +
  geom_line(data = bfac_6nt7_holo_db_rest_nofit, aes(x= bfac_6nt7_holo_db_rest_nofit$ATOM, y= bfac_6nt7_holo_db_rest_nofit$NBFAC), color = "slateblue4", size = 0.8) +
  geom_line(data = bfac_6nt6_db_rest_nofit, aes(x = bfac_6nt6_db_rest_nofit$ATOM, y = bfac_6nt6_db_rest_nofit$NBFAC), color = "yellowgreen", size = 0.8 ) +
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

################################## 8- ANALISIS DE COMPONENTES PRINCIPALES (PCA) ###########################################
#6NT6
ca_ndx_6nt6_rest_nofit <- atom.select(pdb_6nt6_rest_nofit, elety = "CA") 
pca_6nt6_rest_nofit <- pca.xyz(dcd_6n6_rest_nofit[,ca_ndx_6nt6_rest_nofit$xyz])
plot(pca_6nt6_rest_nofit, col=topo.colors(nrow(dcd_6n6_rest_nofit)))

pca_6nt6_db_rest_nofit <- data.frame("PC" = pca_6nt6_rest_nofit$z)
pca_6nt6_db_rest_nofit$time <- c(1:2002)
pca_6nt6_db_rest_nofit$time[(1:2002)] = rms_6nt6_rest_nofit$`Time (ns)`[c(500:2501)]

A <- ggplot(pca_6nt6_db_rest_nofit, aes(x=pca_6nt6_db_rest_nofit$PC.1, y=pca_6nt6_db_rest_nofit$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (20.97%)") + xlab("CP1 (60.59%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
B <- ggplot(pca_6nt6_db_rest_nofit, aes(x=pca_6nt6_db_rest_nofit$PC.3, y=pca_6nt6_db_rest_nofit$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (20.97%)") + xlab("CP3 (6.84%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt6_db_rest_nofit, aes(x=pca_6nt6_db_rest_nofit$PC.1, y=pca_6nt6_db_rest_nofit$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (60.59%)") + ylab("CP3 (6.84%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

eigen_bd_6nt6_rest_nofit <- data.frame("eigen" = pca_6nt6_rest_nofit$L)
sum_eig_6nt6_rest_nofit = sum(eigen_bd_6nt6_rest_nofit$eigen)
eigen_bd_6nt6_rest_nofit$percent_6nt6 <- (eigen_bd_6nt6_rest_nofit$eigen/sum_eig_6nt6_rest_nofit)*100
eigen_bd_6nt6_rest_nofit$num <- c(1:1206)
eigen_bd_6nt6_rest_nofit$sum_percent = 0
eigen_bd_6nt6_rest_nofit$sum_percent[1] = eigen_bd_6nt6_rest_nofit$percent_6nt6[1]

for (i in 2:1206){
  eigen_bd_6nt6_rest_nofit$sum_percent[i] = eigen_bd_6nt6_rest_nofit$percent_6nt6[i] + eigen_bd_6nt6_rest_nofit$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt6_rest_nofit[c(1:7),], aes(x = eigen_bd_6nt6_rest_nofit$num[c(1:7)], y=eigen_bd_6nt6_rest_nofit$percent_6nt6[c(1:7)])) + 
  geom_line(color="yellowgreen", size=0.8) + 
  geom_point(color="yellowgreen") + 
  geom_text(aes(y= eigen_bd_6nt6_rest_nofit$percent_6nt6[c(1:7)] + 1, 
                label = round(eigen_bd_6nt6_rest_nofit$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

pc1_6nt6_rest_nofit <- readPNG("pics/6nt6_pc1_rest_nofitime.png")
g <- rasterGrob(pc1_6nt6_rest_nofit, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

grid.arrange(arrangeGrob(B, ncol = 1, nrow = 1),
             arrangeGrob(A, ncol = 1, nrow = 1), widths = c(2,3))

mktrj.pca(pca_6nt6_rest_nofit, pc=1, b=pca_6nt6_rest_nofit$au[,1], file="ABIERTA/NO_FIT/6nt6_pc1_rest_nofit.pdb")
mktrj.pca(pca_6nt6_rest_nofit, pc=2, b=pca_6nt6_rest_nofit$au[,2], file="ABIERTA/NO_FIT/6nt6_pc2_rest_nofit.pdb")
mktrj.pca(pca_6nt6_rest_nofit, pc=3, b=pca_6nt6_rest_nofit$au[,3], file="ABIERTA/NO_FIT/6nt6_pc3_rest_nofit.pdb")

contribution_6nt6_rest_nofit <- data.frame("num" = pdb_6nt6_rest_nofit$atom$resno[ca_ndx_6nt6_rest_nofit$atom],
                                     "resnum" = pdb_6nt6_rest_nofit$atom$resno[ca_ndx_6nt6_rest_nofit$atom],
                                     "res" = pdb_6nt6_rest_nofit$atom$resid[ca_ndx_6nt6_rest_nofit$atom], 
                                     "pc1" = pca_6nt6_rest_nofit$au[,1], 
                                     "pc2" = pca_6nt6_rest_nofit$au[,2], 
                                     "pc3" = pca_6nt6_rest_nofit$au[,3]) 
i = 146
for (j in 1:402){
  contribution_6nt6_rest_nofit$num[j] = i
  i = i+1
}

#threshold (max-min)/2
threshold_6nt6_pc1_rest_nofit = (max(contribution_6nt6_rest_nofit$pc1) + min(contribution_6nt6_rest_nofit$pc1))/2
contribution_6nt6_rest_nofit[c(which(contribution_6nt6_rest_nofit$pc1>threshold_6nt6_pc1_rest_nofit)),]
threshold_6nt6_pc2_rest_nofit = (max(contribution_6nt6_rest_nofit$pc2) + min(contribution_6nt6_rest_nofit$pc2))/2
contribution_6nt6_rest_nofit[c(which(contribution_6nt6_rest_nofit$pc2>threshold_6nt6_pc2_rest_nofit)),]
threshold_6nt6_pc3_rest_nofit = (max(contribution_6nt6_rest_nofit$pc3) + min(contribution_6nt6_rest_nofit$pc3))/2
contribution_6nt6_rest_nofit[c(which(contribution_6nt6_rest_nofit$pc3>threshold_6nt6_pc3_rest_nofit)),]

ggplot(contribution_6nt6_rest_nofit, aes(x = contribution_6nt6_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt6_rest_nofit$pc1), color = "chartreuse3", size = 0.4)  + 
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc1_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)


ggplot(contribution_6nt6_rest_nofit, aes(x = contribution_6nt6_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt6_rest_nofit$pc2), color = "darkolivegreen3", size = 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc2_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt6_rest_nofit, aes(x = contribution_6nt6_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt6_rest_nofit$pc3), color = "springgreen1", size = 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc3_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.1, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

#6NT7 HOLO
ca_ndx_6nt7_holo_rest_nofit <- atom.select(pdb_6nt7_holo_rest_nofit, elety = "CA")
pca_6nt7_holo_rest_nofit <- pca.xyz(dcd_6nt7_holo_rest_nofit[,ca_ndx_6nt7_holo_rest_nofit$xyz])
plot(pca_6nt7_holo_rest_nofit, col=topo.colors(nrow(dcd_6nt7_holo_rest_nofit)))

pca_6nt7_holo_db_rest_nofit <- data.frame("PC" = pca_6nt7_holo_rest_nofit$z)
pca_6nt7_holo_db_rest_nofit$time <- c(1:2002)
pca_6nt7_holo_db_rest_nofit$time[(1:2002)] = hbonds_rest_nofit$TIME

A <- ggplot(pca_6nt7_holo_db_rest_nofit, aes(x=pca_6nt7_holo_db_rest_nofit$PC.1, y=pca_6nt7_holo_db_rest_nofit$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (4.09%)") + xlab("CP1 (91.64%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
B <- ggplot(pca_6nt7_holo_db_rest_nofit, aes(x=pca_6nt7_holo_db_rest_nofit$PC.1, y=pca_6nt7_holo_db_rest_nofit$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP3 (1.79%)") + xlab("CP1 (91.64%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt7_holo_db_rest_nofit, aes(x=pca_6nt7_holo_db_rest_nofit$PC.3, y=pca_6nt7_holo_db_rest_nofit$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (4.09%)") + xlab("CP3 (1.79%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))

pc1_6nt7_holo_rest_nofit <- readPNG("pics/6nt7_holo_pc1_rest_nofitime.png")
g <- rasterGrob(pc1_6nt7_holo_rest_nofit, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

pc2_6nt7_holo_rest_nofit <- readPNG("pics/6nt7_holo_pc2_rest_nofitime.png")
h <- rasterGrob(pc2_6nt7_holo_rest_nofit, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(h, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

grid.arrange(arrangeGrob(A, ncol = 1, nrow = 1),
             arrangeGrob(B, ncol = 1, nrow = 1), widths = c(2,3))

eigen_bd_6nt7_holo_rest_nofit <- data.frame("eigen" = pca_6nt7_holo_rest_nofit$L)
sum_eig_6nt7_holo_rest_nofit = sum(eigen_bd_6nt7_holo_rest_nofit$eigen)
eigen_bd_6nt7_holo_rest_nofit$percent <- (eigen_bd_6nt7_holo_rest_nofit$eigen/sum_eig_6nt7_holo_rest_nofit)*100
eigen_bd_6nt7_holo_rest_nofit$num <- c(1:1182)
eigen_bd_6nt7_holo_rest_nofit$sum_percent = 0
eigen_bd_6nt7_holo_rest_nofit$sum_percent[1] = eigen_bd_6nt7_holo_rest_nofit$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_holo_rest_nofit$sum_percent[i] = eigen_bd_6nt7_holo_rest_nofit$percent[i] + eigen_bd_6nt7_holo_rest_nofit$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt7_holo_rest_nofit[c(1:7),], aes(x = eigen_bd_6nt7_holo_rest_nofit$num[c(1:7)], 
                                                  y=eigen_bd_6nt7_holo_rest_nofit$percent[c(1:7)])) + 
  geom_line(color="slateblue4", size=0.8) + 
  geom_point(color="slateblue4") + 
  geom_text(aes(y= eigen_bd_6nt7_holo_rest_nofit$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_holo_rest_nofit$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente") + 
  ylab("Acumulated percentaje of variance") + 
  a + scale_x_continuous(breaks = c(1:8)) + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

mktrj.pca(pca_6nt7_holo_rest_nofit, pc=1, b=pca_6nt7_holo_rest_nofit$au[,1], file="CERRADA_LIGANDO/NO_FIT/6nt7_holo_pc1_rest_nofit.pdb")
mktrj.pca(pca_6nt7_holo_rest_nofit, pc=2, b=pca_6nt7_holo_rest_nofit$au[,2], file="CERRADA_LIGANDO/NO_FIT/6nt7_holo_pc2_rest_nofit.pdb")
mktrj.pca(pca_6nt7_holo_rest_nofit, pc=3, b=pca_6nt7_holo_rest_nofit$au[,3], file="CERRADA_LIGANDO/NO_FIT/6nt7_holo_pc3_rest_nofit.pdb")

contribution_6nt7_holo_rest_nofit <- data.frame("num" = pdb_6nt7_holo_rest_nofit$atom$resno[ca_ndx_6nt7_holo_rest_nofit$atom],
                                          "resnum" = pdb_6nt7_holo_rest_nofit$atom$resno[ca_ndx_6nt7_holo_rest_nofit$atom],
                                          "res" = pdb_6nt7_holo_rest_nofit$atom$resid[ca_ndx_6nt7_holo_rest_nofit$atom], 
                                          "pc1" = pca_6nt7_holo_rest_nofit$au[,1], 
                                          "pc2" = pca_6nt7_holo_rest_nofit$au[,2], 
                                          "pc3" = pca_6nt7_holo_rest_nofit$au[,3]) 

i = 146
for (j in 1:392){
  contribution_6nt7_holo_rest_nofit$num[j] = i
  i = i+1
}

threshold_6nt7_holo_pc1_rest_nofit = (max(contribution_6nt7_holo_rest_nofit$pc1) + min(contribution_6nt7_holo_rest_nofit$pc1))/2
contribution_6nt7_holo_rest_nofit[c(which(contribution_6nt7_holo_rest_nofit$pc1>threshold_6nt7_holo_pc1_rest_nofit)),]
threshold_6nt7_holo_pc2_rest_nofit = (max(contribution_6nt7_holo_rest_nofit$pc2) + min(contribution_6nt7_holo_rest_nofit$pc2))/2
contribution_6nt7_holo_rest_nofit[c(which(contribution_6nt7_holo_rest_nofit$pc2>threshold_6nt7_holo_pc2_rest_nofit)),]
threshold_6nt7_holo_pc3_rest_nofit = (max(contribution_6nt7_holo_rest_nofit$pc3) + min(contribution_6nt7_holo_rest_nofit$pc3))/2
contribution_6nt7_holo_rest_nofit[c(which(contribution_6nt7_holo_rest_nofit$pc3>threshold_6nt7_holo_pc3_rest_nofit)),]

ggplot(contribution_6nt7_holo_rest_nofit, aes(x = contribution_6nt7_holo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest_nofit$pc1), color = "mediumorchid", size = 0.4)  + 
  ylab("Contribution") + xlab("Residue") + ggtitle("FIRST PRINCIPAL COMPONENT") + a +
  geom_line(aes(y = threshold_6nt7_holo_pc1_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_holo_rest_nofit, aes(x = contribution_6nt7_holo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest_nofit$pc2), color = "darkorchid4", size = 0.4)  + 
  ylab("Contribution") + xlab("Residue") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc2_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_holo_rest_nofit, aes(x = contribution_6nt7_holo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest_nofit$pc3), color = "mediumpurple1", size = 0.4)   +
  ylab("Contribution") + xlab("Residue") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc3_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.1, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

#6NT7 APO
ca_ndx_6nt7_apo_rest_nofit <- atom.select(pdb_6nt7_apo_rest_nofit, elety = "CA") 
pca_6nt7_apo_rest_nofit <- pca.xyz(dcd_6nt7_apo_rest_nofit[,ca_ndx_6nt7_apo_rest_nofit$xyz])
plot(pca_6nt7_apo_rest_nofit, col=topo.colors(nrow(dcd_6nt7_apo_rest_nofit)))

pca_6nt7_apo_db_rest_nofit <- data.frame("PC" = pca_6nt7_apo_rest_nofit$z)
pca_6nt7_apo_db_rest_nofit$time <- c(1:2002)
pca_6nt7_apo_db_rest_nofit$time[1] = 100
pca_6nt7_apo_db_rest_nofit$time[(2:2002)] = rms_6nt7_apo_rest_nofit$`Time (ns)`[c(501:2501)]

A <- ggplot(pca_6nt7_apo_db_rest_nofit, aes(x=pca_6nt7_apo_db_rest_nofit$PC.1, y=pca_6nt7_apo_db_rest_nofit$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (18.08%)") + xlab("CP1 (69.23%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
B <- ggplot(pca_6nt7_apo_db_rest_nofit, aes(x=pca_6nt7_apo_db_rest_nofit$PC.3, y=pca_6nt7_apo_db_rest_nofit$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (18.08%)") + xlab("CP3 (9.17%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt7_apo_db_rest_nofit, aes(x=pca_6nt7_apo_db_rest_nofit$PC.1, y=pca_6nt7_apo_db_rest_nofit$PC.3)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP3 (9.17%)") + xlab("CP1 (69.23%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))

pc1_6nt7_apo_rest_nofit <- readPNG("pics/cerrada_sin_PC1_rest_nofitime.png")
pc2_6nt7_apo_rest_nofit  <- readPNG("pics/cerrada_sin_PC2_rest_nofitime.png")
g <- rasterGrob(pc1_6nt7_apo_rest_nofit, interpolate=TRUE)
i <- rasterGrob(pc2_6nt7_apo_rest_nofit, interpolate=TRUE)
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

eigen_bd_6nt7_apo_rest_nofit <- data.frame("eigen" = pca_6nt7_apo_rest_nofit$L)
sum_eig_6nt7_apo_rest_nofit = sum(eigen_bd_6nt7_apo_rest_nofit$eigen)
eigen_bd_6nt7_apo_rest_nofit$percent <- (eigen_bd_6nt7_apo_rest_nofit$eigen/sum_eig_6nt7_apo_rest_nofit)*100
eigen_bd_6nt7_apo_rest_nofit$num <- c(1:1182)
eigen_bd_6nt7_apo_rest_nofit$sum_percent = 0
eigen_bd_6nt7_apo_rest_nofit$sum_percent[1] = eigen_bd_6nt7_apo_rest_nofit$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_apo_rest_nofit$sum_percent[i] = eigen_bd_6nt7_apo_rest_nofit$percent[i] + eigen_bd_6nt7_apo_rest_nofit$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt7_apo_rest_nofit[c(1:7),], aes(x = eigen_bd_6nt7_apo_rest_nofit$num[c(1:7)], 
                                                 y=eigen_bd_6nt7_apo_rest_nofit$percent[c(1:7)])) + 
  geom_line(color="royalblue", size=0.8) + 
  geom_point(color="royalblue") +
  geom_text(aes(y= eigen_bd_6nt7_apo_rest_nofit$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_apo_rest_nofit$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente Principal") + scale_x_continuous(breaks = c(1:7)) + 
  ylab("Acumulated percentaje of variance") + a + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35))

mktrj.pca(pca_6nt7_apo_rest_nofit, pc=1, b=pca_6nt7_apo_rest_nofit$au[,1], file="CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_pc1_rest_nofit.pdb")
mktrj.pca(pca_6nt7_apo_rest_nofit, pc=2, b=pca_6nt7_apo_rest_nofit$au[,2], file="CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_pc2_rest_nofit.pdb")
mktrj.pca(pca_6nt7_apo_rest_nofit, pc=3, b=pca_6nt7_apo_rest_nofit$au[,3], file="CERRADA_SIN_LIGANDO/NO_FIT/6nt7_apo_pc3_rest_nofit.pdb")

contribution_6nt7_apo_rest_nofit <- data.frame("num" = pdb_6nt7_apo_rest_nofit$atom$resno[ca_ndx_6nt7_apo_rest_nofit$atom], 
                                         "res" = pdb_6nt7_apo_rest_nofit$atom$resid[ca_ndx_6nt7_apo_rest_nofit$atom],
                                         "resnum" = pdb_6nt7_apo_rest_nofit$atom$resno[ca_ndx_6nt7_apo_rest_nofit$atom],
                                         "pc1" = pca_6nt7_apo_rest_nofit$au[,1], 
                                         "pc2" = pca_6nt7_apo_rest_nofit$au[,2], 
                                         "pc3" = pca_6nt7_apo_rest_nofit$au[,3]) 

i = 146
for (j in 1:392){
  contribution_6nt7_apo_rest_nofit$num[j] = i
  i = i+1
}

threshold_6nt7_apo_pc1_rest_nofit = (max(contribution_6nt7_apo_rest_nofit$pc1) + min(contribution_6nt7_apo_rest_nofit$pc1))/2
contribution_6nt7_apo_rest_nofit[c(which(contribution_6nt7_apo_rest_nofit$pc1>threshold_6nt7_apo_pc1_rest_nofit)),]
threshold_6nt7_apo_pc2_rest_nofit = (max(contribution_6nt7_apo_rest_nofit$pc2) + min(contribution_6nt7_apo_rest_nofit$pc2))/2
contribution_6nt7_apo_rest_nofit[c(which(contribution_6nt7_apo_rest_nofit$pc2>threshold_6nt7_apo_pc2_rest_nofit)),]
threshold_6nt7_apo_pc3_rest_nofit = (max(contribution_6nt7_apo_rest_nofit$pc3) + min(contribution_6nt7_apo_rest_nofit$pc3))/2
contribution_6nt7_apo_rest_nofit[c(which(contribution_6nt7_apo_rest_nofit$pc3>threshold_6nt7_apo_pc3_rest_nofit)),]

ggplot(contribution_6nt7_apo_rest_nofit, aes(x = contribution_6nt7_apo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest_nofit$pc1), color = "dodgerblue3", size = 0.4)  + 
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_apo_pc1_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_apo_rest_nofit, aes(x = contribution_6nt7_apo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest_nofit$pc2), color = "steelblue", size = 0.4)  +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("SECOND PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc2_rest_nofit), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt7_apo_rest_nofit, aes(x = contribution_6nt7_apo_rest_nofit$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest_nofit$pc3), color = "mediumblue", size = 0.4)   +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("THIRD PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc3_rest_nofit), color = "black", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
  annotate("rect", xmin = 191, xmax=196, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=0, ymax= 0.15, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

############################################### 9- PUENTES DE HIDRÓGENO ###################################################
summary(hbonds_rest_nofit)

ggplot(hbonds_rest_nofit, aes(x = HBOND)) + geom_histogram(binwidth = 0.5, fill = "black") + a +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  ylab("Ocurrencia") + xlab("Cantidad de puentes de hidrógeno") 


ggplot(hbonds_rest_nofit, aes(x=TIME)) + 
  geom_point(aes(y = hbonds_rest_nofit$HBOND), color = "slateblue4", size = 0.6) + 
  geom_line(aes(y = hbonds_rest_nofit$HBOND), color = "black", size = 0.2, alpha = 0.5) + 
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  a + xlab("Tiempo (ns)") + ylab("Cantidad de puentes de hidrógeno") +
  geom_line(data = df_rest_nofit, aes(x= df_rest_nofit$Time, y=df_rest_nofit$avg), color = "blue", size = 0.8)

df_rest_nofit <- data.frame("Time" = c(0:100), "hbond" = 0, "count" = 0)
df_rest_nofit$Time = df_rest_nofit$Time*5
df_rest_nofit$hbond <- as.integer(df_rest_nofit$hbond)
j = 1
df_rest_nofit <- df_rest_nofit[-(1:21),]
for(i in 1: 2001){
  if(hbonds_rest_nofit$TIME[i] == df_rest_nofit$Time[j]){
    df_rest_nofit$hbond[j] = df_rest_nofit$hbond[j] + hbonds_rest_nofit$HBOND[i]
    df_rest_nofit$count[j] = df_rest_nofit$count[j] + 1
  }
  else{
    if(hbonds_rest_nofit$TIME[i] < df_rest_nofit$Time[j]){
      df_rest_nofit$hbond[j] = df_rest_nofit$hbond[j] + hbonds_rest_nofit$HBOND[i]
      df_rest_nofit$count[j] = df_rest_nofit$count[j] + 1
    }
    else{
      j = j + 1
      df_rest_nofit$hbond[j] = df_rest_nofit$hbond[j] + hbonds_rest_nofit$HBOND[i]
      df_rest_nofit$count[j] = df_rest_nofit$count[j] + 1
    }
  }
}
df_rest_nofit$avg <- df_rest_nofit$hbond/df_rest_nofit$count
ggplot(df_rest_nofit, aes(Time, avg)) + geom_line() + geom_point()
###########################################################################################################################

plots  <- list(A,B,C)
do.call(grid.arrange, plots)

grid.arrange(A,B,C,D, ncol=2, nrow=2)