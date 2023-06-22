###########################################################################################################################
###########################################################################################################################
#####                                                                                                                 #####
#####                                                SCRIPT TESINA                                                    #####
#####                                                                                                                 #####
##### 0- Archivos y cosas generales                                                                                   #####
##### 1- Analisis STING abierta                                                                                       #####
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

############################################# 0- ARCHIVOS Y COSAS GENERALES ###############################################
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
rms_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_rmsd.xvg", sep="")
names(rms_6nt6_rest)[names(rms_6nt6_rest)=="Time"] <- "Time (ns)"
names(rms_6nt6_rest)[names(rms_6nt6_rest)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_rmsf.xvg", sep="")
i=146
for (j in 1:402){
  rmsf_6nt6_rest$ATOM[j]=i
  i=i+1
}
gr_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_gr.xvg", sep="")
gr_6nt6_rest$Time = gr_6nt6_rest$Time/1000
pdb_6nt6_rest <- read.pdb("ABIERTA/6nt6_rest_fit_400ns.pdb")
dcd_6n6_rest <- read.dcd("ABIERTA/6nt6_rest_fit_400ns.dcd")  

rms_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_restringida_rmsd.xvg", sep="")
names(rms_6nt7_holo_rest)[names(rms_6nt7_holo_rest)=="Time"] <- "Time (ns)"
names(rms_6nt7_holo_rest)[names(rms_6nt7_holo_rest)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_restringida_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_holo_rest$ATOM[j]=i
  i=i+1
}
gr_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_restringida_gr.xvg", sep="")
gr_6nt7_holo_rest$Time = gr_6nt7_holo_rest$Time/1000
pdb_6nt7_holo_rest <- read.pdb("CERRADA_LIGANDO/6nt7_holo_res_400.pdb")
dcd_6nt7_holo_rest <- read.dcd("CERRADA_LIGANDO/6nt7_holo_res_400.dcd")

rms_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_rms.xvg", sep="")
names(rms_6nt7_apo_rest)[names(rms_6nt7_apo_rest)=="Time"] <- "Time (ns)"
names(rms_6nt7_apo_rest)[names(rms_6nt7_apo_rest)=="RMSD"] <- "RMSD (nm)"
rmsf_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_rmsf.xvg", sep="")
i=146
for (j in 1:394){
  rmsf_6nt7_apo_rest$ATOM[j]=i
  i=i+1
}
gr_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_gr.xvg", sep="")
gr_6nt7_apo_rest$Time = gr_6nt7_apo_rest$Time/1000
pdb_6nt7_apo_rest <- read.pdb("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_400ns.pdb")
dcd_6nt7_apo_rest <- read.dcd("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_400ns.dcd") 

pol_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_pol_rmsd.xvg", sep="")
pol_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_rms_pol.xvg", sep="")
pol_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_rest_fit_pol_rmsd.xvg", sep="")

helix_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_helix_rmsd.xvg", sep="")
helix_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_rms_helix.xvg", sep="")
helix_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_rest_fit_helix_rmsd.xvg", sep="")

loop_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_loop_rmsd.xvg", sep="")
loop_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_rms_loop.xvg", sep="")
loop_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_rest_fit_loop_rmsd.xvg", sep="")

gly_6nt6_rest <- read.csv("ABIERTA/6nt6_rest_fit_gly.xvg", sep="")
gly_6nt7_apo_rest <- read.csv("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_gly.xvg", sep="")
gly_6nt7_holo_rest <- read.csv("CERRADA_LIGANDO/6nt7_holo_rest_fit_distgly163.xvg", sep="")
gly_6nt7_holo_rest$TIME <- gly_6nt7_holo_rest$TIME/1000

bfac_6nt6_rest <- read.pdb("ABIERTA/6nt6_rest_fit_bfac.pdb")
bfac_6nt7_holo_rest <- read.pdb("CERRADA_LIGANDO/6nt7_holo_restringida_bfac.pdb")
bfac_6nt7_apo_rest <- read.pdb("CERRADA_SIN_LIGANDO/6nt7_apo_rest_fit_bfac.pdb")

hbonds_rest <- read.table("CERRADA_LIGANDO/hbonds.dat", quote="\"", comment.char="")
names(hbonds_rest)[names(hbonds_rest)=="V1"] <- "TIME"
names(hbonds_rest)[names(hbonds_rest)=="V2"] <- "HBOND"
hbonds_rest$TIME <- rms_6nt7_holo_rest$`Time (ns)`[c(501:2501)]

############################################## 1- ANALISIS STING ABIERTA ##################################################
#RMSD
ggplot(rms_6nt6_rest, aes(x=rms_6nt6_rest$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt6_rest$`RMSD (nm)`), color = "springgreen3", size = 0.8) + 
  ggtitle("RMSD OPEN STING (XYZ RESTRICTION) FITTED") + ylim(0,0.6) +
  xlab("Time (ns)") + 
  ylab("RMSD (nm)") + 
  a

#RMSF
ggplot(rmsf_6nt6_rest, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "springgreen3", size=0.8 ) + ggtitle("RMSF OPEN STING (XYZ RESTRICTION) FITTED") +
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
ggplot(gr_6nt6_rest, aes(x=Time)) + 
  geom_line(aes(y=GR), color="springgreen3", size=0.8) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (XYZ RESTRICTION) FITTED") + 
  xlab("Time (ns)") +   ylab("RG (nm)") + a + ylim(2.1, 2.4)

ggplot(gr_6nt6_rest, aes(x=Time)) + 
  geom_line(aes(y=GR), color="springgreen3", size=0.8) + 
  geom_line(aes(y=GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=GRY), color="magenta2", size=0.8) + 
  geom_line(aes(y=GRZ), color="orange", size=0.8) + 
  ggtitle("RADIUS OF GYRATION OPEN STING (XYZ RESTRICTION) FITTED") + xlab("Time (ns)") + 
  ylab("GR (nm)") + a +
  annotate("text", x=10, y=1.6, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.85, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.15, label="Y-RG", color="magenta2")


#70-105
#280-400
#400-500
############################################# 2- ANALISIS STING CERRADA CON LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_holo_rest, aes(x=rms_6nt7_holo_rest$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_holo_rest$`RMSD (nm)`), color = "darkmagenta", size=0.8) + 
  ggtitle("RMSD CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") + ylim(0,0.6)+
  xlab("Time (ns)") + ylab("RMSD (nm)") + a 

#RMSF
ggplot(rmsf_6nt7_holo_rest, aes(x=ATOM, y=nm)) + ylim(0,1) +
  geom_line(color = "darkmagenta", size=0.8) + ggtitle("RMSF CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") + 
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
ggplot(gr_6nt7_holo_rest, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo_rest$GR), color="darkmagenta", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") + ylim(2.1, 2.4) +
  xlab("Time (ns)") + ylab("RG (nm)") + a

ggplot(gr_6nt7_holo_rest, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_holo_rest$GR), color="darkmagenta", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest$GRY), color="brown", size=0.8) + 
  geom_line(aes(y=gr_6nt7_holo_rest$GRZ), color="orange", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") + 
  xlab("Time (ns)") + ylab("RG (nm)") + a +
  annotate("text", x=10, y=1.6, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.8, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.1, label="Y-RG", color="brown")
############################################# 3- ANALISIS STING CERRADA SIN LIGANDO #######################################
#RMSD
ggplot(rms_6nt7_apo_rest, aes(x=rms_6nt7_apo_rest$`Time (ns)`)) + 
  geom_line(aes(y=rms_6nt7_apo_rest$`RMSD (nm)`), color = "dodgerblue3", size=0.8) + 
  ggtitle("RMSD CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") + a + ylim(0,0.6) +
  xlab("Time (ns)") + ylab("RMSD (nm)") 

#RMSF
ggplot(rmsf_6nt7_apo_rest, aes(x=ATOM, y=nm)) + 
  geom_line(color = "dodgerblue3", size=0.8) + ggtitle("RMSD CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") + 
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
ggplot(gr_6nt7_apo_rest, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo_rest$GR), color="dodgerblue3", size=0.8) + 
  ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") + a +
  xlab("Time (ns)") + ylab("RG (nm)") + ylim(2.1, 2.4)

ggplot(gr_6nt7_apo_rest, aes(x=Time)) + 
  geom_line(aes(y=gr_6nt7_apo_rest$GR), color="dodgerblue3", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest$GRX), color="#DC143C", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest$GRY), color="magenta3", size=0.8) + 
  geom_line(aes(y=gr_6nt7_apo_rest$GRZ), color="orange", size=0.8) + 
  xlab("Time (ns)") + ylab("RG (nm)") + 
  ggtitle("RADIUS OF GYRATION CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") +
  annotate("text", x=10, y=1.65, label="X-RG", color="#DC143C") +
  annotate("text", x=10, y=1.75, label="Z-RG", color="orange") +
  annotate("text", x=22, y=2.05, label="Y-RG", color="magenta3") + a
################################################## 4- COMPARACIONES #######################################################
#RMSD
RMSD_rest <- data.frame("Time(ns)" = rms_6nt6_rest$`Time (ns)`, 
                   "6NT6" = rms_6nt6_rest$`RMSD (nm)`, 
                   "6NT7 APO" = rms_6nt7_apo_rest$`RMSD (nm)`, 
                   "6NT7 HOLO" = rms_6nt7_holo_rest$`RMSD (nm)`)
ggplot(RMSD_rest, aes(x=RMSD_rest$Time.ns.)) + 
  geom_line(aes(y=RMSD_rest$X6NT6), color = "springgreen3", size=0.8) + 
  geom_line(aes(y=RMSD_rest$X6NT7.APO), color = "dodgerblue3", size=0.8) + 
  geom_line(aes(y=RMSD_rest$X6NT7.HOLO), color = "darkmagenta", size=0.8) +
  xlab("Time (ns)") + ylab("RMSD (nm)") + ggtitle("Comparación de RMSD") + 

#RMSF
d <- c(343, 344, 345, 346, 544, 545, 546, 547)
b <- c(1:8)

for (l in 1:8){
  b[l] <- which(rmsf_6nt6_rest$ATOM == d[l])
}

apo_rest <- rmsf_6nt6_rest[c(-b),]
i=146

for (j in 1:394){
  apo_rest$ATOM[j]=i
  i=i+1
}

RMSF_rest <- data.frame("aa" = apo_rest$ATOM, 
                   "6nt6" = apo_rest$nm, 
                   "6nt7-apo" = rmsf_6nt7_apo_rest$nm, 
                   "6nt7-holo"= rmsf_6nt7_holo_rest$nm)
ggplot(RMSF_rest, aes(x=RMSF_rest$aa)) +
  geom_line(aes(y=RMSF_rest$X6nt6), color = 'springgreen3', size=0.8) +
  geom_line(aes(y=RMSF_rest$X6nt7.apo), color ='dodgerblue3', size=0.8) +
  geom_line(aes(y=RMSF_rest$X6nt7.holo), color ='darkmagenta', size=0.8) +
  ylab('RMSF (nm)') + xlab('aminoacid') + a +
  ggtitle('Comparación RSMF') +
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
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

#RADIO DE GIRO
GR_rest <- data.frame("Time (ns)"= gr_6nt6_rest$Time, 
                 "6NT6"=gr_6nt6_rest$GR, 
                 "6NT7-APO"= gr_6nt7_apo_rest$GR, 
                 "6NT7-HOLO"= gr_6nt7_holo_rest$GR)
ggplot(GR_rest, aes(x= GR_rest$Time..ns.)) + 
  geom_line(aes(y=GR_rest$X6NT6), color="springgreen3", size=0.8) + 
  geom_line(aes(y=GR_rest$X6NT7.APO), color="dodgerblue3", size=0.8) + 
  geom_line(aes(y=GR_rest$X6NT7.HOLO), color="darkmagenta", size=0.8) +
  ggtitle("Comparación Radio de Giro") +
  xlab("Time (ns)") + ylab("GR (nm)") + a

############################################### 5- DISTANCIA DE GLICINAS 163 ##############################################
ggplot(gly_6nt6_rest, aes(x=gly_6nt6_rest$TIME)) + 
  geom_line(aes(y=gly_6nt6_rest$DIS), color = "springgreen3", size=0.8, alpha = .7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a  + ylim(0.3, 1.1) +
  ggtitle("GLYCINE 163 DISTANCE OPEN STING (XYZ RESTRICTION) FITTED") +
  geom_line(data = GLY_AVG_rest, aes(x = Time, y = avg_6NT6), color="darkgreen")

ggplot(gly_6nt7_apo_rest, aes(x=gly_6nt7_apo_rest$TIME)) + 
  geom_line(aes(y=gly_6nt7_apo_rest$DIS), color = "dodgerblue3", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a + ylim(0.3, 1.1) +
    ggtitle("GLYCINE 163 DISTANCE CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") +
    geom_line(data = GLY_AVG_rest, aes(x = Time, y = avg_APO), color="blue")

ggplot(gly_6nt7_holo_rest, aes(x=gly_6nt7_holo_rest$TIME)) + 
  geom_line(aes(y=gly_6nt7_holo_rest$DIS), color = "darkmagenta", size=0.8, alpha = 0.7) + 
  xlab("Time (ns)") + ylab("Distance (nm)") + a + ylim(0.3, 1.1) +
  ggtitle("GLYCINE 163 DISTANCE CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") +
  geom_line(data = GLY_AVG_rest, aes(x = Time, y = avg_HOLO), color="slateblue3")

GLY_rest <- data.frame("TIME" = gly_6nt6_rest$TIME, 
                  "6NT6" = gly_6nt6_rest$DIS, 
                  "6NT7 APO" = gly_6nt7_apo_rest$DIS, 
                  "6NT7 HOLO" = gly_6nt7_holo_rest$DIS)

GLY_AVG_rest <- data.frame("Time" = c(0:50), "GLY_6NT6" = 0, "GLY_HOLO"=0, "GLY_APO"=0, "count_6NT6" = 0, "count_APO" = 0, "count_HOLO" = 0)
GLY_AVG_rest$Time = GLY_AVG_rest$Time*10
GLY_AVG_rest$GLY_6NT6 <- as.integer(GLY_AVG_rest$GLY_6NT6)
GLY_AVG_rest$GLY_HOLO <- as.integer(GLY_AVG_rest$GLY_HOLO)
GLY_AVG_rest$GLY_APO <- as.integer(GLY_AVG_rest$GLY_APO)
j = 1
for(i in 1: 2501){
  if(gly_6nt6_rest$TIME[i] == GLY_AVG_rest$Time[j]){
    GLY_AVG_rest$GLY_6NT6[j] = GLY_AVG_rest$GLY_6NT6[j] + gly_6nt6_rest$DIS[i]
    GLY_AVG_rest$count_6NT6[j] = GLY_AVG_rest$count_6NT6[j] + 1
  }
  else{
    if(gly_6nt6_rest$TIME[i] <GLY_AVG_rest$Time[j]){
      GLY_AVG_rest$GLY_6NT6[j] = GLY_AVG_rest$GLY_6NT6[j] + gly_6nt6_rest$DIS[i]
      GLY_AVG_rest$count_6NT6[j] = GLY_AVG_rest$count_6NT6[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest$GLY_6NT6[j] = GLY_AVG_rest$GLY_6NT6[j] + gly_6nt6_rest$DIS[i]
      GLY_AVG_rest$count_6NT6[j] = GLY_AVG_rest$count_6NT6[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_holo_rest$TIME[i] == GLY_AVG_rest$Time[j]){
    GLY_AVG_rest$GLY_HOLO[j] = GLY_AVG_rest$GLY_HOLO[j] + gly_6nt7_holo_rest$DIS[i]
    GLY_AVG_rest$count_HOLO[j] = GLY_AVG_rest$count_HOLO[j] + 1
  }
  else{
    if(gly_6nt7_holo_rest$TIME[i] <GLY_AVG_rest$Time[j]){
      GLY_AVG_rest$GLY_HOLO[j] = GLY_AVG_rest$GLY_HOLO[j] + gly_6nt7_holo_rest$DIS[i]
      GLY_AVG_rest$count_HOLO[j] = GLY_AVG_rest$count_HOLO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest$GLY_HOLO[j] = GLY_AVG_rest$GLY_HOLO[j] + gly_6nt7_holo_rest$DIS[i]
      GLY_AVG_rest$count_HOLO[j] = GLY_AVG_rest$count_HOLO[j] + 1
    }
  }
}
j = 1
for(i in 1: 2501){
  if(gly_6nt7_apo_rest$TIME[i] == GLY_AVG_rest$Time[j]){
    GLY_AVG_rest$GLY_APO[j] = GLY_AVG_rest$GLY_APO[j] + gly_6nt7_apo_rest$DIS[i]
    GLY_AVG_rest$count_APO[j] = GLY_AVG_rest$count_APO[j] + 1
  }
  else{
    if(gly_6nt7_apo_rest$TIME[i] <GLY_AVG_rest$Time[j]){
      GLY_AVG_rest$GLY_APO[j] = GLY_AVG_rest$GLY_APO[j] + gly_6nt7_apo_rest$DIS[i]
      GLY_AVG_rest$count_APO[j] = GLY_AVG_rest$count_APO[j] + 1
    }
    else{
      j = j + 1
      GLY_AVG_rest$GLY_APO[j] = GLY_AVG_rest$GLY_APO[j] + gly_6nt7_apo_rest$DIS[i]
      GLY_AVG_rest$count_APO[j] = GLY_AVG_rest$count_APO[j] + 1
    }
  }
}



GLY_AVG_rest$avg_6NT6 <- GLY_AVG_rest$GLY_6NT6/GLY_AVG_rest$count_6NT6
GLY_AVG_rest$avg_HOLO <- GLY_AVG_rest$GLY_HOLO/GLY_AVG_rest$count_HOLO
GLY_AVG_rest$avg_APO <- GLY_AVG_rest$GLY_APO/GLY_AVG_rest$count_APO


ggplot(GLY_rest, aes(x=GLY_rest$TIME)) +
  geom_line(aes(y= GLY_rest$X6NT6), color = "springgreen3", size=0.8) + 
  geom_line(aes(y= GLY_rest$X6NT7.HOLO), color = "darkmagenta", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") + a 

ggplot(GLY_rest, aes(x=GLY_rest$TIME)) +
  geom_line(aes(y= GLY_rest$X6NT6), color = "springgreen3", size=0.8, alpha = 0.6) + 
  geom_line(aes(y= GLY_rest$X6NT7.HOLO), color = "darkmagenta", size=0.8, alpha=0.6) +
  geom_line(aes(y= GLY_rest$X6NT7.APO), color = "dodgerblue3", size=0.8) +
  xlab("Tiempo (ns)") + ylab("Distancia (nm)") #+ a 

############################################### 6- RMSD DE DISTINTOS DOMINIOS #############################################
#SITIO DE POLIMERIZACION
ggplot(pol_6nt6_rest, aes(x=pol_6nt6_rest$Time)) + 
  geom_line(aes(y=pol_6nt6_rest$RMSD), color = "springgreen3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.7) + 
  ggtitle("RMSD POLYMERIZATION DOMAIN OPEN STING (XYZ RESTRICTION) FITTED")

ggplot(pol_6nt7_apo_rest, aes(x=pol_6nt7_apo_rest$Time)) + 
  geom_line(aes(y=pol_6nt7_apo_rest$RMSD), color = "dodgerblue3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.7)  + 
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED")

ggplot(pol_6nt7_holo_rest, aes(x=pol_6nt7_holo_rest$TIME)) + 
  geom_line(aes(y=pol_6nt7_holo_rest$RMSD), color = "darkmagenta", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.7) + 
  ggtitle("RMSD POLYMERIZATION DOMAIN CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED")

#HELICE CONTECTORA
ggplot(helix_6nt6_rest, aes(x=helix_6nt6_rest$Time)) + 
  geom_line(aes(y=helix_6nt6_rest$RMSD), color = "springgreen3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 1.4) + 
  ggtitle("RMSD CONECTOR HELIX OPEN STING (XYZ RESTRICTION) FITTED")

ggplot(helix_6nt7_apo_rest, aes(x=helix_6nt7_apo_rest$Time)) + 
  geom_line(aes(y=helix_6nt7_apo_rest$RMSD), color = "dodgerblue3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 1.4)  + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED")

ggplot(helix_6nt7_holo_rest, aes(x=helix_6nt7_holo_rest$TIME)) + 
  geom_line(aes(y=helix_6nt7_holo_rest$RMSD), color = "darkmagenta", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 1.4) + 
  ggtitle("RMSD CONECTOR HELIX CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED")

#LOOP CONECTOR
ggplot(loop_6nt6_rest, aes(x=loop_6nt6_rest$Time)) + 
  geom_line(aes(y=loop_6nt6_rest$RMSD), color = "springgreen3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.5) +
  ggtitle("RMSD CONECTOR LOOP OPEN STING (XYZ RESTRICTION) FITTED")

ggplot(loop_6nt7_apo_rest, aes(x=loop_6nt7_apo_rest$Time)) + 
  geom_line(aes(y=loop_6nt7_apo_rest$RMSD), color = "dodgerblue3", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a  + ylim(0, 0.5) +
  ggtitle("RMSD CONECTOR LOOP CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED")

ggplot(loop_6nt7_holo_rest, aes(x=loop_6nt7_holo_rest$TIME)) + 
  geom_line(aes(y=loop_6nt7_holo_rest$RMSD), color = "darkmagenta", size=0.8) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ylim(0, 0.5) +
  ggtitle("RMSD CONECTOR LOOP CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED")

#################################################### 7- B-FACTORS #########################################################
#6nt6
bfac_6nt6_db_rest <- data.frame("NUM" = bfac_6nt6_rest$atom$resno, 
                           "RES" = bfac_6nt6_rest$atom$resid, 
                           "BFAC" = bfac_6nt6_rest$atom$b)
prom_6nt6_rest <- mean(bfac_6nt6_db_rest$BFAC)
sd_6nt6_rest <- sd(bfac_6nt6_db_rest$BFAC)
bfac_6nt6_db_rest$NBFAC <- (bfac_6nt6_db_rest$BFAC - prom_6nt6_rest)/sd_6nt6_rest
bfac_6nt6_db_rest$RMSF <- rmsf_6nt6_rest$nm
bfac_6nt6_db_rest <- bfac_6nt6_db_rest[-c(198:201, 399:402),]

j = 146
for (i in 1:394){
  bfac_6nt6_db_rest$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt6_db_rest, aes(x = bfac_6nt6_db_rest$ATOM)) +
  geom_line(aes(y = bfac_6nt6_db_rest$NBFAC), color = "springgreen3", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-factor") + ggtitle("NORMALIZED B-FACTORS OPEN STING (XYZ RESTRICTION) FITTED") +
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
  annotate("rect", xmin = 189, xmax=200, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=-0.4, ymax= 12.7, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

    bfac_6nt6_75_rest <- which(bfac_6nt6_db_rest$NBFAC > 0.75)
bfac_6nt6_db_rest[bfac_6nt6_75_rest,]
bfac_6nt6_50_rest <- which(bfac_6nt6_db_rest$NBFAC > 0.5)
bfac_6nt6_db_rest[bfac_6nt6_50_rest,]

#6nt7 holo
bfac_6nt7_holo_db_rest <- data.frame("NUM" = bfac_6nt7_holo_rest$atom$resno, 
                                "RES" = bfac_6nt7_holo_rest$atom$resid, 
                                "BFAC" = bfac_6nt7_holo_rest$atom$b)
prom_6nt7_holo_rest <- mean(bfac_6nt7_holo_db_rest$BFAC)
sd_6nt7_holo_rest <- sd(bfac_6nt7_holo_db_rest$BFAC)
bfac_6nt7_holo_db_rest$NBFAC <- (bfac_6nt7_holo_db_rest$BFAC - prom_6nt7_holo_rest)/sd_6nt7_holo_rest
bfac_6nt7_holo_db_rest$RMSF <- rmsf_6nt7_holo_rest$nm

j = 146
for (i in 1:394){
  bfac_6nt7_holo_db_rest$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt7_holo_db_rest, aes(x = bfac_6nt7_holo_db_rest$ATOM)) +
  geom_line(aes(y = bfac_6nt7_holo_db_rest$NBFAC), color = "darkmagenta", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-factor") + ggtitle("NORMALIZED B-FACTORS CLOSE STING + LIGAND (XYZ RESTRICTION) FITTED") + 
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
  annotate("rect", xmin = 191, xmax=196, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=-0.48, ymax= 6, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

bfac_6nt7_holo_75_rest <- which(bfac_6nt7_holo_db_rest$NBFAC > 0.75)
bfac_6nt7_holo_db_rest[bfac_6nt7_holo_75_rest,]
bfac_6nt7_holo_50_rest <- which(bfac_6nt7_holo_db_rest$NBFAC > 0.5)
bfac_6nt7_holo_db_rest[bfac_6nt7_holo_50_rest,]

#6nt7 apo
bfac_6nt7_apo_db_rest <- data.frame("NUM" = bfac_6nt7_apo_rest$atom$resno, 
                               "RES" = bfac_6nt7_apo_rest$atom$resid, 
                               "BFAC" = bfac_6nt7_apo_rest$atom$o) 
bfac_6nt7_apo_db_rest$BFAC <- bfac_6nt7_apo_db_rest$BFAC + bfac_6nt7_apo_rest$atom$b
prom_6nt7_apo_rest <- mean(bfac_6nt7_apo_db_rest$BFAC)
sd_6nt7_apo_rest <- sd(bfac_6nt7_apo_db_rest$BFAC)
bfac_6nt7_apo_db_rest$NBFAC <- (bfac_6nt7_apo_db_rest$BFAC - prom_6nt7_apo_rest)/sd_6nt7_apo_rest
bfac_6nt7_apo_db_rest$RMSF <- rmsf_6nt7_apo_rest$nm

j = 146
for (i in 1:394){
  bfac_6nt7_apo_db_rest$ATOM[i] = j
  j = j+1
}

ggplot(bfac_6nt7_apo_db_rest, aes(x = bfac_6nt7_apo_db_rest$ATOM)) +
  geom_line(aes(y = bfac_6nt7_apo_db_rest$NBFAC), color = "dodgerblue3", size = 0.8) +
  a + xlab("Residue") + ylab("Normalized B-Factor") + 
  ggtitle("NORMALIZED B-FACTORS CLOSE STING NO LIGAND (XYZ RESTRICTION) FITTED") +
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
  annotate("rect", xmin = 191, xmax=196, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 388, xmax=393, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 255, xmax=258, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 452, xmax=455, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 323, xmax=327, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 520, xmax=524, ymin=-0.45, ymax= 12.7, alpha=.2, fill="grey60") +
  theme(panel.grid.major.x = element_blank(),          
        panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(color = "gray95")) + 
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

bfac_6nt7_apo_75_rest <- which(bfac_6nt7_apo_db_rest$NBFAC > 0.75)
bfac_6nt7_apo_db_rest[bfac_6nt7_apo_75_rest,]
bfac_6nt7_apo_50_rest <- which(bfac_6nt7_apo_db_rest$NBFAC > 0.5)
bfac_6nt7_apo_db_rest[bfac_6nt7_apo_50_rest,]

ggplot() +
  geom_line(data = bfac_6nt7_apo_db_rest, aes(x = bfac_6nt7_apo_db_rest$ATOM, y = bfac_6nt7_apo_db_rest$NBFAC), color = "dodgerblue3", size = 0.8) +
  geom_line(data = bfac_6nt7_holo_db_rest, aes(x= bfac_6nt7_holo_db_rest$ATOM, y= bfac_6nt7_holo_db_rest$NBFAC), color = "darkmagenta", size = 0.8) +
  geom_line(data = bfac_6nt6_db_rest, aes(x = bfac_6nt6_db_rest$ATOM, y = bfac_6nt6_db_rest$NBFAC), color = "springgreen3", size = 0.8 ) +
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
ca_ndx_6nt6_rest <- atom.select(pdb_6nt6_rest, elety = "CA") 
pca_6nt6_rest <- pca.xyz(dcd_6n6_rest[,ca_ndx_6nt6_rest$xyz])
plot(pca_6nt6_rest, col=topo.colors(nrow(dcd_6n6_rest)))

pca_6nt6_db_rest <- data.frame("PC" = pca_6nt6_rest$z)
pca_6nt6_db_rest$time <- c(1:2002)
pca_6nt6_db_rest$time[(1:2002)] = rms_6nt6_rest$`Time (ns)`[c(500:2501)]

ggplot(pca_6nt6_db_rest, aes(x=pca_6nt6_db_rest$PC.1, y=pca_6nt6_db_rest$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.77%)") + xlab("CP1 (30.99%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
ggplot(pca_6nt6_db_rest, aes(x=pca_6nt6_db_rest$PC.3, y=pca_6nt6_db_rest$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.77%)") + xlab("CP3 (8.23%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ggplot(pca_6nt6_db_rest, aes(x=pca_6nt6_db_rest$PC.1, y=pca_6nt6_db_rest$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (30.99%)") + ylab("CP3 (8.23%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

eigen_bd_6nt6_rest <- data.frame("eigen" = pca_6nt6_rest$L)
sum_eig_6nt6_rest = sum(eigen_bd_6nt6_rest$eigen)
eigen_bd_6nt6_rest$percent_6nt6 <- (eigen_bd_6nt6_rest$eigen/sum_eig_6nt6_rest)*100
eigen_bd_6nt6_rest$num <- c(1:1206)
eigen_bd_6nt6_rest$sum_percent = 0
eigen_bd_6nt6_rest$sum_percent[1] = eigen_bd_6nt6_rest$percent_6nt6[1]

for (i in 2:1206){
  eigen_bd_6nt6_rest$sum_percent[i] = eigen_bd_6nt6_rest$percent_6nt6[i] + eigen_bd_6nt6_rest$sum_percent[i - 1]
}

ggplot(eigen_bd_6nt6_rest[c(1:7),], aes(x = eigen_bd_6nt6_rest$num[c(1:7)], y=eigen_bd_6nt6_rest$percent_6nt6[c(1:7)])) + 
  geom_line(color="springgreen3", size=0.8) + 
  geom_point(color="springgreen3") + 
  geom_text(aes(y= eigen_bd_6nt6_rest$percent_6nt6[c(1:7)] + 1, 
                label = round(eigen_bd_6nt6_rest$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

pc1_6nt6_rest <- readPNG("pics/6nt6_pc1_restime.png")
g <- rasterGrob(pc1_6nt6_rest, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

grid.arrange(arrangeGrob(B, ncol = 1, nrow = 1),
             arrangeGrob(A, ncol = 1, nrow = 1), widths = c(2,3))

mktrj.pca(pca_6nt6_rest, pc=1, b=pca_6nt6_rest$au[,1], file="ABIERTA/6nt6_pc1_rest.pdb")
mktrj.pca(pca_6nt6_rest, pc=2, b=pca_6nt6_rest$au[,2], file="ABIERTA/6nt6_pc2_rest.pdb")
mktrj.pca(pca_6nt6_rest, pc=3, b=pca_6nt6_rest$au[,3], file="ABIERTA/6nt6_pc3_rest.pdb")

contribution_6nt6_rest <- data.frame("num" = pdb_6nt6_rest$atom$resno[ca_ndx_6nt6_rest$atom],
                                "resnum" = pdb_6nt6_rest$atom$resno[ca_ndx_6nt6_rest$atom],
                                "res" = pdb_6nt6_rest$atom$resid[ca_ndx_6nt6_rest$atom], 
                                "pc1" = pca_6nt6_rest$au[,1], 
                                "pc2" = pca_6nt6_rest$au[,2], 
                                "pc3" = pca_6nt6_rest$au[,3]) 
i = 146
for (j in 1:402){
  contribution_6nt6_rest$num[j] = i
  i = i+1
}

#threshold (max-min)/2
threshold_6nt6_pc1_rest = (max(contribution_6nt6_rest$pc1) + min(contribution_6nt6_rest$pc1))/2
contribution_6nt6_rest[c(which(contribution_6nt6_rest$pc1>threshold_6nt6_pc1_rest)),]
threshold_6nt6_pc2_rest = (max(contribution_6nt6_rest$pc2) + min(contribution_6nt6_rest$pc2))/2
contribution_6nt6_rest[c(which(contribution_6nt6_rest$pc2>threshold_6nt6_pc2_rest)),]
threshold_6nt6_pc3_rest = (max(contribution_6nt6_rest$pc3) + min(contribution_6nt6_rest$pc3))/2
contribution_6nt6_rest[c(which(contribution_6nt6_rest$pc3>threshold_6nt6_pc3_rest)),]

ggplot(contribution_6nt6_rest, aes(x = contribution_6nt6_rest$num)) + 
  geom_col(aes(y = contribution_6nt6_rest$pc1), color = "chartreuse3", size = 0.4)  + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc1_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)


ggplot(contribution_6nt6_rest, aes(x = contribution_6nt6_rest$num)) + 
  geom_col(aes(y = contribution_6nt6_rest$pc2), color = "darkolivegreen3", size = 0.4)  + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc2_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

ggplot(contribution_6nt6_rest, aes(x = contribution_6nt6_rest$num)) + 
  geom_col(aes(y = contribution_6nt6_rest$pc3), color = "springgreen1", size = 0.4)   + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt6_pc3_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
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
  annotate("rect", xmin = 189, xmax=200, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 390, xmax=401, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 232, xmax=244, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 433, xmax=445, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 455, xmax=459, ymin=0, ymax=0.5, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 322, xmax=328, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 523, xmax=529, ymin=0, ymax=0.5, alpha=.2, fill="grey60") +
  geom_vline(xintercept = 347, color = "black", alpha = 0.5, linetype = 5)

#6NT7 HOLO
ca_ndx_6nt7_holo_rest <- atom.select(pdb_6nt7_holo_rest, elety = "CA")
pca_6nt7_holo_rest <- pca.xyz(dcd_6nt7_holo_rest[,ca_ndx_6nt7_holo_rest$xyz])
plot(pca_6nt7_holo_rest, col=topo.colors(nrow(dcd_6nt7_holo_rest)))

pca_6nt7_holo_db_rest <- data.frame("PC" = pca_6nt7_holo_rest$z)
pca_6nt7_holo_db_rest$time <- c(1:2002)
pca_6nt7_holo_db_rest$time[(1:2002)] = rms_6nt7_holo_rest$`Time (ns)`[500:2501]

A <- ggplot(pca_6nt7_holo_db_rest, aes(x=pca_6nt7_holo_db_rest$PC.1, y=pca_6nt7_holo_db_rest$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.16%)") + xlab("CP1 (31.49%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
B <- ggplot(pca_6nt7_holo_db_rest, aes(x=pca_6nt7_holo_db_rest$PC.1, y=pca_6nt7_holo_db_rest$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP3 (8.14%)") + xlab("CP1 (31.49%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt7_holo_db_rest, aes(x=pca_6nt7_holo_db_rest$PC.3, y=pca_6nt7_holo_db_rest$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.16%)") + xlab("CP3 (8.14%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))

pc1_6nt7_holo_rest <- readPNG("pics/6nt7_holo_pc1_restime.png")
g <- rasterGrob(pc1_6nt7_holo_rest, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

pc2_6nt7_holo_rest <- readPNG("pics/6nt7_holo_pc2_restime.png")
h <- rasterGrob(pc2_6nt7_holo_rest, interpolate=TRUE)
qplot(1:10, 1:10, geom="blank", main = "b") +
  annotation_custom(h, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +  
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        axis.ticks = element_blank(), panel.background = element_rect(fill="white"),
        plot.title = element_text(colour = 'black', size =15))

grid.arrange(arrangeGrob(A, ncol = 1, nrow = 1),
             arrangeGrob(B, ncol = 1, nrow = 1), widths = c(2,3))

eigen_bd_6nt7_holo_rest <- data.frame("eigen" = pca_6nt7_holo_rest$L)
sum_eig_6nt7_holo_rest = sum(eigen_bd_6nt7_holo_rest$eigen)
eigen_bd_6nt7_holo_rest$percent <- (eigen_bd_6nt7_holo_rest$eigen/sum_eig_6nt7_holo_rest)*100
eigen_bd_6nt7_holo_rest$num <- c(1:1182)
eigen_bd_6nt7_holo_rest$sum_percent = 0
eigen_bd_6nt7_holo_rest$sum_percent[1] = eigen_bd_6nt7_holo_rest$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_holo_rest$sum_percent[i] = eigen_bd_6nt7_holo_rest$percent[i] + eigen_bd_6nt7_holo_rest$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt7_holo_rest[c(1:7),], aes(x = eigen_bd_6nt7_holo_rest$num[c(1:7)], 
                                             y=eigen_bd_6nt7_holo_rest$percent[c(1:7)])) + 
  geom_line(color="darkmagenta", size=0.8) + 
  geom_point(color="darkmagenta") + 
  geom_text(aes(y= eigen_bd_6nt7_holo_rest$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_holo_rest$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente") + 
  ylab("Acumulated percentaje of variance") + 
  a + scale_x_continuous(breaks = c(1:8)) + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))

mktrj.pca(pca_6nt7_holo_rest, pc=1, b=pca_6nt7_holo_rest$au[,1], file="CERRADA_LIGANDO/6nt7_holo_pc1_rest.pdb")
mktrj.pca(pca_6nt7_holo_rest, pc=2, b=pca_6nt7_holo_rest$au[,2], file="CERRADA_LIGANDO/6nt7_holo_pc2_rest.pdb")
mktrj.pca(pca_6nt7_holo_rest, pc=3, b=pca_6nt7_holo_rest$au[,3], file="CERRADA_LIGANDO/6nt7_holo_pc3_rest.pdb")

contribution_6nt7_holo_rest <- data.frame("num" = pdb_6nt7_holo_rest$atom$resno[ca_ndx_6nt7_holo_rest$atom],
                                     "resnum" = pdb_6nt7_holo_rest$atom$resno[ca_ndx_6nt7_holo_rest$atom],
                                     "res" = pdb_6nt7_holo_rest$atom$resid[ca_ndx_6nt7_holo_rest$atom], 
                                     "pc1" = pca_6nt7_holo_rest$au[,1], 
                                     "pc2" = pca_6nt7_holo_rest$au[,2], 
                                     "pc3" = pca_6nt7_holo_rest$au[,3]) 

i = 146
for (j in 1:392){
  contribution_6nt7_holo_rest$num[j] = i
  i = i+1
}

threshold_6nt7_holo_pc1_rest = (max(contribution_6nt7_holo_rest$pc1) + min(contribution_6nt7_holo_rest$pc1))/2
contribution_6nt7_holo_rest[c(which(contribution_6nt7_holo_rest$pc1>threshold_6nt7_holo_pc1_rest)),]
threshold_6nt7_holo_pc2_rest = (max(contribution_6nt7_holo_rest$pc2) + min(contribution_6nt7_holo_rest$pc2))/2
contribution_6nt7_holo_rest[c(which(contribution_6nt7_holo_rest$pc2>threshold_6nt7_holo_pc2_rest)),]
threshold_6nt7_holo_pc3_rest = (max(contribution_6nt7_holo_rest$pc3) + min(contribution_6nt7_holo_rest$pc3))/2
contribution_6nt7_holo_rest[c(which(contribution_6nt7_holo_rest$pc3>threshold_6nt7_holo_pc3_rest)),]

ggplot(contribution_6nt7_holo_rest, aes(x = contribution_6nt7_holo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest$pc1), color = "mediumorchid", size = 0.4)  + ylim(0, 0.5) +
  ylab("Contribution") + xlab("Residue") + ggtitle("FIRST PRINCIPAL COMPONENT") + a +
  geom_line(aes(y = threshold_6nt7_holo_pc1_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

ggplot(contribution_6nt7_holo_rest, aes(x = contribution_6nt7_holo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest$pc2), color = "darkorchid4", size = 0.4)  + ylim(0, 0.5) +
  ylab("Contribution") + xlab("Residue") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc2_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

ggplot(contribution_6nt7_holo_rest, aes(x = contribution_6nt7_holo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_holo_rest$pc3), color = "mediumpurple1", size = 0.4)   + ylim(0, 0.5) +
  ylab("Contribution") + xlab("Residue") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_holo_pc3_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

#6NT7 APO
ca_ndx_6nt7_apo_rest <- atom.select(pdb_6nt7_apo_rest, elety = "CA") 
pca_6nt7_apo_rest <- pca.xyz(dcd_6nt7_apo_rest[,ca_ndx_6nt7_apo_rest$xyz])
plot(pca_6nt7_apo_rest, col=topo.colors(nrow(dcd_6nt7_apo_rest)))

pca_6nt7_apo_db_rest <- data.frame("PC" = pca_6nt7_apo_rest$z)
pca_6nt7_apo_db_rest$time <- c(1:2002)
pca_6nt7_apo_db_rest$time[1] = 100
pca_6nt7_apo_db_rest$time[(2:2002)] = rms_6nt7_apo_rest$`Time (ns)`[c(501:2501)]

A <- ggplot(pca_6nt7_apo_db_rest, aes(x=pca_6nt7_apo_db_rest$PC.1, y=pca_6nt7_apo_db_rest$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (14.1%)") + xlab("CP1 (22.23%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
B <- ggplot(pca_6nt7_apo_db_rest, aes(x=pca_6nt7_apo_db_rest$PC.3, y=pca_6nt7_apo_db_rest$PC.2)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP2 (14.1%)") + xlab("CP3 (6.98%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
C <- ggplot(pca_6nt7_apo_db_rest, aes(x=pca_6nt7_apo_db_rest$PC.1, y=pca_6nt7_apo_db_rest$PC.3)) + 
  geom_point(aes(colour=time), size=2) + 
  ylab("CP3 (6.98%)") + xlab("CP1 (22.23%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))

pc1_6nt7_apo_rest <- readPNG("pics/cerrada_sin_PC1_restime.png")
pc2_6nt7_apo_rest  <- readPNG("pics/cerrada_sin_PC2_restime.png")
g <- rasterGrob(pc1_6nt7_apo_rest, interpolate=TRUE)
i <- rasterGrob(pc2_6nt7_apo_rest, interpolate=TRUE)
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

eigen_bd_6nt7_apo_rest <- data.frame("eigen" = pca_6nt7_apo_rest$L)
sum_eig_6nt7_apo_rest = sum(eigen_bd_6nt7_apo_rest$eigen)
eigen_bd_6nt7_apo_rest$percent <- (eigen_bd_6nt7_apo_rest$eigen/sum_eig_6nt7_apo_rest)*100
eigen_bd_6nt7_apo_rest$num <- c(1:1182)
eigen_bd_6nt7_apo_rest$sum_percent = 0
eigen_bd_6nt7_apo_rest$sum_percent[1] = eigen_bd_6nt7_apo_rest$percent[1]

for (i in 2:1182){
  eigen_bd_6nt7_apo_rest$sum_percent[i] = eigen_bd_6nt7_apo_rest$percent[i] + eigen_bd_6nt7_apo_rest$sum_percent[i - 1]
}

D <- ggplot(eigen_bd_6nt7_apo_rest[c(1:7),], aes(x = eigen_bd_6nt7_apo_rest$num[c(1:7)], 
                                            y=eigen_bd_6nt7_apo_rest$percent[c(1:7)])) + 
  geom_line(color="dodgerblue3", size=0.8) + 
  geom_point(color="dodgerblue3") +
  geom_text(aes(y= eigen_bd_6nt7_apo_rest$percent[c(1:7)] + 1, 
                label = round(eigen_bd_6nt7_apo_rest$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal Componente Principal") + scale_x_continuous(breaks = c(1:7)) + 
  ylab("Acumulated percentaje of variance") + a + scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35))

mktrj.pca(pca_6nt7_apo_rest, pc=1, b=pca_6nt7_apo_rest$au[,1], file="CERRADA_SIN_LIGANDO/6nt7_apo_pc1_rest.pdb")
mktrj.pca(pca_6nt7_apo_rest, pc=2, b=pca_6nt7_apo_rest$au[,2], file="CERRADA_SIN_LIGANDO/6nt7_apo_pc2_rest.pdb")
mktrj.pca(pca_6nt7_apo_rest, pc=3, b=pca_6nt7_apo_rest$au[,3], file="CERRADA_SIN_LIGANDO/6nt7_apo_pc3_rest.pdb")

contribution_6nt7_apo_rest <- data.frame("num" = pdb_6nt7_apo_rest$atom$resno[ca_ndx_6nt7_apo_rest$atom], 
                                    "res" = pdb_6nt7_apo_rest$atom$resid[ca_ndx_6nt7_apo_rest$atom],
                                    "resnum" = pdb_6nt7_apo_rest$atom$resno[ca_ndx_6nt7_apo_rest$atom],
                                    "pc1" = pca_6nt7_apo_rest$au[,1], 
                                    "pc2" = pca_6nt7_apo_rest$au[,2], 
                                    "pc3" = pca_6nt7_apo_rest$au[,3]) 

i = 146
for (j in 1:392){
  contribution_6nt7_apo_rest$num[j] = i
  i = i+1
}

threshold_6nt7_apo_pc1_rest = (max(contribution_6nt7_apo_rest$pc1) + min(contribution_6nt7_apo_rest$pc1))/2
contribution_6nt7_apo_rest[c(which(contribution_6nt7_apo_rest$pc1>threshold_6nt7_apo_pc1_rest)),]
threshold_6nt7_apo_pc2_rest = (max(contribution_6nt7_apo_rest$pc2) + min(contribution_6nt7_apo_rest$pc2))/2
contribution_6nt7_apo_rest[c(which(contribution_6nt7_apo_rest$pc2>threshold_6nt7_apo_pc2_rest)),]
threshold_6nt7_apo_pc3_rest = (max(contribution_6nt7_apo_rest$pc3) + min(contribution_6nt7_apo_rest$pc3))/2
contribution_6nt7_apo_rest[c(which(contribution_6nt7_apo_rest$pc3>threshold_6nt7_apo_pc3_rest)),]

ggplot(contribution_6nt7_apo_rest, aes(x = contribution_6nt7_apo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest$pc1), color = "dodgerblue3", size = 0.4)  + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = threshold_6nt7_apo_pc1_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540) + 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

ggplot(contribution_6nt7_apo_rest, aes(x = contribution_6nt7_apo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest$pc2), color = "steelblue", size = 0.4)  + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("SECOND PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc2_rest), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 540)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

ggplot(contribution_6nt7_apo_rest, aes(x = contribution_6nt7_apo_rest$num)) + 
  geom_col(aes(y = contribution_6nt7_apo_rest$pc3), color = "mediumblue", size = 0.4)   + ylim(0, 0.5) +
  xlab("Residue") + ylab("Contribution") + a + ggtitle("THIRD PRINCIPAL COMPONENT") + 
  geom_line(aes(y = threshold_6nt7_apo_pc3_rest), color = "black", alpha = 0.5, linetype = 5) + xlim(145, 549)+ 
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
  theme(axis.text.x = element_text(angle = 90, face = "bold")) +
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

############################################### 9- PUENTES DE HIDRÓGENO ###################################################
summary(hbonds_rest)

ggplot(hbonds_rest, aes(x = HBOND)) + geom_histogram(binwidth = 0.5, fill = "black") + a +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  ylab("Ocurrencia") + xlab("Cantidad de puentes de hidrógeno") 


ggplot(hbonds_rest, aes(x=TIME)) + 
  geom_point(aes(y = hbonds_rest$HBOND), color = "darkmagenta", size = 0.6) + 
  geom_line(aes(y = hbonds_rest$HBOND), color = "black", size = 0.2, alpha = 0.5) + 
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
  a + xlab("Tiempo (ns)") + ylab("Cantidad de puentes de hidrógeno") +
  geom_line(data = df_rest, aes(x= df_rest$Time, y=df_rest$avg), color = "blue", size = 0.8)

df_rest <- data.frame("Time" = c(0:100), "hbond" = 0, "count" = 0)
df_rest$Time = df_rest$Time*5
df_rest$hbond <- as.integer(df_rest$hbond)
j = 1
df_rest <- df_rest[-(1:21),]
for(i in 1: 2001){
  if(hbonds_rest$TIME[i] == df_rest$Time[j]){
    df_rest$hbond[j] = df_rest$hbond[j] + hbonds_rest$HBOND[i]
    df_rest$count[j] = df_rest$count[j] + 1
  }
  else{
    if(hbonds_rest$TIME[i] < df_rest$Time[j]){
      df_rest$hbond[j] = df_rest$hbond[j] + hbonds_rest$HBOND[i]
      df_rest$count[j] = df_rest$count[j] + 1
    }
    else{
      j = j + 1
      df_rest$hbond[j] = df_rest$hbond[j] + hbonds_rest$HBOND[i]
      df_rest$count[j] = df_rest$count[j] + 1
    }
  }
}
df_rest$avg <- df_rest$hbond/df_rest$count
ggplot(df_rest, aes(Time, avg)) + geom_line() + geom_point()
###########################################################################################################################

plots  <- list(A,B,C)
do.call(grid.arrange, plots)

grid.arrange(A,B,C,D, ncol=2, nrow=2)