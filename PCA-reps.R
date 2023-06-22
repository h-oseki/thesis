time <-seq(from = 100, to = 200.2, by = 0.2)

##### ABIERTA 1 #####
open1_pdb <- read.pdb("abierta-rep1/crop-200/pca/abierta1.pdb")
open1_dcd <- read.dcd("abierta-rep1/crop-200/pca/abierta1.dcd")

open1_ca <- atom.select(open1_pdb, elety = "CA") 
open1_pca <- pca.xyz(open1_dcd[,open1_ca$xyz])
plot(open1_pca, col=topo.colors(nrow(open1_dcd)))

open1_pca_db <- data.frame("PC" = open1_pca$z)
open1_pca_db$time <- c(1:502)
open1_pca_db$time[(1:502)] = time

open1_pc1 <- ggplot(open1_pca_db, aes(x=open1_pca_db$PC.1, y=open1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2)  + 
  ylab("CP2 (15.65%)") + xlab("CP1 (27.05%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
open1_pc2 <- ggplot(open1_pca_db, aes(x=open1_pca_db$PC.3, y=open1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.65%)") + xlab("CP3 (7.91%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
open1_pc3 <- ggplot(open1_pca_db, aes(x=open1_pca_db$PC.1, y=open1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (27.05%)") + ylab("CP3 (7.91%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

open1_eigen_bd <- data.frame("eigen" = open1_pca$L)
open1_sum_eig = sum(open1_eigen_bd$eigen)
open1_eigen_bd$percent_6nt6 <- (open1_eigen_bd$eigen/open1_sum_eig)*100
open1_eigen_bd$num <- c(1:1182)
open1_eigen_bd$sum_percent = 0
open1_eigen_bd$sum_percent[1] = open1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  open1_eigen_bd$sum_percent[i] = open1_eigen_bd$percent_6nt6[i] + open1_eigen_bd$sum_percent[i - 1]
}

open1_eigens <- ggplot(open1_eigen_bd[c(1:7),], aes(x = open1_eigen_bd$num[c(1:7)], y=open1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="palegreen", size=0.8) + 
  geom_point(color="palegreen") + 
  geom_text(aes(y= open1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(open1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("abierta-rep1/crop-200/pca/open1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(open1_pc1, open1_pc2, open1_pc3, open1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(open1_pca, pc=1, b=open1_pca$au[,1], file="abierta-rep1/crop-200/pca/open1_pca1.pdb")
mktrj.pca(open1_pca, pc=2, b=open1_pca$au[,2], file="abierta-rep1/crop-200/pca/open1_pca2.pdb")
mktrj.pca(open1_pca, pc=3, b=open1_pca$au[,3], file="abierta-rep1/crop-200/pca/open1_pca3.pdb")

open1_contribution <- data.frame("num" = open1_pdb$atom$resno[open1_ca$atom],
                                 "resnum" = open1_pdb$atom$resno[open1_ca$atom],
                                 "res" = open1_pdb$atom$resid[open1_ca$atom], 
                                 "pc1" = open1_pca$au[,1], 
                                 "pc2" = open1_pca$au[,2], 
                                 "pc3" = open1_pca$au[,3]) 
i = 146
for (j in 1:394){
  open1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
open1_threshold_pc1 = (max(open1_contribution$pc1) + min(open1_contribution$pc1))/2
open1_contribution[c(which(open1_contribution$pc1>open1_threshold_pc1)),]
open1_threshold_pc2 = (max(open1_contribution$pc2) + min(open1_contribution$pc2))/2
open1_contribution[c(which(open1_contribution$pc2>open1_threshold_pc2)),]
open1_threshold_pc3 = (max(open1_contribution$pc3) + min(open1_contribution$pc3))/2
open1_contribution[c(which(open1_contribution$pc3>open1_threshold_pc3)),]

open1_cont1 <- ggplot(open1_contribution, aes(x = open1_contribution$num)) + 
  geom_col(aes(y = open1_contribution$pc1), color = "palegreen", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


open1_cont2 <- ggplot(open1_contribution, aes(x = open1_contribution$num)) + 
  geom_col(aes(y = open1_contribution$pc2), color = "palegreen", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

open1_cont3 <- ggplot(open1_contribution, aes(x = open1_contribution$num)) + 
  geom_col(aes(y = open1_contribution$pc3), color = "palegreen", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("abierta-rep1/crop-200/pca/open1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(open1_cont1, open1_cont2, open1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(open1_ca, open1_cont1, open1_cont2, open1_cont3, open1_contribution, open1_dcd, open1_eigen_bd,
   open1_eigens, open1_pc1, open1_pc2, open1_pc3, open1_pca, open1_pca_db, open1_pdb, open1_sum_eig,
   open1_threshold_pc1, open1_threshold_pc2, open1_threshold_pc3,i, j)



##### ABIERTA 2 #####
open2_pdb <- read.pdb("abierta-rep2/pca/abierta2.pdb")
open2_dcd <- read.dcd("abierta-rep2/pca/abierta2.dcd")

open2_ca <- atom.select(open2_pdb, elety = "CA") 
open2_pca <- pca.xyz(open2_dcd[,open2_ca$xyz])
plot(open2_pca, col=topo.colors(nrow(open2_dcd)))

open2_pca_db <- data.frame("PC" = open2_pca$z)
open2_pca_db$time <- c(1:502)
open2_pca_db$time[(1:502)] = time

open2_pc1 <- ggplot(open2_pca_db, aes(x=open2_pca_db$PC.1, y=open2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.66%)") + xlab("CP1 (28.93%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
open2_pc2 <- ggplot(open2_pca_db, aes(x=open2_pca_db$PC.3, y=open2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.66%)") + xlab("CP3 (7.72%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
open2_pc3 <- ggplot(open2_pca_db, aes(x=open2_pca_db$PC.1, y=open2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (28.93%)") + ylab("CP3 (7.72%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

open2_eigen_bd <- data.frame("eigen" = open2_pca$L)
open2_sum_eig = sum(open2_eigen_bd$eigen)
open2_eigen_bd$percent_6nt6 <- (open2_eigen_bd$eigen/open2_sum_eig)*100
open2_eigen_bd$num <- c(1:1182)
open2_eigen_bd$sum_percent = 0
open2_eigen_bd$sum_percent[1] = open2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  open2_eigen_bd$sum_percent[i] = open2_eigen_bd$percent_6nt6[i] + open2_eigen_bd$sum_percent[i - 1]
}

open2_eigens <- ggplot(open2_eigen_bd[c(1:7),], aes(x = open2_eigen_bd$num[c(1:7)], y=open2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="limegreen", size=0.8) + 
  geom_point(color="limegreen") + 
  geom_text(aes(y= open2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(open2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("abierta-rep2/pca/open2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(open2_pc1, open2_pc2, open2_pc3, open2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(open2_pca, pc=1, b=open2_pca$au[,1], file="abierta-rep2/pca/open2_pca1.pdb")
mktrj.pca(open2_pca, pc=2, b=open2_pca$au[,2], file="abierta-rep2/pca/open2_pca2.pdb")
mktrj.pca(open2_pca, pc=3, b=open2_pca$au[,3], file="abierta-rep2/pca/open2_pca3.pdb")

open2_contribution <- data.frame("num" = open2_pdb$atom$resno[open2_ca$atom],
                                 "resnum" = open2_pdb$atom$resno[open2_ca$atom],
                                 "res" = open2_pdb$atom$resid[open2_ca$atom], 
                                 "pc1" = open2_pca$au[,1], 
                                 "pc2" = open2_pca$au[,2], 
                                 "pc3" = open2_pca$au[,3]) 
i = 146
for (j in 1:394){
  open2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
open2_threshold_pc1 = (max(open2_contribution$pc1) + min(open2_contribution$pc1))/2
open2_contribution[c(which(open2_contribution$pc1>open2_threshold_pc1)),]
open2_threshold_pc2 = (max(open2_contribution$pc2) + min(open2_contribution$pc2))/2
open2_contribution[c(which(open2_contribution$pc2>open2_threshold_pc2)),]
open2_threshold_pc3 = (max(open2_contribution$pc3) + min(open2_contribution$pc3))/2
open2_contribution[c(which(open2_contribution$pc3>open2_threshold_pc3)),]

open2_cont1 <- ggplot(open2_contribution, aes(x = open2_contribution$num)) + 
  geom_col(aes(y = open2_contribution$pc1), color = "limegreen", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


open2_cont2 <- ggplot(open2_contribution, aes(x = open2_contribution$num)) + 
  geom_col(aes(y = open2_contribution$pc2), color = "limegreen", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

open2_cont3 <- ggplot(open2_contribution, aes(x = open2_contribution$num)) + 
  geom_col(aes(y = open2_contribution$pc3), color = "limegreen", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = open2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("abierta-rep2/pca/open2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(open2_cont1, open2_cont2, open2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(open2_ca, open2_cont1, open2_cont2, open2_cont3, open2_contribution, open2_dcd, open2_eigen_bd,
   open2_eigens, open2_pc1, open2_pc2, open2_pc3, open2_pca, open2_pca_db, open2_pdb, open2_sum_eig,
   open2_threshold_pc1, open2_threshold_pc2, open2_threshold_pc3,i, j)
##### CERRADA 1 #####
close1_pdb <- read.pdb("cerrada-rep1/prod/pca/cerrada1.pdb")
close1_dcd <- read.dcd("cerrada-rep1/prod/pca/cerrada1.dcd")

close1_ca <- atom.select(close1_pdb, elety = "CA") 
close1_pca <- pca.xyz(close1_dcd[,close1_ca$xyz])
plot(close1_pca, col=topo.colors(nrow(close1_dcd)))

close1_pca_db <- data.frame("PC" = close1_pca$z)
close1_pca_db$time <- c(1:502)
close1_pca_db$time[(1:502)] = time

close1_pc1 <- ggplot(close1_pca_db, aes(x=close1_pca_db$PC.1, y=close1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.1%)") + xlab("CP1 (26%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
close1_pc2 <- ggplot(close1_pca_db, aes(x=close1_pca_db$PC.3, y=close1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.1%)") + xlab("CP3 (5.63%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
close1_pc3 <- ggplot(close1_pca_db, aes(x=close1_pca_db$PC.1, y=close1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (26%)") + ylab("CP3 (5.63%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

close1_eigen_bd <- data.frame("eigen" = close1_pca$L)
close1_sum_eig = sum(close1_eigen_bd$eigen)
close1_eigen_bd$percent_6nt6 <- (close1_eigen_bd$eigen/close1_sum_eig)*100
close1_eigen_bd$num <- c(1:1182)
close1_eigen_bd$sum_percent = 0
close1_eigen_bd$sum_percent[1] = close1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  close1_eigen_bd$sum_percent[i] = close1_eigen_bd$percent_6nt6[i] + close1_eigen_bd$sum_percent[i - 1]
}

close1_eigens <- ggplot(close1_eigen_bd[c(1:7),], aes(x = close1_eigen_bd$num[c(1:7)], y=close1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="blue", size=0.8) + 
  geom_point(color="blue") + 
  geom_text(aes(y= close1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(close1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-rep1/prod/pca/close1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close1_pc1, close1_pc2, close1_pc3, close1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(close1_pca, pc=1, b=close1_pca$au[,1], file="cerrada-rep1/prod/pca/close1_pca1.pdb")
mktrj.pca(close1_pca, pc=2, b=close1_pca$au[,2], file="cerrada-rep1/prod/pca/close1_pca2.pdb")
mktrj.pca(close1_pca, pc=3, b=close1_pca$au[,3], file="cerrada-rep1/prod/pca/close1_pca3.pdb")

close1_contribution <- data.frame("num" = close1_pdb$atom$resno[close1_ca$atom],
                                  "resnum" = close1_pdb$atom$resno[close1_ca$atom],
                                  "res" = close1_pdb$atom$resid[close1_ca$atom], 
                                  "pc1" = close1_pca$au[,1], 
                                  "pc2" = close1_pca$au[,2], 
                                  "pc3" = close1_pca$au[,3]) 
i = 146
for (j in 1:394){
  close1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
close1_threshold_pc1 = (max(close1_contribution$pc1) + min(close1_contribution$pc1))/2
close1_contribution[c(which(close1_contribution$pc1>close1_threshold_pc1)),]
close1_threshold_pc2 = (max(close1_contribution$pc2) + min(close1_contribution$pc2))/2
close1_contribution[c(which(close1_contribution$pc2>close1_threshold_pc2)),]
close1_threshold_pc3 = (max(close1_contribution$pc3) + min(close1_contribution$pc3))/2
close1_contribution[c(which(close1_contribution$pc3>close1_threshold_pc3)),]

close1_cont1 <- ggplot(close1_contribution, aes(x = close1_contribution$num)) + 
  geom_col(aes(y = close1_contribution$pc1), color = "blue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


close1_cont2 <- ggplot(close1_contribution, aes(x = close1_contribution$num)) + 
  geom_col(aes(y = close1_contribution$pc2), color = "blue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

close1_cont3 <- ggplot(close1_contribution, aes(x = close1_contribution$num)) + 
  geom_col(aes(y = close1_contribution$pc3), color = "blue", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-rep1/prod/pca/close1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close1_cont1, close1_cont2, close1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(close1_ca, close1_cont1, close1_cont2, close1_cont3, close1_contribution, close1_dcd, close1_eigen_bd,
   close1_eigens, close1_pc1, close1_pc2, close1_pc3, close1_pca, close1_pca_db, close1_pdb, close1_sum_eig,
   close1_threshold_pc1, close1_threshold_pc2, close1_threshold_pc3,i, j)

##### CERRADA 2 #####
close2_pdb <- read.pdb("cerrada-rep2/pca/cerrada2.pdb")
close2_dcd <- read.dcd("cerrada-rep2/pca/cerrada2.dcd")

close2_ca <- atom.select(close2_pdb, elety = "CA") 
close2_pca <- pca.xyz(close2_dcd[,close2_ca$xyz])
plot(close2_pca, col=topo.colors(nrow(close2_dcd)))

close2_pca_db <- data.frame("PC" = close2_pca$z)
close2_pca_db$time <- c(1:502)
close2_pca_db$time[(1:502)] = time

close2_pc1 <- ggplot(close2_pca_db, aes(x=close2_pca_db$PC.1, y=close2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.69%)") + xlab("CP1 (28.07%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
close2_pc2 <- ggplot(close2_pca_db, aes(x=close2_pca_db$PC.3, y=close2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.69%)") + xlab("CP3 (5.19%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
close2_pc3 <- ggplot(close2_pca_db, aes(x=close2_pca_db$PC.1, y=close2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (28.07)") + ylab("CP3 (5.19%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

close2_eigen_bd <- data.frame("eigen" = close2_pca$L)
close2_sum_eig = sum(close2_eigen_bd$eigen)
close2_eigen_bd$percent_6nt6 <- (close2_eigen_bd$eigen/close2_sum_eig)*100
close2_eigen_bd$num <- c(1:1182)
close2_eigen_bd$sum_percent = 0
close2_eigen_bd$sum_percent[1] = close2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  close2_eigen_bd$sum_percent[i] = close2_eigen_bd$percent_6nt6[i] + close2_eigen_bd$sum_percent[i - 1]
}

close2_eigens <- ggplot(close2_eigen_bd[c(1:7),], aes(x = close2_eigen_bd$num[c(1:7)], y=close2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="dodgerblue", size=0.8) + 
  geom_point(color="dodgerblue") + 
  geom_text(aes(y= close2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(close2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-rep2/pca/close2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close2_pc1, close2_pc2, close2_pc3, close2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(close2_pca, pc=1, b=close2_pca$au[,1], file="cerrada-rep2/pca/close2_pca1.pdb")
mktrj.pca(close2_pca, pc=2, b=close2_pca$au[,2], file="cerrada-rep2/pca/close2_pca2.pdb")
mktrj.pca(close2_pca, pc=3, b=close2_pca$au[,3], file="cerrada-rep2/pca/close2_pca3.pdb")

close2_contribution <- data.frame("num" = close2_pdb$atom$resno[close2_ca$atom],
                                  "resnum" = close2_pdb$atom$resno[close2_ca$atom],
                                  "res" = close2_pdb$atom$resid[close2_ca$atom], 
                                  "pc1" = close2_pca$au[,1], 
                                  "pc2" = close2_pca$au[,2], 
                                  "pc3" = close2_pca$au[,3]) 
i = 146
for (j in 1:394){
  close2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
close2_threshold_pc1 = (max(close2_contribution$pc1) + min(close2_contribution$pc1))/2
close2_contribution[c(which(close2_contribution$pc1>close2_threshold_pc1)),]
close2_threshold_pc2 = (max(close2_contribution$pc2) + min(close2_contribution$pc2))/2
close2_contribution[c(which(close2_contribution$pc2>close2_threshold_pc2)),]
close2_threshold_pc3 = (max(close2_contribution$pc3) + min(close2_contribution$pc3))/2
close2_contribution[c(which(close2_contribution$pc3>close2_threshold_pc3)),]

close2_cont1 <- ggplot(close2_contribution, aes(x = close2_contribution$num)) + 
  geom_col(aes(y = close2_contribution$pc1), color = "dodgerblue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


close2_cont2 <- ggplot(close2_contribution, aes(x = close2_contribution$num)) + 
  geom_col(aes(y = close2_contribution$pc2), color = "dodgerblue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

close2_cont3 <- ggplot(close2_contribution, aes(x = close2_contribution$num)) + 
  geom_col(aes(y = close2_contribution$pc3), color = "dodgerblue", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-rep2/pca/close2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close2_cont1, close2_cont2, close2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(close2_ca, close2_cont1, close2_cont2, close2_cont3, close2_contribution, close2_dcd, close2_eigen_bd,
   close2_eigens, close2_pc1, close2_pc2, close2_pc3, close2_pca, close2_pca_db, close2_pdb, close2_sum_eig,
   close2_threshold_pc1, close2_threshold_pc2, close2_threshold_pc3,i, j)

##### CERRADA 3 #####
close3_pdb <- read.pdb("cerrada-rep3/pca/cerrada3.pdb")
close3_dcd <- read.dcd("cerrada-rep3/pca/cerrada3.dcd")

close3_ca <- atom.select(close3_pdb, elety = "CA") 
close3_pca <- pca.xyz(close3_dcd[,close3_ca$xyz])
plot(close3_pca, col=topo.colors(nrow(close3_dcd)))

close3_pca_db <- data.frame("PC" = close3_pca$z)
close3_pca_db$time <- c(1:502)
close3_pca_db$time[(1:502)] = time

close3_pc1 <- ggplot(close3_pca_db, aes(x=close3_pca_db$PC.1, y=close3_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (16.39%)") + xlab("CP1 (34.84%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
close3_pc2 <- ggplot(close3_pca_db, aes(x=close3_pca_db$PC.3, y=close3_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (16.39%)") + xlab("CP3 (5.76%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
close3_pc3 <- ggplot(close3_pca_db, aes(x=close3_pca_db$PC.1, y=close3_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (34.84%)") + ylab("CP3 (5.76%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

close3_eigen_bd <- data.frame("eigen" = close3_pca$L)
close3_sum_eig = sum(close3_eigen_bd$eigen)
close3_eigen_bd$percent_6nt6 <- (close3_eigen_bd$eigen/close3_sum_eig)*100
close3_eigen_bd$num <- c(1:1182)
close3_eigen_bd$sum_percent = 0
close3_eigen_bd$sum_percent[1] = close3_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  close3_eigen_bd$sum_percent[i] = close3_eigen_bd$percent_6nt6[i] + close3_eigen_bd$sum_percent[i - 1]
}

close3_eigens <- ggplot(close3_eigen_bd[c(1:7),], aes(x = close3_eigen_bd$num[c(1:7)], y=close3_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="deepskyblue", size=0.8) + 
  geom_point(color="deepskyblue") + 
  geom_text(aes(y= close3_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(close3_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-rep3/pca/close3_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close3_pc1, close3_pc2, close3_pc3, close3_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(close3_pca, pc=1, b=close3_pca$au[,1], file="cerrada-rep3/pca/close3_pca1.pdb")
mktrj.pca(close3_pca, pc=2, b=close3_pca$au[,2], file="cerrada-rep3/pca/close3_pca2.pdb")
mktrj.pca(close3_pca, pc=3, b=close3_pca$au[,3], file="cerrada-rep3/pca/close3_pca3.pdb")

close3_contribution <- data.frame("num" = close3_pdb$atom$resno[close3_ca$atom],
                                  "resnum" = close3_pdb$atom$resno[close3_ca$atom],
                                  "res" = close3_pdb$atom$resid[close3_ca$atom], 
                                  "pc1" = close3_pca$au[,1], 
                                  "pc2" = close3_pca$au[,2], 
                                  "pc3" = close3_pca$au[,3]) 
i = 146
for (j in 1:394){
  close3_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
close3_threshold_pc1 = (max(close3_contribution$pc1) + min(close3_contribution$pc1))/2
close3_contribution[c(which(close3_contribution$pc1>close3_threshold_pc1)),]
close3_threshold_pc2 = (max(close3_contribution$pc2) + min(close3_contribution$pc2))/2
close3_contribution[c(which(close3_contribution$pc2>close3_threshold_pc2)),]
close3_threshold_pc3 = (max(close3_contribution$pc3) + min(close3_contribution$pc3))/2
close3_contribution[c(which(close3_contribution$pc3>close3_threshold_pc3)),]

close3_cont1 <- ggplot(close3_contribution, aes(x = close3_contribution$num)) + 
  geom_col(aes(y = close3_contribution$pc1), color = "deepskyblue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close3_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


close3_cont2 <- ggplot(close3_contribution, aes(x = close3_contribution$num)) + 
  geom_col(aes(y = close3_contribution$pc2), color = "deepskyblue", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close3_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

close3_cont3 <- ggplot(close3_contribution, aes(x = close3_contribution$num)) + 
  geom_col(aes(y = close3_contribution$pc3), color = "deepskyblue", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = close3_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-rep3/pca/close3_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(close3_cont1, close3_cont2, close3_cont3, ncol = 1, nrow = 3)
dev.off()

rm(close3_ca, close3_cont1, close3_cont2, close3_cont3, close3_contribution, close3_dcd, close3_eigen_bd,
   close3_eigens, close3_pc1, close3_pc2, close3_pc3, close3_pca, close3_pca_db, close3_pdb, close3_sum_eig,
   close3_threshold_pc1, close3_threshold_pc2, close3_threshold_pc3,i, j)
##### LIGROT 1  #####
ligrot1_pdb <- read.pdb("cerrada-ligrot-rep1/prod/pca/ligrot1.pdb")
ligrot1_dcd <- read.dcd("cerrada-ligrot-rep1/prod/pca/ligrot1.dcd")

ligrot1_ca <- atom.select(ligrot1_pdb, elety = "CA") 
ligrot1_pca <- pca.xyz(ligrot1_dcd[,ligrot1_ca$xyz])
plot(ligrot1_pca, col=topo.colors(nrow(ligrot1_dcd)))

ligrot1_pca_db <- data.frame("PC" = ligrot1_pca$z)
ligrot1_pca_db$time <- c(1:502)
ligrot1_pca_db$time[(1:502)] = time

ligrot1_pc1 <- ggplot(ligrot1_pca_db, aes(x=ligrot1_pca_db$PC.1, y=ligrot1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.78%)") + xlab("CP1 (25.52%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
ligrot1_pc2 <- ggplot(ligrot1_pca_db, aes(x=ligrot1_pca_db$PC.3, y=ligrot1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.78%)") + xlab("CP3 (5.83%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ligrot1_pc3 <- ggplot(ligrot1_pca_db, aes(x=ligrot1_pca_db$PC.1, y=ligrot1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (25.52%)") + ylab("CP3 (5.83%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

ligrot1_eigen_bd <- data.frame("eigen" = ligrot1_pca$L)
ligrot1_sum_eig = sum(ligrot1_eigen_bd$eigen)
ligrot1_eigen_bd$percent_6nt6 <- (ligrot1_eigen_bd$eigen/ligrot1_sum_eig)*100
ligrot1_eigen_bd$num <- c(1:1182)
ligrot1_eigen_bd$sum_percent = 0
ligrot1_eigen_bd$sum_percent[1] = ligrot1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  ligrot1_eigen_bd$sum_percent[i] = ligrot1_eigen_bd$percent_6nt6[i] + ligrot1_eigen_bd$sum_percent[i - 1]
}

ligrot1_eigens <- ggplot(ligrot1_eigen_bd[c(1:7),], aes(x = ligrot1_eigen_bd$num[c(1:7)], y=ligrot1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="black", size=0.8) + 
  geom_point(color="black") + 
  geom_text(aes(y= ligrot1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(ligrot1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-ligrot-rep2/pca/ligrot1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ligrot1_pc1, ligrot1_pc2, ligrot1_pc3, ligrot1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(ligrot1_pca, pc=1, b=ligrot1_pca$au[,1], file="cerrada-ligrot-rep2/pca/ligrot1_pca1.pdb")
mktrj.pca(ligrot1_pca, pc=2, b=ligrot1_pca$au[,2], file="cerrada-ligrot-rep2/pca/ligrot1_pca2.pdb")
mktrj.pca(ligrot1_pca, pc=3, b=ligrot1_pca$au[,3], file="cerrada-ligrot-rep2/pca/ligrot1_pca3.pdb")

ligrot1_contribution <- data.frame("num" = ligrot1_pdb$atom$resno[ligrot1_ca$atom],
                                   "resnum" = ligrot1_pdb$atom$resno[ligrot1_ca$atom],
                                   "res" = ligrot1_pdb$atom$resid[ligrot1_ca$atom], 
                                   "pc1" = ligrot1_pca$au[,1], 
                                   "pc2" = ligrot1_pca$au[,2], 
                                   "pc3" = ligrot1_pca$au[,3]) 
i = 146
for (j in 1:394){
  ligrot1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
ligrot1_threshold_pc1 = (max(ligrot1_contribution$pc1) + min(ligrot1_contribution$pc1))/2
ligrot1_contribution[c(which(ligrot1_contribution$pc1>ligrot1_threshold_pc1)),]
ligrot1_threshold_pc2 = (max(ligrot1_contribution$pc2) + min(ligrot1_contribution$pc2))/2
ligrot1_contribution[c(which(ligrot1_contribution$pc2>ligrot1_threshold_pc2)),]
ligrot1_threshold_pc3 = (max(ligrot1_contribution$pc3) + min(ligrot1_contribution$pc3))/2
ligrot1_contribution[c(which(ligrot1_contribution$pc3>ligrot1_threshold_pc3)),]

ligrot1_cont1 <- ggplot(ligrot1_contribution, aes(x = ligrot1_contribution$num)) + 
  geom_col(aes(y = ligrot1_contribution$pc1), color = "black", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


ligrot1_cont2 <- ggplot(ligrot1_contribution, aes(x = ligrot1_contribution$num)) + 
  geom_col(aes(y = ligrot1_contribution$pc2), color = "black", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ligrot1_cont3 <- ggplot(ligrot1_contribution, aes(x = ligrot1_contribution$num)) + 
  geom_col(aes(y = ligrot1_contribution$pc3), color = "black", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-ligrot-rep1/pca/ligrot1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ligrot1_cont1, ligrot1_cont2, ligrot1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(ligrot1_ca, ligrot1_cont1, ligrot1_cont2, ligrot1_cont3, ligrot1_contribution, ligrot1_dcd, ligrot1_eigen_bd,
   ligrot1_eigens, ligrot1_pc1, ligrot1_pc2, ligrot1_pc3, ligrot1_pca, ligrot1_pca_db, ligrot1_pdb, ligrot1_sum_eig,
   ligrot1_threshold_pc1, ligrot1_threshold_pc2, ligrot1_threshold_pc3,i, j)

##### LIGROT 2  #####
ligrot2_pdb <- read.pdb("cerrada-ligrot-rep2/pca/ligrot2.pdb")
ligrot2_dcd <- read.dcd("cerrada-ligrot-rep2/pca/ligrot2.dcd")

ligrot2_ca <- atom.select(ligrot2_pdb, elety = "CA") 
ligrot2_pca <- pca.xyz(ligrot2_dcd[,ligrot2_ca$xyz])
plot(ligrot2_pca, col=topo.colors(nrow(ligrot2_dcd)))

ligrot2_pca_db <- data.frame("PC" = ligrot2_pca$z)
ligrot2_pca_db$time <- c(1:502)
ligrot2_pca_db$time[(1:502)] = time

ligrot2_pc1 <- ggplot(ligrot2_pca_db, aes(x=ligrot2_pca_db$PC.1, y=ligrot2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.78%)") + xlab("CP1 (25.52%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
ligrot2_pc2 <- ggplot(ligrot2_pca_db, aes(x=ligrot2_pca_db$PC.3, y=ligrot2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (8.78%)") + xlab("CP3 (5.83%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
ligrot2_pc3 <- ggplot(ligrot2_pca_db, aes(x=ligrot2_pca_db$PC.1, y=ligrot2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (25.52%)") + ylab("CP3 (5.83%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

ligrot2_eigen_bd <- data.frame("eigen" = ligrot2_pca$L)
ligrot2_sum_eig = sum(ligrot2_eigen_bd$eigen)
ligrot2_eigen_bd$percent_6nt6 <- (ligrot2_eigen_bd$eigen/ligrot2_sum_eig)*100
ligrot2_eigen_bd$num <- c(1:1182)
ligrot2_eigen_bd$sum_percent = 0
ligrot2_eigen_bd$sum_percent[1] = ligrot2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  ligrot2_eigen_bd$sum_percent[i] = ligrot2_eigen_bd$percent_6nt6[i] + ligrot2_eigen_bd$sum_percent[i - 1]
}

ligrot2_eigens <- ggplot(ligrot2_eigen_bd[c(1:7),], aes(x = ligrot2_eigen_bd$num[c(1:7)], y=ligrot2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="black", size=0.8) + 
  geom_point(color="black") + 
  geom_text(aes(y= ligrot2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(ligrot2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-ligrot-rep2/pca/ligrot2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ligrot2_pc1, ligrot2_pc2, ligrot2_pc3, ligrot2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(ligrot2_pca, pc=1, b=ligrot2_pca$au[,1], file="cerrada-ligrot-rep2/pca/ligrot2_pca1.pdb")
mktrj.pca(ligrot2_pca, pc=2, b=ligrot2_pca$au[,2], file="cerrada-ligrot-rep2/pca/ligrot2_pca2.pdb")
mktrj.pca(ligrot2_pca, pc=3, b=ligrot2_pca$au[,3], file="cerrada-ligrot-rep2/pca/ligrot2_pca3.pdb")

ligrot2_contribution <- data.frame("num" = ligrot2_pdb$atom$resno[ligrot2_ca$atom],
                                   "resnum" = ligrot2_pdb$atom$resno[ligrot2_ca$atom],
                                   "res" = ligrot2_pdb$atom$resid[ligrot2_ca$atom], 
                                   "pc1" = ligrot2_pca$au[,1], 
                                   "pc2" = ligrot2_pca$au[,2], 
                                   "pc3" = ligrot2_pca$au[,3]) 
i = 146
for (j in 1:394){
  ligrot2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
ligrot2_threshold_pc1 = (max(ligrot2_contribution$pc1) + min(ligrot2_contribution$pc1))/2
ligrot2_contribution[c(which(ligrot2_contribution$pc1>ligrot2_threshold_pc1)),]
ligrot2_threshold_pc2 = (max(ligrot2_contribution$pc2) + min(ligrot2_contribution$pc2))/2
ligrot2_contribution[c(which(ligrot2_contribution$pc2>ligrot2_threshold_pc2)),]
ligrot2_threshold_pc3 = (max(ligrot2_contribution$pc3) + min(ligrot2_contribution$pc3))/2
ligrot2_contribution[c(which(ligrot2_contribution$pc3>ligrot2_threshold_pc3)),]

ligrot2_cont1 <- ggplot(ligrot2_contribution, aes(x = ligrot2_contribution$num)) + 
  geom_col(aes(y = ligrot2_contribution$pc1), color = "slategray4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


ligrot2_cont2 <- ggplot(ligrot2_contribution, aes(x = ligrot2_contribution$num)) + 
  geom_col(aes(y = ligrot2_contribution$pc2), color = "slategray4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

ligrot2_cont3 <- ggplot(ligrot2_contribution, aes(x = ligrot2_contribution$num)) + 
  geom_col(aes(y = ligrot2_contribution$pc3), color = "slategray4", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = ligrot2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-ligrot-rep2/pca/ligrot2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ligrot2_cont1, ligrot2_cont2, ligrot2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(ligrot2_ca, ligrot2_cont1, ligrot2_cont2, ligrot2_cont3, ligrot2_contribution, ligrot2_dcd, ligrot2_eigen_bd,
   ligrot2_eigens, ligrot2_pc1, ligrot2_pc2, ligrot2_pc3, ligrot2_pca, ligrot2_pca_db, ligrot2_pdb, ligrot2_sum_eig,
   ligrot2_threshold_pc1, ligrot2_threshold_pc2, ligrot2_threshold_pc3,i, j)
##### 33CGAMP 1 #####
s33cgamp1_pdb <- read.pdb("cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1.pdb")
s33cgamp1_dcd <- read.dcd("cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1.dcd")

s33cgamp1_ca <- atom.select(s33cgamp1_pdb, elety = "CA") 
s33cgamp1_pca <- pca.xyz(s33cgamp1_dcd[,s33cgamp1_ca$xyz])
plot(s33cgamp1_pca, col=topo.colors(nrow(s33cgamp1_dcd)))

s33cgamp1_pca_db <- data.frame("PC" = s33cgamp1_pca$z)
s33cgamp1_pca_db$time <- c(1:502)
s33cgamp1_pca_db$time[(1:502)] = time

s33cgamp1_pc1 <- ggplot(s33cgamp1_pca_db, aes(x=s33cgamp1_pca_db$PC.1, y=s33cgamp1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.32%)") + xlab("CP1 (25.69%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
s33cgamp1_pc2 <- ggplot(s33cgamp1_pca_db, aes(x=s33cgamp1_pca_db$PC.3, y=s33cgamp1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.32%)") + xlab("CP3 (6.66%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
s33cgamp1_pc3 <- ggplot(s33cgamp1_pca_db, aes(x=s33cgamp1_pca_db$PC.1, y=s33cgamp1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (25.69%)") + ylab("CP3 (6.66%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

s33cgamp1_eigen_bd <- data.frame("eigen" = s33cgamp1_pca$L)
s33cgamp1_sum_eig = sum(s33cgamp1_eigen_bd$eigen)
s33cgamp1_eigen_bd$percent_6nt6 <- (s33cgamp1_eigen_bd$eigen/s33cgamp1_sum_eig)*100
s33cgamp1_eigen_bd$num <- c(1:1182)
s33cgamp1_eigen_bd$sum_percent = 0
s33cgamp1_eigen_bd$sum_percent[1] = s33cgamp1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  s33cgamp1_eigen_bd$sum_percent[i] = s33cgamp1_eigen_bd$percent_6nt6[i] + s33cgamp1_eigen_bd$sum_percent[i - 1]
}

s33cgamp1_eigens <- ggplot(s33cgamp1_eigen_bd[c(1:7),], aes(x = s33cgamp1_eigen_bd$num[c(1:7)], y=s33cgamp1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="black", size=0.8) + 
  geom_point(color="black") + 
  geom_text(aes(y= s33cgamp1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(s33cgamp1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "black", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(s33cgamp1_pc1, s33cgamp1_pc2, s33cgamp1_pc3, s33cgamp1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(s33cgamp1_pca, pc=1, b=s33cgamp1_pca$au[,1], file="cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1_pca1.pdb")
mktrj.pca(s33cgamp1_pca, pc=2, b=s33cgamp1_pca$au[,2], file="cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1_pca2.pdb")
mktrj.pca(s33cgamp1_pca, pc=3, b=s33cgamp1_pca$au[,3], file="cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1_pca3.pdb")

s33cgamp1_contribution <- data.frame("num" = s33cgamp1_pdb$atom$resno[s33cgamp1_ca$atom],
                                     "resnum" = s33cgamp1_pdb$atom$resno[s33cgamp1_ca$atom],
                                     "res" = s33cgamp1_pdb$atom$resid[s33cgamp1_ca$atom], 
                                     "pc1" = s33cgamp1_pca$au[,1], 
                                     "pc2" = s33cgamp1_pca$au[,2], 
                                     "pc3" = s33cgamp1_pca$au[,3]) 
i = 146
for (j in 1:394){
  s33cgamp1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
s33cgamp1_threshold_pc1 = (max(s33cgamp1_contribution$pc1) + min(s33cgamp1_contribution$pc1))/2
s33cgamp1_contribution[c(which(s33cgamp1_contribution$pc1>s33cgamp1_threshold_pc1)),]
s33cgamp1_threshold_pc2 = (max(s33cgamp1_contribution$pc2) + min(s33cgamp1_contribution$pc2))/2
s33cgamp1_contribution[c(which(s33cgamp1_contribution$pc2>s33cgamp1_threshold_pc2)),]
s33cgamp1_threshold_pc3 = (max(s33cgamp1_contribution$pc3) + min(s33cgamp1_contribution$pc3))/2
s33cgamp1_contribution[c(which(s33cgamp1_contribution$pc3>s33cgamp1_threshold_pc3)),]

s33cgamp1_cont1 <- ggplot(s33cgamp1_contribution, aes(x = s33cgamp1_contribution$num)) + 
  geom_col(aes(y = s33cgamp1_contribution$pc1), color = "slategray4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


s33cgamp1_cont2 <- ggplot(s33cgamp1_contribution, aes(x = s33cgamp1_contribution$num)) + 
  geom_col(aes(y = s33cgamp1_contribution$pc2), color = "slategray4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

s33cgamp1_cont3 <- ggplot(s33cgamp1_contribution, aes(x = s33cgamp1_contribution$num)) + 
  geom_col(aes(y = s33cgamp1_contribution$pc3), color = "slategray4", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada+33cGAMP-rep1/prod_200/pca/s33cgamp1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(s33cgamp1_cont1, s33cgamp1_cont2, s33cgamp1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(s33cgamp1_ca, s33cgamp1_cont1, s33cgamp1_cont2, s33cgamp1_cont3, s33cgamp1_contribution, s33cgamp1_dcd, s33cgamp1_eigen_bd,
   s33cgamp1_eigens, s33cgamp1_pc1, s33cgamp1_pc2, s33cgamp1_pc3, s33cgamp1_pca, s33cgamp1_pca_db, s33cgamp1_pdb, s33cgamp1_sum_eig,
   s33cgamp1_threshold_pc1, s33cgamp1_threshold_pc2, s33cgamp1_threshold_pc3,i, j)
##### 33CGAMP 2 #####
s33cgamp2_pdb <- read.pdb("cerrada+33cGAMP-rep2/pca/s33cgamp2.pdb")
s33cgamp2_dcd <- read.dcd("cerrada+33cGAMP-rep2/pca/s33cgamp2.dcd")

s33cgamp2_ca <- atom.select(s33cgamp2_pdb, elety = "CA") 
s33cgamp2_pca <- pca.xyz(s33cgamp2_dcd[,s33cgamp2_ca$xyz])
plot(s33cgamp2_pca, col=topo.colors(nrow(s33cgamp2_dcd)))

s33cgamp2_pca_db <- data.frame("PC" = s33cgamp2_pca$z)
s33cgamp2_pca_db$time <- c(1:502)
s33cgamp2_pca_db$time[(1:502)] = time

s33cgamp2_pc1 <- ggplot(s33cgamp2_pca_db, aes(x=s33cgamp2_pca_db$PC.1, y=s33cgamp2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (13.21%)") + xlab("CP1 (14.95%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
s33cgamp2_pc2 <- ggplot(s33cgamp2_pca_db, aes(x=s33cgamp2_pca_db$PC.3, y=s33cgamp2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (13.21%)") + xlab("CP3 (7.94%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
s33cgamp2_pc3 <- ggplot(s33cgamp2_pca_db, aes(x=s33cgamp2_pca_db$PC.1, y=s33cgamp2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (14.95%)") + ylab("CP3 (7.94%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

s33cgamp2_eigen_bd <- data.frame("eigen" = s33cgamp2_pca$L)
s33cgamp2_sum_eig = sum(s33cgamp2_eigen_bd$eigen)
s33cgamp2_eigen_bd$percent_6nt6 <- (s33cgamp2_eigen_bd$eigen/s33cgamp2_sum_eig)*100
s33cgamp2_eigen_bd$num <- c(1:1182)
s33cgamp2_eigen_bd$sum_percent = 0
s33cgamp2_eigen_bd$sum_percent[1] = s33cgamp2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  s33cgamp2_eigen_bd$sum_percent[i] = s33cgamp2_eigen_bd$percent_6nt6[i] + s33cgamp2_eigen_bd$sum_percent[i - 1]
}

s33cgamp2_eigens <- ggplot(s33cgamp2_eigen_bd[c(1:7),], aes(x = s33cgamp2_eigen_bd$num[c(1:7)], y=s33cgamp2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="firebrick1", size=0.8) + 
  geom_point(color="firebrick1") + 
  geom_text(aes(y= s33cgamp2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(s33cgamp2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "firebrick1", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada+33cGAMP-rep2/pca/s33cgamp2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(s33cgamp2_pc1, s33cgamp2_pc2, s33cgamp2_pc3, s33cgamp2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(s33cgamp2_pca, pc=1, b=s33cgamp2_pca$au[,1], file="cerrada+33cGAMP-rep2/pca/s33cgamp2_pca1.pdb")
mktrj.pca(s33cgamp2_pca, pc=2, b=s33cgamp2_pca$au[,2], file="cerrada+33cGAMP-rep2/pca/s33cgamp2_pca2.pdb")
mktrj.pca(s33cgamp2_pca, pc=3, b=s33cgamp2_pca$au[,3], file="cerrada+33cGAMP-rep2/pca/s33cgamp2_pca3.pdb")

s33cgamp2_contribution <- data.frame("num" = s33cgamp2_pdb$atom$resno[s33cgamp2_ca$atom],
                                     "resnum" = s33cgamp2_pdb$atom$resno[s33cgamp2_ca$atom],
                                     "res" = s33cgamp2_pdb$atom$resid[s33cgamp2_ca$atom], 
                                     "pc1" = s33cgamp2_pca$au[,1], 
                                     "pc2" = s33cgamp2_pca$au[,2], 
                                     "pc3" = s33cgamp2_pca$au[,3]) 
i = 146
for (j in 1:394){
  s33cgamp2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
s33cgamp2_threshold_pc1 = (max(s33cgamp2_contribution$pc1) + min(s33cgamp2_contribution$pc1))/2
s33cgamp2_contribution[c(which(s33cgamp2_contribution$pc1>s33cgamp2_threshold_pc1)),]
s33cgamp2_threshold_pc2 = (max(s33cgamp2_contribution$pc2) + min(s33cgamp2_contribution$pc2))/2
s33cgamp2_contribution[c(which(s33cgamp2_contribution$pc2>s33cgamp2_threshold_pc2)),]
s33cgamp2_threshold_pc3 = (max(s33cgamp2_contribution$pc3) + min(s33cgamp2_contribution$pc3))/2
s33cgamp2_contribution[c(which(s33cgamp2_contribution$pc3>s33cgamp2_threshold_pc3)),]

s33cgamp2_cont1 <- ggplot(s33cgamp2_contribution, aes(x = s33cgamp2_contribution$num)) + 
  geom_col(aes(y = s33cgamp2_contribution$pc1), color = "firebrick1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


s33cgamp2_cont2 <- ggplot(s33cgamp2_contribution, aes(x = s33cgamp2_contribution$num)) + 
  geom_col(aes(y = s33cgamp2_contribution$pc2), color = "firebrick1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

s33cgamp2_cont3 <- ggplot(s33cgamp2_contribution, aes(x = s33cgamp2_contribution$num)) + 
  geom_col(aes(y = s33cgamp2_contribution$pc3), color = "firebrick1", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = s33cgamp2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada+33cGAMP-rep2/pca/s33cgamp2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(s33cgamp2_cont1, s33cgamp2_cont2, s33cgamp2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(s33cgamp2_ca, s33cgamp2_cont1, s33cgamp2_cont2, s33cgamp2_cont3, s33cgamp2_contribution, s33cgamp2_dcd, s33cgamp2_eigen_bd,
   s33cgamp2_eigens, s33cgamp2_pc1, s33cgamp2_pc2, s33cgamp2_pc3, s33cgamp2_pca, s33cgamp2_pca_db, s33cgamp2_pdb, s33cgamp2_sum_eig,
   s33cgamp2_threshold_pc1, s33cgamp2_threshold_pc2, s33cgamp2_threshold_pc3,i, j)
##### CDIGMP 1  #####
digmp1_pdb <- read.pdb("cerrada+cdiGMP-rep1/prod/pca/digmp1.pdb")
digmp1_dcd <- read.dcd("cerrada+cdiGMP-rep1/prod/pca/digmp1.dcd")

digmp1_ca <- atom.select(digmp1_pdb, elety = "CA") 
digmp1_pca <- pca.xyz(digmp1_dcd[,digmp1_ca$xyz])
plot(digmp1_pca, col=topo.colors(nrow(digmp1_dcd)))

digmp1_pca_db <- data.frame("PC" = digmp1_pca$z)
digmp1_pca_db$time <- c(1:502)
digmp1_pca_db$time[(1:502)] = time

digmp1_pc1 <- ggplot(digmp1_pca_db, aes(x=digmp1_pca_db$PC.1, y=digmp1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (7.56%)") + xlab("CP1 (30.35%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
digmp1_pc2 <- ggplot(digmp1_pca_db, aes(x=digmp1_pca_db$PC.3, y=digmp1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (7.56%)") + xlab("CP3 (4.82%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
digmp1_pc3 <- ggplot(digmp1_pca_db, aes(x=digmp1_pca_db$PC.1, y=digmp1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (30.35)") + ylab("CP3 (4.82%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

digmp1_eigen_bd <- data.frame("eigen" = digmp1_pca$L)
digmp1_sum_eig = sum(digmp1_eigen_bd$eigen)
digmp1_eigen_bd$percent_6nt6 <- (digmp1_eigen_bd$eigen/digmp1_sum_eig)*100
digmp1_eigen_bd$num <- c(1:1182)
digmp1_eigen_bd$sum_percent = 0
digmp1_eigen_bd$sum_percent[1] = digmp1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  digmp1_eigen_bd$sum_percent[i] = digmp1_eigen_bd$percent_6nt6[i] + digmp1_eigen_bd$sum_percent[i - 1]
}

digmp1_eigens <- ggplot(digmp1_eigen_bd[c(1:7),], aes(x = digmp1_eigen_bd$num[c(1:7)], y=digmp1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="darkorchid2", size=0.8) + 
  geom_point(color="darkorchid2") + 
  geom_text(aes(y= digmp1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(digmp1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "darkorchid2", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada+cdiGMP-rep1/prod/pca/digmp1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(digmp1_pc1, digmp1_pc2, digmp1_pc3, digmp1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(digmp1_pca, pc=1, b=digmp1_pca$au[,1], file="cerrada+cdiGMP-rep1/prod/pca/digmp1_pca1.pdb")
mktrj.pca(digmp1_pca, pc=2, b=digmp1_pca$au[,2], file="cerrada+cdiGMP-rep1/prod/pca/digmp1_pca2.pdb")
mktrj.pca(digmp1_pca, pc=3, b=digmp1_pca$au[,3], file="cerrada+cdiGMP-rep1/prod/pca/digmp1_pca3.pdb")

digmp1_contribution <- data.frame("num" = digmp1_pdb$atom$resno[digmp1_ca$atom],
                                  "resnum" = digmp1_pdb$atom$resno[digmp1_ca$atom],
                                  "res" = digmp1_pdb$atom$resid[digmp1_ca$atom], 
                                  "pc1" = digmp1_pca$au[,1], 
                                  "pc2" = digmp1_pca$au[,2], 
                                  "pc3" = digmp1_pca$au[,3]) 
i = 146
for (j in 1:394){
  digmp1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
digmp1_threshold_pc1 = (max(digmp1_contribution$pc1) + min(digmp1_contribution$pc1))/2
digmp1_contribution[c(which(digmp1_contribution$pc1>digmp1_threshold_pc1)),]
digmp1_threshold_pc2 = (max(digmp1_contribution$pc2) + min(digmp1_contribution$pc2))/2
digmp1_contribution[c(which(digmp1_contribution$pc2>digmp1_threshold_pc2)),]
digmp1_threshold_pc3 = (max(digmp1_contribution$pc3) + min(digmp1_contribution$pc3))/2
digmp1_contribution[c(which(digmp1_contribution$pc3>digmp1_threshold_pc3)),]

digmp1_cont1 <- ggplot(digmp1_contribution, aes(x = digmp1_contribution$num)) + 
  geom_col(aes(y = digmp1_contribution$pc1), color = "darkorchid2", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


digmp1_cont2 <- ggplot(digmp1_contribution, aes(x = digmp1_contribution$num)) + 
  geom_col(aes(y = digmp1_contribution$pc2), color = "darkorchid2", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

digmp1_cont3 <- ggplot(digmp1_contribution, aes(x = digmp1_contribution$num)) + 
  geom_col(aes(y = digmp1_contribution$pc3), color = "darkorchid2", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada+cdiGMP-rep1/prod/pca/digmp1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(digmp1_cont1, digmp1_cont2, digmp1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(digmp1_ca, digmp1_cont1, digmp1_cont2, digmp1_cont3, digmp1_contribution, digmp1_dcd, digmp1_eigen_bd,
   digmp1_eigens, digmp1_pc1, digmp1_pc2, digmp1_pc3, digmp1_pca, digmp1_pca_db, digmp1_pdb, digmp1_sum_eig,
   digmp1_threshold_pc1, digmp1_threshold_pc2, digmp1_threshold_pc3,i, j)
##### CDIGMP 2  #####
digmp2_pdb <- read.pdb("cerrada+cdiGMP-rep2/pca/digmp2.pdb")
digmp2_dcd <- read.dcd("cerrada+cdiGMP-rep2/pca/digmp2.dcd")

digmp2_ca <- atom.select(digmp2_pdb, elety = "CA") 
digmp2_pca <- pca.xyz(digmp2_dcd[,digmp2_ca$xyz])
plot(digmp2_pca, col=topo.colors(nrow(digmp2_dcd)))

digmp2_pca_db <- data.frame("PC" = digmp2_pca$z)
digmp2_pca_db$time <- c(1:502)
digmp2_pca_db$time[(1:502)] = time

digmp2_pc1 <- ggplot(digmp2_pca_db, aes(x=digmp2_pca_db$PC.1, y=digmp2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.83%)") + xlab("CP1 (28.14%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
digmp2_pc2 <- ggplot(digmp2_pca_db, aes(x=digmp2_pca_db$PC.3, y=digmp2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.83%)") + xlab("CP3 (5.52%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
digmp2_pc3 <- ggplot(digmp2_pca_db, aes(x=digmp2_pca_db$PC.1, y=digmp2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (28.14)") + ylab("CP3 (5.52%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

digmp2_eigen_bd <- data.frame("eigen" = digmp2_pca$L)
digmp2_sum_eig = sum(digmp2_eigen_bd$eigen)
digmp2_eigen_bd$percent_6nt6 <- (digmp2_eigen_bd$eigen/digmp2_sum_eig)*100
digmp2_eigen_bd$num <- c(1:1182)
digmp2_eigen_bd$sum_percent = 0
digmp2_eigen_bd$sum_percent[1] = digmp2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  digmp2_eigen_bd$sum_percent[i] = digmp2_eigen_bd$percent_6nt6[i] + digmp2_eigen_bd$sum_percent[i - 1]
}

digmp2_eigens <- ggplot(digmp2_eigen_bd[c(1:7),], aes(x = digmp2_eigen_bd$num[c(1:7)], y=digmp2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="magenta4", size=0.8) + 
  geom_point(color="magenta4") + 
  geom_text(aes(y= digmp2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(digmp2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "magenta4", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada+cdiGMP-rep2/pca/digmp2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(digmp2_pc1, digmp2_pc2, digmp2_pc3, digmp2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(digmp2_pca, pc=1, b=digmp2_pca$au[,1], file="cerrada+cdiGMP-rep2/pca/digmp2_pca1.pdb")
mktrj.pca(digmp2_pca, pc=2, b=digmp2_pca$au[,2], file="cerrada+cdiGMP-rep2/pca/digmp2_pca2.pdb")
mktrj.pca(digmp2_pca, pc=3, b=digmp2_pca$au[,3], file="cerrada+cdiGMP-rep2/pca/digmp2_pca3.pdb")

digmp2_contribution <- data.frame("num" = digmp2_pdb$atom$resno[digmp2_ca$atom],
                                  "resnum" = digmp2_pdb$atom$resno[digmp2_ca$atom],
                                  "res" = digmp2_pdb$atom$resid[digmp2_ca$atom], 
                                  "pc1" = digmp2_pca$au[,1], 
                                  "pc2" = digmp2_pca$au[,2], 
                                  "pc3" = digmp2_pca$au[,3]) 
i = 146
for (j in 1:394){
  digmp2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
digmp2_threshold_pc1 = (max(digmp2_contribution$pc1) + min(digmp2_contribution$pc1))/2
digmp2_contribution[c(which(digmp2_contribution$pc1>digmp2_threshold_pc1)),]
digmp2_threshold_pc2 = (max(digmp2_contribution$pc2) + min(digmp2_contribution$pc2))/2
digmp2_contribution[c(which(digmp2_contribution$pc2>digmp2_threshold_pc2)),]
digmp2_threshold_pc3 = (max(digmp2_contribution$pc3) + min(digmp2_contribution$pc3))/2
digmp2_contribution[c(which(digmp2_contribution$pc3>digmp2_threshold_pc3)),]

digmp2_cont1 <- ggplot(digmp2_contribution, aes(x = digmp2_contribution$num)) + 
  geom_col(aes(y = digmp2_contribution$pc1), color = "magenta4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


digmp2_cont2 <- ggplot(digmp2_contribution, aes(x = digmp2_contribution$num)) + 
  geom_col(aes(y = digmp2_contribution$pc2), color = "magenta4", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

digmp2_cont3 <- ggplot(digmp2_contribution, aes(x = digmp2_contribution$num)) + 
  geom_col(aes(y = digmp2_contribution$pc3), color = "magenta4", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = digmp2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada+cdiGMP-rep2/pca/digmp2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(digmp2_cont1, digmp2_cont2, digmp2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(digmp2_ca, digmp2_cont1, digmp2_cont2, digmp2_cont3, digmp2_contribution, digmp2_dcd, digmp2_eigen_bd,
   digmp2_eigens, digmp2_pc1, digmp2_pc2, digmp2_pc3, digmp2_pca, digmp2_pca_db, digmp2_pdb, digmp2_sum_eig,
   digmp2_threshold_pc1, digmp2_threshold_pc2, digmp2_threshold_pc3,i, j)
dev.off()
##### V160M  1  #####
v160m1_pdb <- read.pdb("v160m-rep1/prod_200ns/pca/v160m1.pdb")
v160m1_dcd <- read.dcd("v160m-rep1/prod_200ns/pca/v160m1.dcd")

v160m1_ca <- atom.select(v160m1_pdb, elety = "CA") 
v160m1_pca <- pca.xyz(v160m1_dcd[,v160m1_ca$xyz])
plot(v160m1_pca, col=topo.colors(nrow(v160m1_dcd)))

v160m1_pca_db <- data.frame("PC" = v160m1_pca$z)
v160m1_pca_db$time <- c(1:502)
v160m1_pca_db$time[(1:502)] = time

v160m1_pc1 <- ggplot(v160m1_pca_db, aes(x=v160m1_pca_db$PC.1, y=v160m1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.75%)") + xlab("CP1 (20.91%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
v160m1_pc2 <- ggplot(v160m1_pca_db, aes(x=v160m1_pca_db$PC.3, y=v160m1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (11.75%)") + xlab("CP3 (8.63%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
v160m1_pc3 <- ggplot(v160m1_pca_db, aes(x=v160m1_pca_db$PC.1, y=v160m1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (20.91%)") + ylab("CP3 (8.63%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

v160m1_eigen_bd <- data.frame("eigen" = v160m1_pca$L)
v160m1_sum_eig = sum(v160m1_eigen_bd$eigen)
v160m1_eigen_bd$percent_6nt6 <- (v160m1_eigen_bd$eigen/v160m1_sum_eig)*100
v160m1_eigen_bd$num <- c(1:1182)
v160m1_eigen_bd$sum_percent = 0
v160m1_eigen_bd$sum_percent[1] = v160m1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  v160m1_eigen_bd$sum_percent[i] = v160m1_eigen_bd$percent_6nt6[i] + v160m1_eigen_bd$sum_percent[i - 1]
}

v160m1_eigens <- ggplot(v160m1_eigen_bd[c(1:7),], aes(x = v160m1_eigen_bd$num[c(1:7)], y=v160m1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="deeppink", size=0.8) + 
  geom_point(color="deeppink") + 
  geom_text(aes(y= v160m1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(v160m1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "deeppink", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("v160m-rep1/prod_200ns/pca/v160m1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(v160m1_pc1, v160m1_pc2, v160m1_pc3, v160m1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(v160m1_pca, pc=1, b=v160m1_pca$au[,1], file="v160m-rep1/prod_200ns/pca/v160m1_pca1.pdb")
mktrj.pca(v160m1_pca, pc=2, b=v160m1_pca$au[,2], file="v160m-rep1/prod_200ns/pca/v160m1_pca2.pdb")
mktrj.pca(v160m1_pca, pc=3, b=v160m1_pca$au[,3], file="v160m-rep1/prod_200ns/pca/v160m1_pca3.pdb")

v160m1_contribution <- data.frame("num" = v160m1_pdb$atom$resno[v160m1_ca$atom],
                                  "resnum" = v160m1_pdb$atom$resno[v160m1_ca$atom],
                                  "res" = v160m1_pdb$atom$resid[v160m1_ca$atom], 
                                  "pc1" = v160m1_pca$au[,1], 
                                  "pc2" = v160m1_pca$au[,2], 
                                  "pc3" = v160m1_pca$au[,3]) 
i = 146
for (j in 1:394){
  v160m1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
v160m1_threshold_pc1 = (max(v160m1_contribution$pc1) + min(v160m1_contribution$pc1))/2
v160m1_contribution[c(which(v160m1_contribution$pc1>v160m1_threshold_pc1)),]
v160m1_threshold_pc2 = (max(v160m1_contribution$pc2) + min(v160m1_contribution$pc2))/2
v160m1_contribution[c(which(v160m1_contribution$pc2>v160m1_threshold_pc2)),]
v160m1_threshold_pc3 = (max(v160m1_contribution$pc3) + min(v160m1_contribution$pc3))/2
v160m1_contribution[c(which(v160m1_contribution$pc3>v160m1_threshold_pc3)),]

v160m1_cont1 <- ggplot(v160m1_contribution, aes(x = v160m1_contribution$num)) + 
  geom_col(aes(y = v160m1_contribution$pc1), color = "deeppink", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


v160m1_cont2 <- ggplot(v160m1_contribution, aes(x = v160m1_contribution$num)) + 
  geom_col(aes(y = v160m1_contribution$pc2), color = "deeppink", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

v160m1_cont3 <- ggplot(v160m1_contribution, aes(x = v160m1_contribution$num)) + 
  geom_col(aes(y = v160m1_contribution$pc3), color = "deeppink", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("v160m-rep1/prod_200ns/pca/v160m1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(v160m1_cont1, v160m1_cont2, v160m1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(v160m1_ca, v160m1_cont1, v160m1_cont2, v160m1_cont3, v160m1_contribution, v160m1_dcd, v160m1_eigen_bd,
   v160m1_eigens, v160m1_pc1, v160m1_pc2, v160m1_pc3, v160m1_pca, v160m1_pca_db, v160m1_pdb, v160m1_sum_eig,
   v160m1_threshold_pc1, v160m1_threshold_pc2, v160m1_threshold_pc3,i, j)
dev.off()
##### V160M  2  #####
v160m2_pdb <- read.pdb("v160m-rep2/pca/v160m2.pdb")
v160m2_dcd <- read.dcd("v160m-rep2/pca/v160m2.dcd")

v160m2_ca <- atom.select(v160m2_pdb, elety = "CA") 
v160m2_pca <- pca.xyz(v160m2_dcd[,v160m2_ca$xyz])
plot(v160m2_pca, col=topo.colors(nrow(v160m2_dcd)))

v160m2_pca_db <- data.frame("PC" = v160m2_pca$z)
v160m2_pca_db$time <- c(1:502)
v160m2_pca_db$time[(1:502)] = time

v160m2_pc1 <- ggplot(v160m2_pca_db, aes(x=v160m2_pca_db$PC.1, y=v160m2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.71%)") + xlab("CP1 (28.9%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
v160m2_pc2 <- ggplot(v160m2_pca_db, aes(x=v160m2_pca_db$PC.3, y=v160m2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (15.71%)") + xlab("CP3 (6.17%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
v160m2_pc3 <- ggplot(v160m2_pca_db, aes(x=v160m2_pca_db$PC.1, y=v160m2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (28.9%)") + ylab("CP3 (6.17%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

v160m2_eigen_bd <- data.frame("eigen" = v160m2_pca$L)
v160m2_sum_eig = sum(v160m2_eigen_bd$eigen)
v160m2_eigen_bd$percent_6nt6 <- (v160m2_eigen_bd$eigen/v160m2_sum_eig)*100
v160m2_eigen_bd$num <- c(1:1182)
v160m2_eigen_bd$sum_percent = 0
v160m2_eigen_bd$sum_percent[1] = v160m2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  v160m2_eigen_bd$sum_percent[i] = v160m2_eigen_bd$percent_6nt6[i] + v160m2_eigen_bd$sum_percent[i - 1]
}

v160m2_eigens <- ggplot(v160m2_eigen_bd[c(1:7),], aes(x = v160m2_eigen_bd$num[c(1:7)], y=v160m2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="hotpink1", size=0.8) + 
  geom_point(color="hotpink1") + 
  geom_text(aes(y= v160m2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(v160m2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "hotpink1", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("v160m-rep2/pca/v160m2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(v160m2_pc1, v160m2_pc2, v160m2_pc3, v160m2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(v160m2_pca, pc=1, b=v160m2_pca$au[,1], file="v160m-rep2/pca/v160m2_pca1.pdb")
mktrj.pca(v160m2_pca, pc=2, b=v160m2_pca$au[,2], file="v160m-rep2/pca/v160m2_pca2.pdb")
mktrj.pca(v160m2_pca, pc=3, b=v160m2_pca$au[,3], file="v160m-rep2/pca/v160m2_pca3.pdb")

v160m2_contribution <- data.frame("num" = v160m2_pdb$atom$resno[v160m2_ca$atom],
                                  "resnum" = v160m2_pdb$atom$resno[v160m2_ca$atom],
                                  "res" = v160m2_pdb$atom$resid[v160m2_ca$atom], 
                                  "pc1" = v160m2_pca$au[,1], 
                                  "pc2" = v160m2_pca$au[,2], 
                                  "pc3" = v160m2_pca$au[,3]) 
i = 146
for (j in 1:394){
  v160m2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
v160m2_threshold_pc1 = (max(v160m2_contribution$pc1) + min(v160m2_contribution$pc1))/2
v160m2_contribution[c(which(v160m2_contribution$pc1>v160m2_threshold_pc1)),]
v160m2_threshold_pc2 = (max(v160m2_contribution$pc2) + min(v160m2_contribution$pc2))/2
v160m2_contribution[c(which(v160m2_contribution$pc2>v160m2_threshold_pc2)),]
v160m2_threshold_pc3 = (max(v160m2_contribution$pc3) + min(v160m2_contribution$pc3))/2
v160m2_contribution[c(which(v160m2_contribution$pc3>v160m2_threshold_pc3)),]

v160m2_cont1 <- ggplot(v160m2_contribution, aes(x = v160m2_contribution$num)) + 
  geom_col(aes(y = v160m2_contribution$pc1), color = "hotpink1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


v160m2_cont2 <- ggplot(v160m2_contribution, aes(x = v160m2_contribution$num)) + 
  geom_col(aes(y = v160m2_contribution$pc2), color = "hotpink1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

v160m2_cont3 <- ggplot(v160m2_contribution, aes(x = v160m2_contribution$num)) + 
  geom_col(aes(y = v160m2_contribution$pc3), color = "hotpink1", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = v160m2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("v160m-rep2/pca/v160m2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(v160m2_cont1, v160m2_cont2, v160m2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(v160m2_ca, v160m2_cont1, v160m2_cont2, v160m2_cont3, v160m2_contribution, v160m2_dcd, v160m2_eigen_bd,
   v160m2_eigens, v160m2_pc1, v160m2_pc2, v160m2_pc3, v160m2_pca, v160m2_pca_db, v160m2_pdb, v160m2_sum_eig,
   v160m2_threshold_pc1, v160m2_threshold_pc2, v160m2_threshold_pc3,i, j)
dev.off()
##### SINLIG  1 #####
sinlig1_pdb <- read.pdb("cerrada-sinlig-rep1/prod_200/pca/sinlig1.pdb")
sinlig1_dcd <- read.dcd("cerrada-sinlig-rep1/prod_200/pca/sinlig1.dcd")

sinlig1_ca <- atom.select(sinlig1_pdb, elety = "CA") 
sinlig1_pca <- pca.xyz(sinlig1_dcd[,sinlig1_ca$xyz])
plot(sinlig1_pca, col=topo.colors(nrow(sinlig1_dcd)))

sinlig1_pca_db <- data.frame("PC" = sinlig1_pca$z)
sinlig1_pca_db$time <- c(1:502)
sinlig1_pca_db$time[(1:502)] = time

sinlig1_pc1 <- ggplot(sinlig1_pca_db, aes(x=sinlig1_pca_db$PC.1, y=sinlig1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (12.09%)") + xlab("CP1 (23.27%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
sinlig1_pc2 <- ggplot(sinlig1_pca_db, aes(x=sinlig1_pca_db$PC.3, y=sinlig1_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (12.09%)") + xlab("CP3 (6.85%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
sinlig1_pc3 <- ggplot(sinlig1_pca_db, aes(x=sinlig1_pca_db$PC.1, y=sinlig1_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (23.27%)") + ylab("CP3 (6.85%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

sinlig1_eigen_bd <- data.frame("eigen" = sinlig1_pca$L)
sinlig1_sum_eig = sum(sinlig1_eigen_bd$eigen)
sinlig1_eigen_bd$percent_6nt6 <- (sinlig1_eigen_bd$eigen/sinlig1_sum_eig)*100
sinlig1_eigen_bd$num <- c(1:1182)
sinlig1_eigen_bd$sum_percent = 0
sinlig1_eigen_bd$sum_percent[1] = sinlig1_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  sinlig1_eigen_bd$sum_percent[i] = sinlig1_eigen_bd$percent_6nt6[i] + sinlig1_eigen_bd$sum_percent[i - 1]
}

sinlig1_eigens <- ggplot(sinlig1_eigen_bd[c(1:7),], aes(x = sinlig1_eigen_bd$num[c(1:7)], y=sinlig1_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="chocolate1", size=0.8) + 
  geom_point(color="chocolate1") + 
  geom_text(aes(y= sinlig1_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(sinlig1_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "chocolate1", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-sinlig-rep1/prod_200/pca/sinlig1_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(sinlig1_pc1, sinlig1_pc2, sinlig1_pc3, sinlig1_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(sinlig1_pca, pc=1, b=sinlig1_pca$au[,1], file="cerrada-sinlig-rep1/prod_200/pca/sinlig1_pca1.pdb")
mktrj.pca(sinlig1_pca, pc=2, b=sinlig1_pca$au[,2], file="cerrada-sinlig-rep1/prod_200/pca/sinlig1_pca2.pdb")
mktrj.pca(sinlig1_pca, pc=3, b=sinlig1_pca$au[,3], file="cerrada-sinlig-rep1/prod_200/pca/sinlig1_pca3.pdb")

sinlig1_contribution <- data.frame("num" = sinlig1_pdb$atom$resno[sinlig1_ca$atom],
                                   "resnum" = sinlig1_pdb$atom$resno[sinlig1_ca$atom],
                                   "res" = sinlig1_pdb$atom$resid[sinlig1_ca$atom], 
                                   "pc1" = sinlig1_pca$au[,1], 
                                   "pc2" = sinlig1_pca$au[,2], 
                                   "pc3" = sinlig1_pca$au[,3]) 
i = 146
for (j in 1:394){
  sinlig1_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
sinlig1_threshold_pc1 = (max(sinlig1_contribution$pc1) + min(sinlig1_contribution$pc1))/2
sinlig1_contribution[c(which(sinlig1_contribution$pc1>sinlig1_threshold_pc1)),]
sinlig1_threshold_pc2 = (max(sinlig1_contribution$pc2) + min(sinlig1_contribution$pc2))/2
sinlig1_contribution[c(which(sinlig1_contribution$pc2>sinlig1_threshold_pc2)),]
sinlig1_threshold_pc3 = (max(sinlig1_contribution$pc3) + min(sinlig1_contribution$pc3))/2
sinlig1_contribution[c(which(sinlig1_contribution$pc3>sinlig1_threshold_pc3)),]

sinlig1_cont1 <- ggplot(sinlig1_contribution, aes(x = sinlig1_contribution$num)) + 
  geom_col(aes(y = sinlig1_contribution$pc1), color = "chocolate1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig1_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


sinlig1_cont2 <- ggplot(sinlig1_contribution, aes(x = sinlig1_contribution$num)) + 
  geom_col(aes(y = sinlig1_contribution$pc2), color = "chocolate1", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig1_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

sinlig1_cont3 <- ggplot(sinlig1_contribution, aes(x = sinlig1_contribution$num)) + 
  geom_col(aes(y = sinlig1_contribution$pc3), color = "chocolate1", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig1_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-sinlig-rep1/prod_200/pca/sinlig1_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(sinlig1_cont1, sinlig1_cont2, sinlig1_cont3, ncol = 1, nrow = 3)
dev.off()

rm(sinlig1_ca, sinlig1_cont1, sinlig1_cont2, sinlig1_cont3, sinlig1_contribution, sinlig1_dcd, sinlig1_eigen_bd,
   sinlig1_eigens, sinlig1_pc1, sinlig1_pc2, sinlig1_pc3, sinlig1_pca, sinlig1_pca_db, sinlig1_pdb, sinlig1_sum_eig,
   sinlig1_threshold_pc1, sinlig1_threshold_pc2, sinlig1_threshold_pc3,i, j)
dev.off()
##### SINLIG  2 #####
sinlig2_pdb <- read.pdb("cerrada-sinlig-rep2/pca/sinlig2.pdb")
sinlig2_dcd <- read.dcd("cerrada-sinlig-rep2/pca/sinlig2.dcd")

sinlig2_ca <- atom.select(sinlig2_pdb, elety = "CA") 
sinlig2_pca <- pca.xyz(sinlig2_dcd[,sinlig2_ca$xyz])
plot(sinlig2_pca, col=topo.colors(nrow(sinlig2_dcd)))

sinlig2_pca_db <- data.frame("PC" = sinlig2_pca$z)
sinlig2_pca_db$time <- c(1:502)
sinlig2_pca_db$time[(1:502)] = time

sinlig2_pc1 <- ggplot(sinlig2_pca_db, aes(x=sinlig2_pca_db$PC.1, y=sinlig2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.67%)") + xlab("CP1 (15.57%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 
sinlig2_pc2 <- ggplot(sinlig2_pca_db, aes(x=sinlig2_pca_db$PC.3, y=sinlig2_pca_db$PC.2, colour=time)) + 
  geom_point(size = 2) + 
  ylab("CP2 (9.67%)") + xlab("CP3 (6.73%)") + a +
  scale_color_gradientn(colours = c("orange", "red", "magenta"))
sinlig2_pc3 <- ggplot(sinlig2_pca_db, aes(x=sinlig2_pca_db$PC.1, y=sinlig2_pca_db$PC.3, colour=time)) + 
  geom_point(size = 2) + 
  xlab("CP1 (15.57%)") + ylab("CP3 (6.73%)") + a + 
  scale_color_gradientn(colours = c("orange", "red", "magenta")) 

sinlig2_eigen_bd <- data.frame("eigen" = sinlig2_pca$L)
sinlig2_sum_eig = sum(sinlig2_eigen_bd$eigen)
sinlig2_eigen_bd$percent_6nt6 <- (sinlig2_eigen_bd$eigen/sinlig2_sum_eig)*100
sinlig2_eigen_bd$num <- c(1:1182)
sinlig2_eigen_bd$sum_percent = 0
sinlig2_eigen_bd$sum_percent[1] = sinlig2_eigen_bd$percent_6nt6[1]

for (i in 2:1182){
  sinlig2_eigen_bd$sum_percent[i] = sinlig2_eigen_bd$percent_6nt6[i] + sinlig2_eigen_bd$sum_percent[i - 1]
}

sinlig2_eigens <- ggplot(sinlig2_eigen_bd[c(1:7),], aes(x = sinlig2_eigen_bd$num[c(1:7)], y=sinlig2_eigen_bd$percent_6nt6[c(1:7)])) + 
  geom_line(color="gold", size=0.8) + 
  geom_point(color="gold") + 
  geom_text(aes(y= sinlig2_eigen_bd$percent_6nt6[c(1:7)] + 1, 
                label = round(sinlig2_eigen_bd$sum_percent[c(1:7)], 2)), 
            color = "gold", size = 4, hjust=0) +
  xlab("Principal component") +
  ylab("Acumulated percentaje of variance") + a + scale_x_continuous(breaks = c(1:8)) + 
  scale_y_continuous(breaks = c(0,5,10,15,20,25,30,35,40,45))



png("cerrada-sinlig-rep2/pca/sinlig2_pca.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(sinlig2_pc1, sinlig2_pc2, sinlig2_pc3, sinlig2_eigens,
          ncol = 2, nrow = 2)
dev.off()


mktrj.pca(sinlig2_pca, pc=1, b=sinlig2_pca$au[,1], file="cerrada-sinlig-rep2/pca/sinlig2_pca1.pdb")
mktrj.pca(sinlig2_pca, pc=2, b=sinlig2_pca$au[,2], file="cerrada-sinlig-rep2/pca/sinlig2_pca2.pdb")
mktrj.pca(sinlig2_pca, pc=3, b=sinlig2_pca$au[,3], file="cerrada-sinlig-rep2/pca/sinlig2_pca3.pdb")

sinlig2_contribution <- data.frame("num" = sinlig2_pdb$atom$resno[sinlig2_ca$atom],
                                   "resnum" = sinlig2_pdb$atom$resno[sinlig2_ca$atom],
                                   "res" = sinlig2_pdb$atom$resid[sinlig2_ca$atom], 
                                   "pc1" = sinlig2_pca$au[,1], 
                                   "pc2" = sinlig2_pca$au[,2], 
                                   "pc3" = sinlig2_pca$au[,3]) 
i = 146
for (j in 1:394){
  sinlig2_contribution$num[j] = i
  i = i+1
}

#threshold (max-min)/2
sinlig2_threshold_pc1 = (max(sinlig2_contribution$pc1) + min(sinlig2_contribution$pc1))/2
sinlig2_contribution[c(which(sinlig2_contribution$pc1>sinlig2_threshold_pc1)),]
sinlig2_threshold_pc2 = (max(sinlig2_contribution$pc2) + min(sinlig2_contribution$pc2))/2
sinlig2_contribution[c(which(sinlig2_contribution$pc2>sinlig2_threshold_pc2)),]
sinlig2_threshold_pc3 = (max(sinlig2_contribution$pc3) + min(sinlig2_contribution$pc3))/2
sinlig2_contribution[c(which(sinlig2_contribution$pc3>sinlig2_threshold_pc3)),]

sinlig2_cont1 <- ggplot(sinlig2_contribution, aes(x = sinlig2_contribution$num)) + 
  geom_col(aes(y = sinlig2_contribution$pc1), color = "gold", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("FIRST PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig2_threshold_pc1), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)


sinlig2_cont2 <- ggplot(sinlig2_contribution, aes(x = sinlig2_contribution$num)) + 
  geom_col(aes(y = sinlig2_contribution$pc2), color = "gold", size = 0.4)  + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("SECOND PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig2_threshold_pc2), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

sinlig2_cont3 <- ggplot(sinlig2_contribution, aes(x = sinlig2_contribution$num)) + 
  geom_col(aes(y = sinlig2_contribution$pc3), color = "gold", size = 0.4)   + ylim(0, 0.4) +
  xlab("Residue") + ylab("Contribution") + ggtitle("THIRD PRINCIPAL COMPONENT") + a + 
  geom_line(aes(y = sinlig2_threshold_pc3), color = "red", alpha = 0.5, linetype = 5) + xlim(145, 549) +
  scale_x_continuous(breaks = c(146, 166, 186, 206, 226, 246, 266, 286, 306, 326,  
                                346, 366, 386, 406, 426, 446, 466, 486, 506, 526),
                     labels = paste(c("ALA146A", "TRP166A", "GLU186A", "LEU206A", "TYR226A",  
                                      "LYS246A", "PHE266A", "ARG286A", "SER306A", "GLU326A",  
                                      "VAL149B", "TYR169B", "ARG189B", "LEU209B", "ASP229B",
                                      "LEU249B", "PRO269B", "ARG289B", "GLU309B", "PHE329B"))) +
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
  geom_vline(xintercept = 343, color = "black", alpha = 0.5, linetype = 5)

png("cerrada-sinlig-rep2/pca/sinlig2_pca_cont.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(sinlig2_cont1, sinlig2_cont2, sinlig2_cont3, ncol = 1, nrow = 3)
dev.off()

rm(sinlig2_ca, sinlig2_cont1, sinlig2_cont2, sinlig2_cont3, sinlig2_contribution, sinlig2_dcd, sinlig2_eigen_bd,
   sinlig2_eigens, sinlig2_pc1, sinlig2_pc2, sinlig2_pc3, sinlig2_pca, sinlig2_pca_db, sinlig2_pdb, sinlig2_sum_eig,
   sinlig2_threshold_pc1, sinlig2_threshold_pc2, sinlig2_threshold_pc3,i, j)
dev.off()