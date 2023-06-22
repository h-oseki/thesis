#######################################################################################
######                                                                           ######
######                               FIGURAS PAPER                               ######
######                                                                           ######
#######################################################################################

#######################################################################################
setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 12),
  axis.title = element_text(colour = 'black', face = 'italic', size = 10),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)
##
##################################### Figura 1 ########################################
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

bfacs_f1 <- data.frame("res" = c(open2_db$ATOM, close2_db$ATOM),
                    "nbfac" = c(open2_db$NBFAC,close2_db$NBFAC),
                    "Structure" = "label")

bfacs_f1$Structure[1:394] <- "Open STING" #palegreen
bfacs_f1$Structure[395:788] <- "Closed STING" #limegreen

bfacs_f1$Structure <- as.factor(bfacs_f1$Structure)
bfacs_f1$Structure <- ordered(bfacs_f1$Structure, levels=c("Open STING", "Closed STING"))

rm(open2, open2_db, close2, close2_db)

bfacs_f1_A <- ggplot(bfacs_f1[c(1:197, 395:591),]
                  , aes(x = res, y = nbfac)) + ylim(-1, 12.4) + a +
  geom_line(aes(color = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("forestgreen", "blue1")) +
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
  annotate("rect", xmin = 276, xmax=290, ymin= -1, ymax= 8, fill=NA, colour="black", size = .4) +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_f1_B <- ggplot(bfacs_f1[c(198:394, 592:788),], 
                  aes(x = res, y = nbfac)) + a + ylim(-1, 12.4) +
  geom_line(aes(color = Structure), size = 1) + xlab("Residues") + 
  scale_color_manual(values = c("forestgreen", "blue1")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax=487, ymin= -1, ymax= 8, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

png("figuras_paper/figura_1/F1.png", height = 90, width = 184, units = "mm", res = 300)
ggarrange(bfacs_f1_A, bfacs_f1_B, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "bottom")
dev.off()


##################################### Figura 4 ########################################
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

dists <- data.frame("time" = c(open2$time, open2$time, close2$time, close2$time),
                    "dist" = c(open2$chainA_norm, open2$chainB_norm, close2$chainA_norm, close2$chainB_norm),
                    "Structure" = "label")

dists$Structure[1:1001] <- "Open STING - Chain A"
dists$Structure[1002:2002] <- "Open STING - Chain B"
dists$Structure[2003:3003] <- "Closed STING - Chain A"
dists$Structure[3004:4004] <- "Closed STING - Chain B"

dists$Structure <- as.factor(dists$Structure)
dists$Structure <- ordered(dists$Structure, levels=c("Open STING - Chain A",
                                                     "Open STING - Chain B",
                                                     "Closed STING - Chain A", 
                                                     "Closed STING - Chain B"))


dist_open <- ggplot() +
  geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  geom_line(data = dists[1002:2002,], aes(x = time, y = dist), size = .8, color = "gray60") +
  geom_line(data = dists[1:1001,], aes(x = time, y = dist), size = .8, color = "forestgreen") + 
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) +
  geom_line(size = .8) + xlab("Tiempo (ns)") + ylab("Distancia (nm)")  +  a + ylim(-0.5, 0.5) 

dist_close <- ggplot() +
  geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  geom_line(data = dists[3004:4004,], aes(x = time, y = dist), size = .8, color = "gray50") +
  geom_line(data = dists[2003:3003,], aes(x = time, y = dist), size = .8, color = "blue1") + 
  ylab("Distancia (nm)") + xlab("Tiempo (ns)") + 
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) + a + ylim(-0.5, 0.5) 

dist_both <- ggplot(dists, aes(x = time, y = dist, color = Structure)) + 
  scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) + 
  geom_hline(yintercept = 0) +
  #geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7),
        axis.title.y = element_blank(), 
        axis.title.x = element_blank()) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)")  +  a + ylim(-0.5, 0.5) 

dist_boxplot <- ggplot(dists[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(x = time, y = dist, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray59")) +
  a +  theme(axis.title.y = element_blank(), 
             axis.title.x = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank()) 

figd <- ggplot(dists[c(501:1001, 1502:2002, 2503:3003,3504:4004),], 
               aes(dist, fill = Structure, color = Structure)) + coord_flip() +
  geom_density(alpha = .7) + a + ylab("Density") + theme(axis.title.y = element_blank()) + 
  scale_fill_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) +
  theme(axis.title = element_blank()) + 
  scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) + xlim(-0.5, 0.5)
  


abierta_rep2_sasa <- readXVG("abierta-rep2/sasa/abierta-rep2-sasa.xvg")
abierta_rep2_sasa <- as.data.frame(sapply(abierta_rep2_sasa, as.numeric))
cerrada_rep2_sasa <- readXVG("cerrada-rep2/sasa/cerrada-rep2-sasa.xvg")
cerrada_rep2_sasa <- as.data.frame(sapply(cerrada_rep2_sasa, as.numeric))

abierta_rep2_sasa$normA <- abierta_rep2_sasa$polsiteA - abierta_rep2_sasa$polsiteA[1]
abierta_rep2_sasa$normB <- abierta_rep2_sasa$polsiteB - abierta_rep2_sasa$polsiteB[1]
cerrada_rep2_sasa$normA <- cerrada_rep2_sasa$polsiteA - cerrada_rep2_sasa$polsiteA[1]
cerrada_rep2_sasa$normB <- cerrada_rep2_sasa$polsiteB - cerrada_rep2_sasa$polsiteB[1]

sasas <- data.frame("time" = c(abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                               cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time),
                    "sasa" = c(abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB,
                              cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB),
                    "Structure" = 'label')



sasas$Structure[1:1001] <- "Open STING - Chain A"
sasas$Structure[1002:2002] <- "Open STING - Chain B"
sasas$Structure[2003:3003] <- "Closed STING - Chain A"
sasas$Structure[3004:4004] <- "Closed STING - Chain B"

sasas$Structure <- as.factor(sasas$Structure)
sasas$Structure <- ordered(sasas$Structure, levels=c("Open STING - Chain A",
                                                     "Open STING - Chain B",
                                                     "Closed STING - Chain A", 
                                                     "Closed STING - Chain B"))


sasa_open <- ggplot() +
  geom_line(data = sasas[1002:2002,], aes(x = time, y = sasa), size = .8, color = "gray60") +
  geom_line(data = sasas[1:1001,], aes(x = time, y = sasa), size = .8, color = "forestgreen") + 
  xlab("Time (ns)") + labs(y=expression(SASA (nm^2))) + a  + ylim(6,10) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7),
        axis.title.x = element_blank()) 

sasa_close <- ggplot() +
  geom_line(data = sasas[3004:4004,], aes(x = time, y = sasa), size = .8, color = "gray60") +
  geom_line(data = sasas[2003:3003,], aes(x = time, y = sasa), size = .8, color = "blue1") + 
  geom_line(size = .8) + ylab(bquote(paste("SASA ("~nm^2~")"))) + a  + ylim(6,10) +
   theme(legend.background = element_blank(),
         axis.title.x = element_blank(),
                              legend.title = element_blank(),
                              legend.text = element_text(colour = 'black', size = 7))

sasas_boxplot <- ggplot(sasas[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(x = time, y = sasa, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray59")) +
  a +  theme(axis.title = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())  
  
png("figuras_paper/figura_4/distancias_oswald_alt.png", height = 60, width = 184, units = "mm", res = 300)

ggarrange(dist_both, figd, ncol = 2, nrow = 1, widths = c(2,1),
          common.legend = TRUE, legend = "bottom")

dev.off()


##################################### Figura 6 ########################################
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

bfacs_f6 <- data.frame("res" = c(v160m1_db$ATOM),
                       "nbfac" = c(v160m1_db$NBFAC),
                       "Structure" = "Open STING V160M")

bfacs_f6$Structure <- as.factor(bfacs_f6$Structure)

rm(v160m1_db, v160m1)

bfacs_f6_A <- ggplot(bfacs_f6[c(1:197),]
                     , aes(x = res, y = nbfac)) + ylim(-1, 8) + a +
  geom_line(color = "purple", size = 1) + ylab("Normalized B-factor") + 
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax=290, ymin= -1, ymax= 5, fill=NA, colour="black", size = .4) +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_f6_B <- ggplot(bfacs_f6[c(198:394),], 
                     aes(x = res, y = nbfac)) + a + ylim(-1, 8) +
  geom_line(color = "purple", size = 1) + xlab("Residues") + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax=487, ymin= -1, ymax= 5, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

png("figuras_paper/figura_6/bfacs.png", height = 70, width = 184, units = "mm", res = 300)
ggarrange(bfacs_f6_A, bfacs_f6_B, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "bottom")
dev.off()

v160m1_chainA <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_polsiteA.xvg")
v160m1_chainB <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_polsiteB.xvg")
v160m1_protein <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_protein.xvg")
v160m1 <- data.frame("time" = v160m1_chainA$time/1000, 
                     "chainA_x" = v160m1_chainA$x,"chainA_y" = v160m1_chainA$y, "chainA_z" = v160m1_chainA$z,
                     "chainB_x" = v160m1_chainB$x, "chainB_y" = v160m1_chainB$y,"chainB_z" = v160m1_chainB$z, 
                     "protein_x" = v160m1_protein$x,"protein_y" = v160m1_protein$y, "protein_z" = v160m1_protein$z)
v160m1$chainA_dist <- sqrt((v160m1$chainA_x - v160m1$protein_x)^2 + 
                             (v160m1$chainA_y - v160m1$protein_y)^2 + 
                             (v160m1$chainA_z - v160m1$protein_z)^2)
v160m1$chainB_dist <- sqrt((v160m1$chainB_x - v160m1$protein_x)^2 + 
                             (v160m1$chainB_y - v160m1$protein_y)^2 + 
                             (v160m1$chainB_z - v160m1$protein_z)^2)
v160m1$chainA_norm <- v160m1$chainA_dist - v160m1$chainA_dist[1]
v160m1$chainB_norm <- v160m1$chainB_dist - v160m1$chainB_dist[1]
rm(v160m1_chainA, v160m1_chainB, v160m1_protein)

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

dists_fig6 <- data.frame("time" = c(open2$time, open2$time, v160m1$time, v160m1$time),
                    "dist" = c(open2$chainA_norm, open2$chainB_norm,
                               v160m1$chainA_norm, v160m1$chainB_norm),
                    "Structure" = "label")

dists_fig6$Structure[1:1001] <- "Open WT STING - Chain A"
dists_fig6$Structure[1002:2002] <- "Open WT STING - Chain B"
dists_fig6$Structure[2003:3003] <- "APO STING V160M Cadena A"
dists_fig6$Structure[3004:4004] <- "APO STING V160M Cadena B"
dists_fig6$Structure <- as.factor(dists_fig6$Structure)
dists_fig6$Structure <- ordered(dists_fig6$Structure, levels=c("Open WT STING - Chain A", 
                                                               "Open WT STING - Chain B",
                                                               "APO STING V160M Cadena A",
                                                               "APO STING V160M Cadena B"))

dist_open_f6 <- ggplot(dists_fig6[c(1:2002),], aes(x = time, y = dist, color = Structure)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8, alpha = .8) + scale_color_manual(values = c("forestgreen", "gray60")) +
  xlab("Time (ns)") + ylab("Distance (nm)")  +  a + ylim(-0.5, 0.5) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) 
dist_v160m <- ggplot(dists_fig6[c(2003:4004),], aes(x = time, y = dist, color = Structure)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8, alpha = .8) + scale_color_manual(values = c("purple", "gray60")) +
  xlab("Time (ns)") + ylab("Distance (nm)")  +  a + ylim(-0.5, 0.5) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) 

dist_boxplot_f6 <- ggplot(dists_fig6[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(x = time, y = dist, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "purple", "gray59")) +
  a +  theme(axis.title.y = element_blank(), 
             axis.title.x = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())

figd <- ggplot(dists_fig6[c(2503:3003,3504:4004),], 
               aes(dist, fill = Structure, color = Structure)) + 
  geom_density(alpha = .7) + a + ylab("Density") + xlab("Distance (nm)") + 
  scale_fill_manual(values = c("purple", "gray60")) + xlim(-0.5, 0.5) +
  scale_color_manual(values = c("purple", "gray60")) + coord_flip() 

png("figuras_tesis/dists_open_mut.png", height =90, width = 184, units = "mm", res = 300)
ggarrange(dist_v160m, dist_open_f6, nrow = 2, ncol = 1, 
          common.legend = TRUE, legend = "bottom")
dev.off()

png("figuras_paper/figura_6/distance_mut1.png", height =60, width = 184, units = "mm", res = 300)
dist_v160m
dev.off()

v160m_rep1_sasa <- readXVG("v160m-rep1/prod_200ns/sasa/v160m-rep1-sasa.xvg")
v160m_rep1_sasa <- as.data.frame(sapply(v160m_rep1_sasa, as.numeric))
abierta_rep2_sasa <- readXVG("abierta-rep2/sasa/abierta-rep2-sasa.xvg")
abierta_rep2_sasa <- as.data.frame(sapply(abierta_rep2_sasa, as.numeric))


sasas_f6 <- data.frame("time" = c(abierta_rep2_sasa$Time, abierta_rep2_sasa$Time, 
                                  v160m_rep1_sasa$Time, v160m_rep1_sasa$Time),
                    "sasa" = c(abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB,
                               v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB),
                    "Structure" = 'label')

sasas_f6$Structure[1:1001] <- "Open STING - Chain A"
sasas_f6$Structure[1002:2002] <- "Open STING - Chain B"
sasas_f6$Structure[2003:3003] <- "Open V160M STING - Chain A"
sasas_f6$Structure[3004:4004] <- "Open WT & V160M STING - Chain B"

sasas_f6$Structure <- as.factor(sasas_f6$Structure)
sasas_f6$Structure <- ordered(sasas_f6$Structure, levels=c("Open STING - Chain A", 
                                                           "Open STING - Chain B",
                                                           "Open V160M STING - Chain A",
                                                           "Open WT & V160M STING - Chain B"))

sasa_v160m <- ggplot(sasas_f6, aes(x = time/1000, y = sasa, color = Structure)) +
  geom_line(size = .8, alpha = .8) + scale_color_manual(values = c("purple", "gray60")) +
  xlab("Time (ns)") + ylab(bquote(paste("SASA ("~nm^2~")"))) +  a +
  theme(axis.title.x = element_blank(),legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) 

sasas_f6_boxplot <- ggplot(sasas_f6[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(y = sasa, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "purple", "gray59")) +
  a +  theme(axis.title = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank()) 

png("figuras_paper/figura_6/dist_sasa.png", height = 120, width = 184, units = "mm", res = 300)
ggarrange(dist_v160m, sasa_v160m, nrow = 2, ncol = 1, common.legend = TRUE, 
          legend = "bottom")
dev.off()

png("figuras_paper/figura_6/sasa.png", height = 60, width = 92, units = "mm", res = 300)
ggarrange(sasas_f6_boxplot, nrow = 1, ncol = 1, common.legend = TRUE, 
          legend = "bottom")
dev.off()

png("figuras_paper/figura_6/for_label.png", height = 60, width = 184, units = "mm", res = 300)
ggarrange(sasas_f6_boxplot, nrow = 1, ncol = 1, common.legend = TRUE, 
          legend = "bottom")
dev.off()

###

##################################### Figura S1 #######################################
open1_rmsd <- read.csv("abierta-rep1/crop-200/abierta-rep1-rmsd.xvg", sep="")
open2_rmsd <- read.csv("abierta-rep2/abierta-rep2-rmsd.xvg", sep="")
close1_rmsd <- read.csv("cerrada-rep1/prod/control_rmsd.xvg", sep="")
close2_rmsd <- read.csv("cerrada-rep2/cerrada-rep2-rmsd.xvg", sep="")
v160m1_rmsd <- read.csv("v160m-rep1/prod_200ns/abierta_v160m_200ns_rmsd.xvg", sep="")
v160m2_rmsd <- read.csv("v160m-rep2/abierta-v160m-rep2-rmsd.xvg", sep="")

open1_rg <- read.csv("abierta-rep1/crop-200/abierta-rep1-gr.xvg", sep="")
open2_rg <- read.csv("abierta-rep2/abierta-rep2-gr.xvg", sep="")
close1_rg <- read.csv("cerrada-rep1/prod/control_rg.xvg", sep="")
close2_rg <- read.csv("cerrada-rep2/cerrada-rep2-rg.xvg", sep="")
v160m1_rg <- read.csv("v160m-rep1/prod_200ns/abierta_v160m_200ns_rg.xvg", sep="")
v160m2_rg <- read.csv("v160m-rep2/abierta-v160m-rep2-gr.xvg", sep="")

rmsds <- data.frame("time" = c(open1_rmsd$time, open2_rmsd$time, 
                              close1_rmsd$time, close2_rmsd$time, 
                              v160m1_rmsd$time, v160m2_rmsd$time),
                   "rmsd" = c(open1_rmsd$rmsd, open2_rmsd$rmsd,
                              close1_rmsd$rmsd, close2_rmsd$rmsd,
                              v160m1_rmsd$rmsd, v160m2_rmsd$rmsd),
                   "Structure" = "label")

rmsds$Structure[1:1001] <- "APO STING rep 1" #palegreen
rmsds$Structure[1002:2002] <- "APO STING rep 2" #limegreen
rmsds$Structure[2003:3003] <- "HOLO STING rep 1" #blue
rmsds$Structure[3004:4004] <- "HOLO STING rep 2" #dodgerblue
rmsds$Structure[4005:5005] <-  "APO STING V160M rep 1" #deeppink
rmsds$Structure[5006:6006] <-  "APO STING V160M rep 2" #hotpink1


corridas <- c("APO STING rep 1", "APO STING rep 2",
              "HOLO STING rep 1", "HOLO STING rep 2",
              "APO STING V160M rep 1", "APO STING V160M rep 2")

rmsds$Structure <- as.factor(rmsds$Structure)
rmsds$Structure <- ordered(rmsds$Structure, levels=corridas)
rm(open1_rmsd, open2_rmsd, close1_rmsd, close2_rmsd, v160m1_rmsd, v160m2_rmsd)

gr <- data.frame("time" = c(open1_rg$time, open2_rg$time,
                            close1_rg$time, close2_rg$time,
                            v160m1_rg$time, v160m2_rg$time),
                 "gr" = c(open1_rg$all, open2_rg$all,
                          close1_rg$all, close2_rg$all,
                          v160m1_rg$all, v160m2_rg$all),
                 "Structure" = "label")

gr$Structure[1:1001] <- "APO STING rep 1" #palegreen
gr$Structure[1002:2002] <- "APO STING rep 2" #limegreen
gr$Structure[2003:3003] <- "HOLO STING rep 1" #blue
gr$Structure[3004:4004] <- "HOLO STING rep 2" #dodgerblue
gr$Structure[4005:5005] <-  "APO STING V160M rep 1" #deeppink
gr$Structure[5006:6006] <-  "APO STING V160M rep 2" #hotpink1

rm(open1_rg, open2_rg,close1_rg, close2_rg, v160m1_rg, v160m2_rg)

gr$Structure <- as.factor(gr$Structure)
gr$Structure <- ordered(gr$Structure, levels=corridas)

grs <- ggplot(gr, aes(x = time/1000, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = Structure), size = 0.8) + 
  scale_color_manual(values = c("mediumseagreen", "forestgreen",
                                "dodgerblue3", "blue1",
                                "magenta3", "magenta1")) +
  xlab("Time (ns)") + ylab(expression(R["g"] (nm))) + a + 
  ggtitle("Radius of gyration") 

rmsdss <- ggplot(rmsds, aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure), size = 0.8) + 
  scale_color_manual(values = c("mediumseagreen", "forestgreen",
                                "dodgerblue3", "blue1",
                                "magenta3", "magenta1")) +
  ylab("RMSD (nm)") + a + ggtitle("Root mean square deviation") + 
  theme(axis.title.x = element_blank())

png("figuras_paper/figura_s1/FS1.png", height = 120, width = 184, units = "mm", res = 300)
ggarrange(rmsdss, grs, nrow = 2, ncol = 1, legend = "bottom", common.legend = TRUE)
dev.off()

rmsd_tesis <- ggplot(rmsds[4005:6006,], aes(x = time, y = rmsd)) + 
  geom_line(aes(color = Structure), size = 0.8) + 
  scale_color_manual(values = c("magenta3", "magenta1")) +
  ylab("RMSD (nm)") + a + theme(axis.title.x = element_blank()) +
  guides(shape = guide_legend(override.aes = list(size = 0.5)))

grs_tesis <- ggplot(gr[4005:6006,], aes(x = time/1000, y = gr)) + ylim(2,2.5) +
  geom_line(aes(color = Structure), size = 0.8) + 
  scale_color_manual(values = c("magenta3", "magenta1")) +
  xlab("Tiempo (ns)") + ylab(expression(R["g"] (nm))) + a + 
  guides(shape = guide_legend(override.aes = list(size = 0.5)))

png("figuras_tesis/RMSD_gr_mut.png", height = 120, width = 184, units = "mm", res = 300)

ggarrange(rmsd_tesis, grs_tesis, nrow = 2, ncol = 1, legend = "bottom", common.legend = TRUE) 
dev.off()



##################################### Figura S2 #######################################
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

bfacs_fs2 <- data.frame("res" = c(open1_db$ATOM, close1_db$ATOM),
                       "nbfac" = c(open1_db$NBFAC,close1_db$NBFAC),
                       "Structure" = "label")

bfacs_fs2$Structure[1:394] <- "Open STING" #palegreen
bfacs_fs2$Structure[395:788] <- "Closed STING" #limegreen

bfacs_fs2$Structure <- as.factor(bfacs_fs2$Structure)
bfacs_fs2$Structure <- ordered(bfacs_fs2$Structure, levels=c("Open STING", "Closed STING"))

rm(open1, open1_db, close1, close1_db)

bfacs_fs2_A <- ggplot(bfacs_fs2[c(1:197, 395:591),]
                     , aes(x = res, y = nbfac)) + ylim(-1, 12.4) + a +
  geom_line(aes(color = Structure), size = 1) + ylab("Normalized B-factor") + 
  scale_color_manual(values = c("forestgreen", "blue1")) +
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
  annotate("rect", xmin = 276, xmax=290, ymin= -1, ymax= 8, fill=NA, colour="black", size = .4) +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_fs2_B <- ggplot(bfacs_fs2[c(198:394, 592:788),], 
                     aes(x = res, y = nbfac)) + a + ylim(-1, 12.4) +
  geom_line(aes(color = Structure), size = 1) + xlab("Residues") + 
  scale_color_manual(values = c("forestgreen", "blue1")) +
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=12.4, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax=487, ymin= -1, ymax= 8, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

png("figuras_paper/figura_s2/Fs2.png", height = 90, width = 184, units = "mm", res = 300)
ggarrange(bfacs_fs2_A, bfacs_fs2_B, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "bottom")
dev.off()


##################################### Figura S4 #######################################
abierta1_pca <- readXVG("abierta-rep1/crop-200/fel/new/eigenval.xvg")
abierta1_pca <- as.data.frame(sapply(abierta1_pca, as.numeric))
colnames(abierta1_pca) <- c("pc", "eigenval")
abierta1_pca$contribution <- abierta1_pca$eigenval/sum(abierta1_pca$eigenval)*100
abierta1_pca$accumulated = 0
abierta1_pca$accumulated[1] = abierta1_pca$contribution[1]
for (i in 2:500){
  abierta1_pca$accumulated[i] = abierta1_pca$contribution[i] + abierta1_pca$accumulated[i - 1]
}
abierta1_pca$text <- paste(as.character(round(abierta1_pca$accumulated, 1)), "%", sep = "")
abierta1_pca_plot <- ggplot(abierta1_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "forestgreen") + a + theme(axis.title = element_blank(),
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank()) +
  geom_text(aes(y= abierta1_pca$contribution[c(1:7)] + 2, 
                label = abierta1_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7")))

abierta2_pca <- readXVG("abierta-rep2/fel/new/eigenval.xvg")
abierta2_pca <- as.data.frame(sapply(abierta2_pca, as.numeric))
colnames(abierta2_pca) <- c("pc", "eigenval")
abierta2_pca$contribution <- abierta2_pca$eigenval/sum(abierta2_pca$eigenval)*100
abierta2_pca$accumulated = 0
abierta2_pca$accumulated[1] = abierta2_pca$contribution[1]
for (i in 2:500){
  abierta2_pca$accumulated[i] = abierta2_pca$contribution[i] + abierta2_pca$accumulated[i - 1]
}
abierta2_pca$text <- paste(as.character(round(abierta2_pca$accumulated, 1)), "%", sep = "")
abierta2_pca_plot <- ggplot(abierta2_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "forestgreen") + a + theme(axis.title = element_blank(), 
                                             panel.grid.major = element_blank(),
                                             panel.grid.minor = element_blank()) +
  geom_text(aes(y= abierta2_pca$contribution[c(1:7)] + 2, 
                label = abierta2_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7")))

cerrada1_pca <- readXVG("cerrada-rep1/prod/fel/new/eigenval.xvg")
cerrada1_pca <- as.data.frame(sapply(cerrada1_pca, as.numeric))
colnames(cerrada1_pca) <- c("pc", "eigenval")
cerrada1_pca$contribution <- cerrada1_pca$eigenval/sum(cerrada1_pca$eigenval)*100
cerrada1_pca$accumulated = 0
cerrada1_pca$accumulated[1] = cerrada1_pca$contribution[1]
for (i in 2:500){
  cerrada1_pca$accumulated[i] = cerrada1_pca$contribution[i] + cerrada1_pca$accumulated[i - 1]
}
cerrada1_pca$text <- paste(as.character(round(cerrada1_pca$accumulated, 1)), "%", sep = "")
cerrada1_pca_plot <- ggplot(cerrada1_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "blue1") + a + theme(axis.title = element_blank(), 
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank()) +
  geom_text(aes(y= cerrada1_pca$contribution[c(1:7)] + 2, 
                label = cerrada1_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7")))

cerrada2_pca <- readXVG("cerrada-rep2/fel/new/eigenval.xvg")
cerrada2_pca <- as.data.frame(sapply(cerrada2_pca, as.numeric))
colnames(cerrada2_pca) <- c("pc", "eigenval")
cerrada2_pca$contribution <- cerrada2_pca$eigenval/sum(cerrada2_pca$eigenval)*100
cerrada2_pca$accumulated = 0
cerrada2_pca$accumulated[1] = cerrada2_pca$contribution[1]
for (i in 2:500){
  cerrada2_pca$accumulated[i] = cerrada2_pca$contribution[i] + cerrada2_pca$accumulated[i - 1]
}
cerrada2_pca$text <- paste(as.character(round(cerrada2_pca$accumulated, 1)), "%", sep = "")
cerrada2_pca_plot <- ggplot(cerrada2_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "blue1") + a + theme(axis.title = element_blank(), 
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank()) +
  geom_text(aes(y= cerrada2_pca$contribution[c(1:7)] + 2, 
                label = cerrada2_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("CP1", "CP2", "CP3", "CP4", "CP5", "CP6", "CP7")))

mutante1_pca <- readXVG("v160m-rep1/prod_200ns/fel/new/eigenval.xvg")
mutante1_pca <- as.data.frame(sapply(mutante1_pca, as.numeric))
colnames(mutante1_pca) <- c("pc", "eigenval")
mutante1_pca$contribution <- mutante1_pca$eigenval/sum(mutante1_pca$eigenval)*100
mutante1_pca$accumulated = 0
mutante1_pca$accumulated[1] = mutante1_pca$contribution[1]
for (i in 2:500){
  mutante1_pca$accumulated[i] = mutante1_pca$contribution[i] + mutante1_pca$accumulated[i - 1]
}
mutante1_pca$text <- paste(as.character(round(mutante1_pca$accumulated, 1)), "%", sep = "")
mutante1_pca_plot <- ggplot(mutante1_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "purple") + a + theme(axis.title = element_blank(), 
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
  geom_text(aes(y= mutante1_pca$contribution[c(1:7)] + 2, 
                label = mutante1_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")))

mutante2_pca <- readXVG("v160m-rep2/fel/new/eigenval.xvg")
mutante2_pca <- as.data.frame(sapply(mutante2_pca, as.numeric))
colnames(mutante2_pca) <- c("pc", "eigenval")
mutante2_pca$contribution <- mutante2_pca$eigenval/sum(mutante2_pca$eigenval)*100
mutante2_pca$accumulated = 0
mutante2_pca$accumulated[1] = mutante2_pca$contribution[1]
for (i in 2:500){
  mutante2_pca$accumulated[i] = mutante2_pca$contribution[i] + mutante2_pca$accumulated[i - 1]
}
mutante2_pca$text <- paste(as.character(round(mutante2_pca$accumulated, 1)), "%", sep = "")
mutante2_pca_plot <- ggplot(mutante2_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "purple") + a + theme(axis.title = element_blank(), 
                                        panel.grid.major = element_blank(),
                                        panel.grid.minor = element_blank()) +
  geom_text(aes(y= mutante2_pca$contribution[c(1:7)] + 2, 
                label = mutante2_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) +
  scale_x_continuous(breaks = c(1:7), labels = paste(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")))

png("figuras_paper/figura_s4/FS4.png", height = 150, width = 184, units = "mm", res = 300)
ggarrange(abierta1_pca_plot, abierta2_pca_plot,
          cerrada1_pca_plot, cerrada2_pca_plot,
          mutante1_pca_plot, mutante2_pca_plot,
          nrow = 3, ncol = 2)
dev.off()
png("figuras_paper/figura_s4/pic_for_axis.png", height = 50, width = 92, units = "mm", res = 300)
ggplot(mutante2_pca[c(1:7),], aes(x = pc, y = contribution)) + 
  geom_col(fill = "purple") + a + xlab("Principal Component") + ylab("Contribution (%)") + 
  geom_text(aes(y= mutante2_pca$contribution[c(1:7)] + 2, 
                label = mutante2_pca$text[c(1:7)]), 
            color = "black", size = 3, hjust=0.5) 
dev.off()


png("figuras_tesis/pca_contribucion.png", height = 60, width = 184, units = "mm", res = 300)
ggarrange(abierta2_pca_plot,cerrada2_pca_plot, nrow = 1, ncol = 2)
dev.off()

png("figuras_tesis/pca_contribucion_supp.png", height = 60, width = 184, units = "mm", res = 300)
ggarrange(abierta1_pca_plot,cerrada1_pca_plot, nrow = 1, ncol = 2)
dev.off()

png("figuras_tesis/pca_contribucion_mut.png", height = 60, width = 184, units = "mm", res = 300)
ggarrange(mutante1_pca_plot,mutante2_pca_plot, nrow = 1, ncol = 2)
dev.off()

####
##################################### Figura S6 #######################################
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

dists2 <- data.frame("time" = c(open1$time, open1$time, close1$time, close1$time),
                    "dist" = c(open1$chainA_norm, open1$chainB_norm, close1$chainA_norm, close1$chainB_norm),
                    "Structure" = "label")

dists2$Structure[1:1001] <- "Open STING - Chain A"
dists2$Structure[1002:2002] <- "Open STING - Chain B"
dists2$Structure[2003:3003] <- "Closed STING - Chain A"
dists2$Structure[3004:4004] <- "Closed STING - Chain B"

dists2$Structure <- as.factor(dists2$Structure)
dists2$Structure <- ordered(dists2$Structure, levels=c("Open STING - Chain A",
                                                     "Open STING - Chain B",
                                                     "Closed STING - Chain A", 
                                                     "Closed STING - Chain B"))


dist_open_s6 <- ggplot() +
  geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  geom_line(data = dists2[1002:2002,], aes(x = time, y = dist), size = .8, color = "gray50") +
  geom_line(data = dists2[1:1001,], aes(x = time, y = dist), size = .8, color = "forestgreen") + 
  theme(axis.title.x = element_blank(),legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab(" ")  +  a + ylim(-0.5, 0.5)

dist_close_s6 <- ggplot() +
  geom_line(aes(x = c(0:200), y = 0), color = "black", alpha = 1, linetype = 1, size = .8) +
  geom_line(data = dists2[3004:4004,], aes(x = time, y = dist), size = .8, color = "gray60") +
  geom_line(data = dists2[2003:3003,], aes(x = time, y = dist), size = .8, color = "blue1") + 
  ylab("Normalized Distance (nm)") + 
  theme(axis.title.x = element_blank(), legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) + a + ylim(-0.5, 0.5) 

dist_boxplot_s6 <- ggplot(dists2[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(x = time, y = dist, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray59")) +
  a +  theme(axis.title.y = element_blank(), 
             axis.title.x = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank()) 

density_dist <- ggplot(dists2[c(501:1001, 1502:2002, 2503:3003,3504:4004),], 
               aes(dist, fill = Structure, color = Structure)) +  coord_flip() +
  geom_density(alpha = .7) + a + ylab("Density") + theme(axis.title.x = element_blank()) + 
  scale_fill_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) +
  scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray50")) + xlim(-0.5, 0.5)


abierta_rep1_sasa <- readXVG("abierta-rep1/crop-200/sasa/abierta-rep1-sasa.xvg")
abierta_rep1_sasa <- as.data.frame(sapply(abierta_rep1_sasa, as.numeric))
cerrada_rep1_sasa <- readXVG("cerrada-rep1/prod/sasa/cerrada-rep1-sasa.xvg")
cerrada_rep1_sasa <- as.data.frame(sapply(cerrada_rep1_sasa, as.numeric))

sasas2 <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                               cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time),
                    "sasa" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                               cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB),
                    "Structure" = 'label')



sasas2$Structure[1:1001] <- "Open STING - Chain A"
sasas2$Structure[1002:2002] <- "Open STING - Chain B"
sasas2$Structure[2003:3003] <- "Closed STING - Chain A  "
sasas2$Structure[3004:4004] <- "Open & Closed STING - Chain B"

sasas2$Structure <- as.factor(sasas2$Structure)
sasas2$Structure <- ordered(sasas2$Structure, levels=c("Open STING - Chain A",
                                                     "Open STING - Chain B",
                                                     "Closed STING - Chain A  ", 
                                                     "Open & Closed STING - Chain B"))


sasa_open_fs6 <- ggplot() +
  geom_line(data = sasas2[1002:2002,], aes(x = time, y = sasa), size = .8, color = "gray60") +
  geom_line(data = sasas2[1:2002,], aes(x = time, y = sasa), size = .8, color = "forestgreen") + 
  xlab("Time (ns)") + labs(y=expression(SASA (nm^2))) + a  + ylim(6,10) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7),
        axis.title.x = element_blank()) 

sasa_close_fs6 <- ggplot() +
  geom_line(data = sasas2[3004:4004,], aes(x = time, y = sasa), size = .8, color = "gray60") +
  geom_line(data = sasas2[2003:3003,], aes(x = time, y = sasa), size = .8, color = "blue1") + 
  geom_line(size = .8) + ylab(bquote(paste("SASA ("~nm^2~")"))) + a  + ylim(6,10) +
  theme(legend.background = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7))

sasas_boxplot_fs6 <- ggplot(sasas2[c(501:1001, 1502:2002, 2503:3003,3504:4004),], aes(x = time, y = sasa, color = Structure)) + 
  geom_boxplot() + scale_color_manual(values = c("forestgreen", "gray60", "blue1", "gray59")) +
  a +  theme(axis.title = element_blank(),
             legend.background = element_blank(),
             legend.title = element_blank(), 
             legend.text = element_text(colour = 'black', size = 7),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())  

png("figuras_paper/figura_s6/FS6.png", height = 180, width = 184, units = "mm", res = 300)

ggarrange(ggarrange(dist_open_s6, dist_close_s6, sasa_open_fs6, sasa_close_fs6,
                    ncol = 1, nrow = 4, common.legend = TRUE, legend = "none"),
          ggarrange(dist_boxplot_s6, sasas_boxplot_fs6, ncol = 1, nrow = 2, legend = "none"),
          common.legend = TRUE, legend = "none", ncol = 2, nrow = 1, widths = c(2,1))

dev.off()

png("figuras_paper/figura_s6/dists_supp.png", height = 60, width = 184, units = "mm", res = 300)
ggarrange(ggarrange(dist_open_s6, dist_close_s6,
                    ncol = 1, nrow = 2, common.legend = TRUE, legend = "none"),
          density_dist,
          common.legend = TRUE, legend = "none", ncol = 2, nrow = 1, widths = c(2,1))

dev.off()

##################################### Figura S8 #######################################
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

bfacs_fs8 <- data.frame("res" = c(v160m2_db$ATOM),
                       "nbfac" = c(v160m2_db$NBFAC),
                       "Structure" = "Open STING V160M")

bfacs_fs8$Structure <- as.factor(bfacs_fs8$Structure)

rm(v160m2_db, v160m2)

bfacs_fs8_A <- ggplot(bfacs_fs8[c(1:197),]
                     , aes(x = res, y = nbfac)) + ylim(-1, 8) + a +
  geom_line(color = "purple", size = 1) + ylab("Normalized B-factor") + 
  scale_x_continuous(breaks = c(146, 156, 166, 176, 186, 196, 206, 216, 226, 236, 
                                246, 256, 266, 276, 286, 296, 306, 316, 326, 336),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=147, y=Inf, label="Chain A", hjust=-.2, vjust=2, fontface="bold") + 
  annotate("rect", xmin = 189, xmax=200, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 232, xmax=244, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 254, xmax=258, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 322, xmax=328, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 276, xmax=290, ymin= -1, ymax= 5, fill=NA, colour="black", size = .4) +
  theme(
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
bfacs_fs8_B <- ggplot(bfacs_fs8[c(198:394),], 
                     aes(x = res, y = nbfac)) + a + ylim(-1, 8) +
  geom_line(color = "purple", size = 1) + xlab("Residues") + 
  scale_x_continuous(breaks = c(343, 353, 363, 373, 383, 393, 403, 413, 423, 433, 
                                443, 453, 463, 473, 483, 493, 503, 513, 523, 533),
                     labels = paste(c("ALA146", "SER156", "TRP166", "VAL176", "GLU186", 
                                      "ALA196", "LEU206", "ASP216", "TYR226", "THR236", 
                                      "LYS246", "ASP256", "PHE266", "MET276", "ARG286", 
                                      "PHE296", "SER306", "LEU316", "GLU326", "TRP336"))) +
  annotate("text", x=-Inf, y=Inf, label="Chain B", hjust=-.2, vjust=2, fontface="bold") +
  annotate("rect", xmin = 386, xmax=397, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 429, xmax=441, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 451, xmax=455, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") + 
  annotate("rect", xmin = 519, xmax=525, ymin= -0.8, ymax=8, alpha=.2, fill="grey60") +
  annotate("rect", xmin = 473, xmax=487, ymin= -1, ymax= 5, fill=NA, colour="black", size = .4) +
  theme(axis.text.x = element_text(angle = 90, face = "bold", size = 7),
        panel.grid.major = element_line(color = NA), 
        panel.grid.minor = element_line(color = NA),
        plot.title = element_text(colour = 'black', size = 15),
        axis.title = element_text(colour = 'black', face = 'italic'),
        panel.background = element_rect(fill = "white", colour = "black",),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 

png("figuras_paper/figura_s8/bfacs.png", height = 70, width = 184, units = "mm", res = 300)
ggarrange(bfacs_fs8_A, bfacs_fs8_B, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "bottom")
dev.off()

v160m2_chainA <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-polA.xvg")
v160m2_chainB <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-polB.xvg")
v160m2_protein <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-protein.xvg")

v160m2 <- data.frame("time" = v160m2_chainA$time/1000, 
                     "chainA_x" = v160m2_chainA$x,"chainA_y" = v160m2_chainA$y, "chainA_z" = v160m2_chainA$z,
                     "chainB_x" = v160m2_chainB$x, "chainB_y" = v160m2_chainB$y,"chainB_z" = v160m2_chainB$z, 
                     "protein_x" = v160m2_protein$x,"protein_y" = v160m2_protein$y, "protein_z" = v160m2_protein$z)
v160m2$chainA_dist <- sqrt((v160m2$chainA_x - v160m2$protein_x)^2 + 
                             (v160m2$chainA_y - v160m2$protein_y)^2 + 
                             (v160m2$chainA_z - v160m2$protein_z)^2)
v160m2$chainB_dist <- sqrt((v160m2$chainB_x - v160m2$protein_x)^2 + 
                             (v160m2$chainB_y - v160m2$protein_y)^2 + 
                             (v160m2$chainB_z - v160m2$protein_z)^2)
v160m2$chainA_norm <- v160m2$chainA_dist - v160m2$chainA_dist[1]
v160m2$chainB_norm <- v160m2$chainB_dist - v160m2$chainB_dist[1]
rm(v160m2_chainA, v160m2_chainB, v160m2_protein)

dists_figs8 <- data.frame("time" = c(v160m2$time, v160m2$time),
                         "dist" = c(v160m2$chainA_norm, v160m2$chainB_norm),
                         "Structure" = "label")

dists_figs8$Structure[1:1001] <- "Chain A"
dists_figs8$Structure[1002:2002] <- "Chain B"
dists_figs8$Structure <- as.factor(dists_figs8$Structure)
dists_figs8$Structure <- ordered(dists_figs8$Structure, levels=c("Chain A", "Chain B"))

dist_v160m_fs8 <- ggplot(dists_figs8, aes(x = time, y = dist, color = Structure)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8, alpha = .8) + scale_color_manual(values = c("purple", "gray60")) +
  xlab("Time (ns)") + ylab("Distance (nm)")  +  a + ylim(-0.5, 0.5) +
  theme(legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) 

figd <- ggplot(dists_figs8[c(501:1001,1502:2002),], 
               aes(dist, fill = Structure, color = Structure)) + 
  geom_density(alpha = .7) + a + ylab("Density") + xlab(" ") + 
  scale_fill_manual(values = c("purple", "gray60")) + coord_flip() +
  scale_color_manual(values = c("purple", "gray60")) + xlim(-0.5, 0.5)

png("figuras_paper/figura_s8/distance_mut2.png", height =60, width = 184, units = "mm", res = 300)
ggarrange(dist_v160m_fs8, figd, nrow =1, ncol = 2, common.legend = TRUE,  widths = c(2,1), legend = "bottom")
dev.off()


v160m_rep2_sasa <- readXVG("v160m-rep2/sasa/v160m-rep2-sasa.xvg")
v160m_rep2_sasa <- as.data.frame(sapply(v160m_rep2_sasa, as.numeric))

sasas_fs8 <- data.frame("time" = c(v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                       "sasa" = c(v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                       "Structure" = 'label')

sasas_fs8$Structure[1:1001] <- "Chain A"
sasas_fs8$Structure[1002:2002] <- "Chain B"

sasas_fs8$Structure <- as.factor(sasas_fs8$Structure)
sasas_fs8$Structure <- ordered(sasas_fs8$Structure, levels=c("Chain A", "Chain B"))

sasa_v160m_fs8 <- ggplot(sasas_fs8, aes(x = time, y = sasa, color = Structure)) +
  geom_line(size = .8, alpha = .8) + scale_color_manual(values = c("purple", "gray60")) +
  xlab("Time (ns)") + ylab(bquote(paste("SASA ("~nm^2~")"))) +  a +
  theme(axis.title.x = element_blank(),legend.background = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black', size = 7)) 

png("figuras_paper/figura_s8/dist_sasa.png", height = 120, width = 184, units = "mm", res = 300)
ggarrange(dist_v160m_fs8, sasa_v160m_fs8, nrow = 2, ncol = 1, common.legend = TRUE, 
          legend = "bottom")
dev.off()


##################################### Figura Tesis Gly163 #############################
abierta_rep1_gly163 <- readXVG("abierta-rep1/crop-200/gly163/abierta_rep1_gly163.xvg")
abierta_rep1_gly163 <- as.data.frame(sapply(abierta_rep1_gly163, as.numeric))

abierta_rep2_gly163 <- readXVG("abierta-rep2/gly163/abierta_rep2_gly163.xvg")
abierta_rep2_gly163 <- as.data.frame(sapply(abierta_rep2_gly163, as.numeric))

cerrada_rep1_gly163 <- readXVG("cerrada-rep1/prod/gly163/cerrada_rep1_gly163.xvg")
cerrada_rep1_gly163 <- as.data.frame(sapply(cerrada_rep1_gly163, as.numeric))

cerrada_rep2_gly163 <- readXVG("cerrada-rep2/gly163/cerrada_rep2_gly163.xvg")
cerrada_rep2_gly163 <- as.data.frame(sapply(cerrada_rep2_gly163, as.numeric))

gly163 <- data.frame("time" = c(abierta_rep1_gly163$Time, abierta_rep2_gly163$Time,
                                cerrada_rep1_gly163$Time, cerrada_rep2_gly163$Time),
                     "dist" =  c(abierta_rep1_gly163[,2] , abierta_rep2_gly163[,2] ,
                                 cerrada_rep1_gly163[,2] , cerrada_rep2_gly163[,2] ),
                     "E" = "label")

gly163$E[1:1001] <- "APO STING rep 1"
gly163$E[1002:2002] <- "APO STING rep 2"
gly163$E[2003:3003] <- "HOLO STING rep 1"
gly163$E[3004:4004] <- "HOLO STING rep 2"

corridas <- c("APO STING rep 1", "APO STING rep 2",
              "HOLO STING rep 1", "HOLO STING rep 2")

gly163$E <- as.factor(gly163$E)
gly163$E <- ordered(gly163$E, levels=corridas)

GLY_tesis <- ggplot(gly163, aes(x = time, y = dist*10)) + 
  geom_line(aes(color = E), size = 0.8, alpha = .8) + xlab("Tiempo (ns)") + 
  scale_color_manual(values = c("mediumseagreen", "forestgreen",
                                "dodgerblue3", "blue1")) +
  ylab("Distancia (Å)") + a + theme(axis.title.x = element_blank()) +
  scale_y_continuous(breaks = c(4, 6, 8, 10)) 

png("figuras_tesis/gly163.png", height = 60, width = 184, units = "mm", res = 300)

ggarrange(GLY_tesis, nrow = 1, ncol = 1, legend = "bottom", common.legend = TRUE) 
dev.off()

##################################### Figura Asimetria ###############################
abierta_rep1_rmsd <- read.table("~/Desktop/Sysbio/cSTING/asymmetry/abierta_rep1_rmsd.txt", quote="\"", comment.char="")
abierta_rep2_rmsd <- read.table("~/Desktop/Sysbio/cSTING/asymmetry/abierta_rep2_rmsd.txt", quote="\"", comment.char="")
close_rep1_rmsd <- read.table("~/Desktop/Sysbio/cSTING/asymmetry/cerrada_rep1_rmsd.txt", quote="\"", comment.char="")
close_rep2_rmsd <- read.table("~/Desktop/Sysbio/cSTING/asymmetry/cerrada_rep2_rmsd.txt", quote="\"", comment.char="")
open2_chainA <-read.delim("abierta-rep2/distance/abierta2_chainA.xvg")

asymmetry <- data.frame("time" = c(open2_chainA$time, open2_chainA$time,
                                   open2_chainA$time, open2_chainA$time),
                        "rmsd" = c(abierta_rep1_rmsd$V9, abierta_rep2_rmsd$V9, 
                                   close_rep1_rmsd$V9, close_rep2_rmsd$V9), 
                        "E" = "label")

asymmetry$E[1:1001] <- "APO STING rep 1"
asymmetry$E[1002:2002] <- "APO STING rep 2"
asymmetry$E[2003:3003] <- "HOLO STING rep 1"
asymmetry$E[3004:4004] <- "HOLO STING rep 2"

corridas <- c("APO STING rep 1", "APO STING rep 2",
              "HOLO STING rep 1", "HOLO STING rep 2")

asymmetry$E <- as.factor(asymmetry$E)
asymmetry$E <- ordered(asymmetry$E, levels=corridas)

asymmetry_plot <- ggplot(asymmetry, aes(x = time/1000, y = rmsd)) + 
  geom_line(aes(color = E), size = 0.8) + ylim(0, 0.55) + 
  scale_color_manual(values = c("mediumseagreen", "forestgreen",
                                "dodgerblue3", "blue1")) +
  ylab("RMSD (nm)") + a + xlab("Tiempo (ns)") 
asymmetry_plot

png("figuras_tesis/asymmetry.png", height = 70, width = 184, units = "mm", res = 300)

ggarrange(asymmetry_plot, nrow = 1, ncol = 1, legend = "bottom", common.legend = TRUE) 
dev.off()
