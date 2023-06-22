setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 30),
  axis.title = element_text(colour = 'black', face = 'italic', size = 15),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

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


dist_closed <- data.frame("time" = c(close1$time, close1$time, close2$time, close2$time),
                          "distance" = c(close1$chainA_norm, close1$chainB_norm, 
                                         close2$chainA_norm, close2$chainB_norm),
                          "label" = "label")
dist_closed$label[1:1001] <- "Closed STING 2'3'cGAMP rep 1 - Chain A"
dist_closed$label[1002:2002] <- "Closed STING 2'3'cGAMP rep 1 - Chain B"
dist_closed$label[2003:3003] <- "Closed STING 2'3'cGAMP rep 2 - Chain A"
dist_closed$label[3004:4004] <- "Closed STING 2'3'cGAMP rep 2 - Chain B"



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

dist_v160m <- data.frame("time" = c(v160m1$time, v160m1$time),
                         "distance" = c(v160m1$chainA_norm, v160m1$chainB_norm,
                                        v160m2$chainA_norm, v160m2$chainB_norm),
                         "label" = "label")
dist_v160m$label[1:1001] <- "Open STING V160M rep 1 - Chain A"
dist_v160m$label[1002:2002] <- "Open STING V160M rep 1 - Chain B"
dist_v160m$label[2003:3003] <- "Open STING V160M rep 2 - Chain A"
dist_v160m$label[3004:4004] <- "Open STING V160M rep 2 - Chain B"

boxes_dist <- data.frame("time" = c(open1$time[501:1001], open1$time[501:1001], open2$time[501:1001], open2$time[501:1001],
                                    close1$time[501:1001], close1$time[501:1001], close2$time[501:1001], close2$time[501:1001],
                                    v160m1$time[501:1001], v160m1$time[501:1001], v160m2$time[501:1001], v160m2$time[501:1001]),
                         "distance" = c(open1$chainA_norm[501:1001], open1$chainB_norm[501:1001], open2$chainA_norm[501:1001], open2$chainB_norm[501:1001],
                                       close1$chainA_norm[501:1001], close1$chainB_norm[501:1001], close2$chainA_norm[501:1001], close2$chainB_norm[501:1001],
                                       v160m1$chainA_norm[501:1001], v160m1$chainB_norm[501:1001], v160m2$chainA_norm[501:1001], v160m2$chainB_norm[501:1001]),
                         "Structure" = "label")
boxes_dist$Structure[1:501] <- "Open STING chain A rep 1"
boxes_dist$Structure[502:1002] <- "Open STING chain B rep 1"
boxes_dist$Structure[1003:1503] <- "Open STING chain A rep 2"
boxes_dist$Structure[1504:2004] <- "Open STING chain B rep 2"
boxes_dist$Structure[2005:2505] <- "Closed STING + 23cGAMP STING chain A rep 1"
boxes_dist$Structure[2506:3006] <- "Closed STING + 23cGAMP STING chain B rep 1"
boxes_dist$Structure[3007:3507] <- "Closed STING + 23cGAMP STING chain A rep 2"
boxes_dist$Structure[3508:4008] <- "Closed STING + 23cGAMP STING chain B rep 2"
boxes_dist$Structure[4009:4509] <- "Open STING V160M chain A rep 1"
boxes_dist$Structure[4510:5010] <- "Open STING V160M chain B rep 1"
boxes_dist$Structure[5011:5511] <- "Open STING V160M chain A rep 2"
boxes_dist$Structure[5512:6012] <- "Open STING V160M chain B rep 2"

orden <- c("Open STING chain A rep 1", "Open STING chain B rep 1", 
           "Open STING chain A rep 2", "Open STING chain B rep 2",
           "Closed STING + 23cGAMP STING chain A rep 1", "Closed STING + 23cGAMP STING chain B rep 1",
           "Closed STING + 23cGAMP STING chain A rep 2", "Closed STING + 23cGAMP STING chain B rep 2",
           "Open STING V160M chain A rep 1","Open STING V160M chain B rep 1",
           "Open STING V160M chain A rep 2","Open STING V160M chain B rep 2")

boxes_dist$Structure <- as.factor(boxes_dist$Structure)
boxes_dist$Structure <- ordered(boxes_dist$Structure, levels=orden)

png("pics-19Sep/dist_100ns_boxplot.png", height = 300, width = 450, units = "mm", res = 300)
ggplot(boxes_dist, aes(x = Structure, y = distance, color = Structure)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen", "limegreen", "black", "slategray4", 
                                "blue", "dodgerblue", "darkorchid2", "magenta4",
                                "red", "firebrick1", "deeppink", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=4.5, ymin= -0.5, ymax=-0.4, alpha=.2, fill="palegreen") +
  annotate("rect", xmin = 4.6, xmax=8.5, ymin= -0.5, ymax=-0.4, alpha=.2, fill="blue") +
  annotate("rect", xmin = 8.6, xmax=12.5, ymin= -0.5, ymax=-0.4, alpha=.2, fill="deeppink") +
  annotate("text", x=2, y=-0.45, label="Open STING", fontface="bold") + 
  annotate("text", x=6, y=-0.45, label="Closed STING + 23cGAMP",fontface="bold") +
  annotate("text", x=10, y=-0.45, label="Open STING V160M", fontface="bold") +
  annotate("text", x = 1, y = 0.125, label=as.character(round(sd(open1$chainA_dist)/mean(open1$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 2, y = 0.25, label=as.character(round(sd(open1$chainB_dist)/mean(open1$chainB_dist)*100)), fontface = "bold") +
  annotate("text", x = 3, y = 0, label=as.character(round(sd(open2$chainA_dist)/mean(open2$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 4, y = 0.25, label=as.character(round(sd(open2$chainB_dist)/mean(open2$chainB_dist)*100)), fontface = "bold") +
  annotate("text", x = 5, y = 0.5, label=as.character(round(sd(close1$chainA_dist)/mean(close1$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 6, y = 0.125, label=as.character(round(sd(close1$chainB_dist)/mean(close1$chainB_dist)*100)), fontface = "bold") +
  annotate("text", x = 7, y = 0.45, label=as.character(round(sd(close2$chainA_dist)/mean(close2$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 8, y = 0.125, label=as.character(round(sd(close2$chainB_dist)/mean(close2$chainB_dist)*100)), fontface = "bold") +
  annotate("text", x = 9, y = 0.35, label=as.character(round(sd(v160m1$chainA_dist)/mean(v160m1$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 10, y = 0.25, label=as.character(round(sd(v160m1$chainB_dist)/mean(v160m1$chainB_dist)*100)), fontface = "bold") +
  annotate("text", x = 11, y = 0.075, label=as.character(round(sd(v160m2$chainA_dist)/mean(v160m2$chainA_dist)*100)), fontface = "bold") + 
  annotate("text", x = 12, y = 0.2, label=as.character(round(sd(v160m2$chainB_dist)/mean(v160m2$chainB_dist)*100)), fontface = "bold") +
  xlab("Corridas") + a  + ylab("Distance (nm)") 
dev.off()


cor(boxes$sasas[1:501], boxes_dist$distance[1:501], method = "spearman") #0.8008917
cor(boxes$sasas[502:1002], boxes_dist$distance[502:1002], method = "spearman") #0.1923806
cor(boxes$sasas[1003:1503], boxes_dist$distance[1003:1503], method = "spearman") #0.2355776
cor(boxes$sasas[1504:2004], boxes_dist$distance[1504:2004], method = "spearman") #0.1683131
cor(boxes$sasas[2005:2505], boxes_dist$distance[2005:2505], method = "spearman") #0.4889512
cor(boxes$sasas[2506:3006], boxes_dist$distance[2506:3006], method = "spearman") #0.1997774
cor(boxes$sasas[3007:3507], boxes_dist$distance[3007:3507], method = "spearman") #-0.2120928
cor(boxes$sasas[3508:4008], boxes_dist$distance[3508:4008], method = "spearman") #0.4628458
cor(boxes$sasas[4009:4509], boxes_dist$distance[4009:4509], method = "spearman") #0.3422186
cor(boxes$sasas[4510:5010], boxes_dist$distance[4510:5010], method = "spearman") #0.2601577
cor(boxes$sasas[5011:5511], boxes_dist$distance[5011:5511], method = "spearman") #0.6108167
cor(boxes$sasas[5512:6012], boxes_dist$distance[5512:6012], method = "spearman") #0.3708378
