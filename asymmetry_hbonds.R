setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)

#### ABIERTA 1 ####
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
  scale_color_manual(values = c("palegreen", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING rep 1")

abierta1_ala161a_gln292a <- read.table("render_hbonds/abierta1_ala161a_gln292a.dat", quote="\"", comment.char="")
abierta1_ala161b_gln292b <- read.table("render_hbonds/abierta1_ala161b_gln292b.dat", quote="\"", comment.char="")
abierta1_arg298a_glu154b <- read.table("render_hbonds/abierta1_arg298a_glu154b.dat", quote="\"", comment.char="")
abierta1_arg298b_glu154a <- read.table("render_hbonds/abierta1_arg298b_glu154a.dat", quote="\"", comment.char="")
abierta1_val160a_gln292a <- read.table("render_hbonds/abierta1_val160a_gln292a.dat", quote="\"", comment.char="")
abierta1_val160b_gln292b <- read.table("render_hbonds/abierta1_val160b_gln292b.dat", quote="\"", comment.char="")

a11 <- ggplot(abierta1_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
a12 <- ggplot(abierta1_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
a13 <- ggplot(abierta1_arg298a_glu154b, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154B") + ylim(0, 2)
a14 <- ggplot(abierta1_arg298b_glu154a, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154A") + ylim(0, 2)
a15 <- ggplot(abierta1_val160a_gln292a, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
a16 <- ggplot(abierta1_val160b_gln292b, aes(V1, V2)) + geom_point(color = "palegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/abierta1_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_open1, a11, a12, a15, a16, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/abierta1_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_open1, a13, a14, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()

#### ABIERTA 2 ####
dist_open2 <- ggplot(dist_opens[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("limegreen", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING rep 2") 

abierta2_ala161a_gln292a <- read.table("render_hbonds/abierta2_ala161a_gln292a.dat", quote="\"", comment.char="")
abierta2_ala161b_gln292b <- read.table("render_hbonds/abierta2_ala161b_gln292b.dat", quote="\"", comment.char="")
abierta2_arg298a_glu154b <- read.table("render_hbonds/abierta2_arg298a_glu154b.dat", quote="\"", comment.char="")
abierta2_arg298b_glu154a <- read.table("render_hbonds/abierta2_arg298b_glu154a.dat", quote="\"", comment.char="")
abierta2_val160a_gln292a <- read.table("render_hbonds/abierta2_val160a_gln292a.dat", quote="\"", comment.char="")
abierta2_val160b_gln292b <- read.table("render_hbonds/abierta2_val160b_gln292b.dat", quote="\"", comment.char="")

a21 <- ggplot(abierta2_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
a22 <- ggplot(abierta2_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
a23 <- ggplot(abierta2_arg298a_glu154b, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154B") + ylim(0, 2)
a24 <- ggplot(abierta2_arg298b_glu154a, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154A") + ylim(0, 2)
a25 <- ggplot(abierta2_val160a_gln292a, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
a26 <- ggplot(abierta2_val160b_gln292b, aes(V1, V2)) + geom_point(color = "limegreen") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/abierta2_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_open2, a21, a22, a25, a26, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/abierta2_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_open2, a23, a24, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()


#### CERRADA 1 ####
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

dist_close1 <- ggplot(dist_closed[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("blue", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 1") 

cerrada1_ala161a_gln292a <- read.table("render_hbonds/cerrada1_ala161a_gln292a.dat", quote="\"", comment.char="")
cerrada1_ala161b_gln292b <- read.table("render_hbonds/cerrada1_ala161b_gln292b.dat", quote="\"", comment.char="")
cerrada1_arg298a_glu154a <- read.table("render_hbonds/cerrada1_arg298a_glu154a.dat", quote="\"", comment.char="")
cerrada1_arg298b_glu154b <- read.table("render_hbonds/cerrada1_arg298b_glu154b.dat", quote="\"", comment.char="")
cerrada1_val160a_gln292a <- read.table("render_hbonds/cerrada1_val160a_gln292a.dat", quote="\"", comment.char="")
cerrada1_val160b_gln292b <- read.table("render_hbonds/cerrada1_val160b_gln292b.dat", quote="\"", comment.char="")

c11 <- ggplot(cerrada1_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
c12 <- ggplot(cerrada1_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
c13 <- ggplot(cerrada1_arg298a_glu154a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154A") + ylim(0, 2)
c14 <- ggplot(cerrada1_arg298b_glu154b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154B") + ylim(0, 2)
c15 <- ggplot(cerrada1_val160a_gln292a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
c16 <- ggplot(cerrada1_val160b_gln292b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/cerrada1_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_close1, c11, c12, c15, c16, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/cerrada1_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_close1, c13, c14, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()

#### CERRADA 2 ####
dist_close2 <- ggplot(dist_closed[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("dodgerblue", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Closed STING rep 2") 

cerrada2_ala161a_gln292a <- read.table("render_hbonds/cerrada2_ala161a_gln292a.dat", quote="\"", comment.char="")
cerrada2_ala161b_gln292b <- read.table("render_hbonds/cerrada2_ala161b_gln292b.dat", quote="\"", comment.char="")
cerrada2_arg298a_glu154a <- read.table("render_hbonds/cerrada2_arg298a_glu154a.dat", quote="\"", comment.char="")
cerrada2_arg298b_glu154b <- read.table("render_hbonds/cerrada2_arg298b_glu154b.dat", quote="\"", comment.char="")
cerrada2_val160a_gln292a <- read.table("render_hbonds/cerrada2_val160a_gln292a.dat", quote="\"", comment.char="")
cerrada2_val160b_gln292b <- read.table("render_hbonds/cerrada2_val160b_gln292b.dat", quote="\"", comment.char="")

c21 <- ggplot(cerrada2_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
c22 <- ggplot(cerrada2_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
c23 <- ggplot(cerrada2_arg298a_glu154a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154A") + ylim(0, 2)
c24 <- ggplot(cerrada2_arg298b_glu154b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154B") + ylim(0, 2)
c25 <- ggplot(cerrada2_val160a_gln292a, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
c26 <- ggplot(cerrada2_val160b_gln292b, aes(V1, V2)) + geom_point(color = "blue") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/cerrada2_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_close2, c21, c22, c25, c26, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/cerrada2_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_close2, c23, c24, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()

#### MUTANTE 1 ####
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

dist_v160m1 <- ggplot(dist_v160m[1:2002,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 1") 
 
mutante1_ala161a_gln292a <- read.table("render_hbonds/mutante1_ala161a_gln292a.dat", quote="\"", comment.char="")
mutante1_ala161b_gln292b <- read.table("render_hbonds/mutante1_ala161b_gln292b.dat", quote="\"", comment.char="")
mutante1_arg298a_glu154b <- read.table("render_hbonds/mutante1_arg298a_glu154b.dat", quote="\"", comment.char="")
mutante1_arg298b_glu154a <- read.table("render_hbonds/mutante1_arg298b_glu154a.dat", quote="\"", comment.char="")
mutante1_val160a_gln292a <- read.table("render_hbonds/mutante1_val160a_gln292a.dat", quote="\"", comment.char="")
mutante1_val160b_gln292b <- read.table("render_hbonds/mutante1_val160b_gln292b.dat", quote="\"", comment.char="")

m11 <- ggplot(mutante1_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
m12 <- ggplot(mutante1_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
m13 <- ggplot(mutante1_arg298a_glu154b, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154B") + ylim(0, 2)
m14 <- ggplot(mutante1_arg298b_glu154a, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154A") + ylim(0, 2)
m15 <- ggplot(mutante1_val160a_gln292a, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
m16 <- ggplot(mutante1_val160b_gln292b, aes(V1, V2)) + geom_point(color = "deeppink") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/mutante1_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m1, m11, m12, m15, m16, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/mutante1_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m1, m13, m14, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()
#### MUTANTE 2 ####
dist_v160m2 <- ggplot(dist_v160m[2003:4004,], aes(x = time, y = distance, color = label)) +
  geom_line(aes(y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("hotpink1", "darkorange1")) + a + ylim(-0.5, 0.5) +
  ggtitle("Normalized distance between the polymerization site and the center of mass of Open STING V160M rep 2")

mutante2_ala161a_gln292a <- read.table("render_hbonds/mutante2_ala161a_gln292a.dat", quote="\"", comment.char="")
mutante2_ala161b_gln292b <- read.table("render_hbonds/mutante2_ala161b_gln292b.dat", quote="\"", comment.char="")
mutante2_arg298a_glu154b <- read.table("render_hbonds/mutante2_arg298a_glu154b.dat", quote="\"", comment.char="")
mutante2_arg298b_glu154a <- read.table("render_hbonds/mutante2_arg298b_glu154a.dat", quote="\"", comment.char="")
mutante2_val160a_gln292a <- read.table("render_hbonds/mutante2_val160a_gln292a.dat", quote="\"", comment.char="")
mutante2_val160b_gln292b <- read.table("render_hbonds/mutante2_val160b_gln292b.dat", quote="\"", comment.char="")

m21 <- ggplot(mutante2_ala161a_gln292a, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161A :: GLN292A") + ylim(0, 2)
m22 <- ggplot(mutante2_ala161b_gln292b, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ALA161B :: GLN292B") + ylim(0, 2)
m23 <- ggplot(mutante2_arg298a_glu154b, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298A :: GLU154B") + ylim(0, 2)
m24 <- ggplot(mutante2_arg298b_glu154a, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("ARG298B :: GLU154A") + ylim(0, 2)
m25 <- ggplot(mutante2_val160a_gln292a, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160A :: GLN292A") + ylim(0, 2)
m26 <- ggplot(mutante2_val160b_gln292b, aes(V1, V2)) + geom_point(color = "hotpink1") + a + 
  ylab("Ocurrencia") + xlab("frame") + ggtitle("VAL160B :: GLN292B") + ylim(0, 2)

png("render_hbonds/mutante2_bonds1.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m2, m21, m22, m25, m26, ncol = 1, nrow = 5, common.legend = TRUE, legend = "bottom")
dev.off()

png("render_hbonds/mutante2_bonds2.png", height = 500, width = 750, units = "mm", res = 300)
ggarrange(dist_v160m2, m23, m24, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()

#### Asimetria ####
open1 <- read.delim("abierta-rep1/crop-200/asymmetry/abierta_rep1_rmsd.txt", sep = "=", header = FALSE)
open2 <- read.delim("abierta-rep2/asymmetry/abierta_rep2_rmsd.txt", sep = "=", header = FALSE)
close1 <- read.delim("cerrada-rep1/prod/asymmetry/cerrada_rep1_rmsd.txt", sep = "=", header = FALSE)
close2 <- read.delim("cerrada-rep2/asymmetry/cerrada_rep2_rmsd.txt", sep = "=", header = FALSE)
mut1 <- read.delim("v160m-rep1/prod_200ns/asymmetry/v160m_rep1_rmsd.txt", sep = "=", header = FALSE)
mut2 <- read.delim("v160m-rep2/asymmetry/v160m_rep2_rmsd.txt", sep = "=", header = FALSE)

open1[] <- lapply(open1, gsub, pattern = "nm", replacement = "", fixed = TRUE)
open1[] <- lapply(open1, gsub, pattern = " ", replacement = "", fixed = TRUE)
open2[] <- lapply(open2, gsub, pattern = "nm", replacement = "", fixed = TRUE)
open2[] <- lapply(open2, gsub, pattern = " ", replacement = "", fixed = TRUE)
close1[] <- lapply(close1, gsub, pattern = "nm", replacement = "", fixed = TRUE)
close1[] <- lapply(close1, gsub, pattern = " ", replacement = "", fixed = TRUE)
close2[] <- lapply(close2, gsub, pattern = "nm", replacement = "", fixed = TRUE)
close2[] <- lapply(close2, gsub, pattern = " ", replacement = "", fixed = TRUE)
mut1[] <- lapply(mut1, gsub, pattern = "nm", replacement = "", fixed = TRUE)
mut1[] <- lapply(mut1, gsub, pattern = " ", replacement = "", fixed = TRUE)
mut2[] <- lapply(mut2, gsub, pattern = "nm", replacement = "", fixed = TRUE)
mut2[] <- lapply(mut2, gsub, pattern = " ", replacement = "", fixed = TRUE)

rmsd <- data.frame("Time" = open1_chainA$time/1000,
                   "Open1" = as.numeric(open1$V2),
                   "Open2" = as.numeric(open2$V2),
                   "Close1" = as.numeric(close1$V2),
                   "Close2" = as.numeric(close2$V2), 
                   "Mut1" = as.numeric(mut1$V2),
                   "Mut2" = as.numeric(mut2$V2))

ggplot(rmsd, aes(x = Time, y = Open1)) + ylim(0, 0.7) + 
  geom_line(color = "palegreen", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Open STING rep 1")
ggplot(rmsd, aes(x = Time, y = Open2)) + ylim(0, 0.7) + 
  geom_line(color = "limegreen", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Open STING rep 2")
ggplot(rmsd, aes(x = Time, y = Close1)) + ylim(0, 0.7) + 
  geom_line(color = "blue", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Closed STING rep 1")
ggplot(rmsd, aes(x = Time, y = Close2)) + ylim(0, 0.7) + 
  geom_line(color = "dodgerblue", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("Closed STING rep 2")
ggplot(rmsd, aes(x = Time, y = Mut1)) + ylim(0, 0.7) + 
  geom_line(color = "deeppink", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("V160M STING rep 1")
ggplot(rmsd, aes(x = Time, y = Mut2)) + ylim(0, 0.7) + 
  geom_line(color = "hotpink1", size = 1) + 
  xlab("Time (ns)") + ylab("RMSD (nm)") + a + ggtitle("V160M STING rep 2")

rmsd[which.max(rmsd$Open1),]
rmsd[which.max(rmsd$Open2),]
rmsd[which.max(rmsd$Close1),]
rmsd[which.max(rmsd$Close2),]
rmsd[which.max(rmsd$Mut1),]
rmsd[which.max(rmsd$Mut2),]
