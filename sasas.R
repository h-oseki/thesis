abierta_rep1_sasa <- readXVG("abierta-rep1/crop-200/sasa/abierta-rep1-sasa.xvg")
abierta_rep1_sasa <- as.data.frame(sapply(abierta_rep1_sasa, as.numeric))
abierta_rep2_sasa <- readXVG("abierta-rep2/sasa/abierta-rep2-sasa.xvg")
abierta_rep2_sasa <- as.data.frame(sapply(abierta_rep2_sasa, as.numeric))
cerrada_rep1_sasa <- readXVG("cerrada-rep1/prod/sasa/cerrada-rep1-sasa.xvg")
cerrada_rep1_sasa <- as.data.frame(sapply(cerrada_rep1_sasa, as.numeric))
cerrada_rep2_sasa <- readXVG("cerrada-rep2/sasa/cerrada-rep2-sasa.xvg")
cerrada_rep2_sasa <- as.data.frame(sapply(cerrada_rep2_sasa, as.numeric))
v160m_rep1_sasa <- readXVG("v160m-rep1/prod_200ns/sasa/v160m-rep1-sasa.xvg")
v160m_rep1_sasa <- as.data.frame(sapply(v160m_rep1_sasa, as.numeric))
v160m_rep1_sasa$Time <- v160m_rep1_sasa$Time/1000
v160m_rep2_sasa <- readXVG("v160m-rep2/sasa/v160m-rep2-sasa.xvg")
v160m_rep2_sasa <- as.data.frame(sapply(v160m_rep2_sasa, as.numeric))
rollmean(abierta_rep1_sasa$polsiteA, 250)

sasas <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                               abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                              cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time,
                              cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time,
                              v160m_rep1_sasa$Time, v160m_rep1_sasa$Time,
                              v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                    "sasa" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                               abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB,
                               cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB,
                               cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB,
                               v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB,
                               v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                    "label" = 'label')
sasas$label[1:1001] <- "Open STING rep 1 - Chain A"
sasas$label[1002:2002] <- "Open STING rep 1 - Chain B"
sasas$label[2003:3003] <- "Open STING rep 2 - Chain A"
sasas$label[3004:4004] <- "Open STING rep 2 - Chain B"
sasas$label[4005:5005] <- "Closed STING 2'3'cGAMP rep 1 - Chain A"
sasas$label[5006:6006] <- "Closed STING 2'3'cGAMP rep 1 - Chain B"
sasas$label[6007:7007] <- "Closed STING 2'3'cGAMP rep 2 - Chain A"
sasas$label[7008:8008] <- "Closed STING 2'3'cGAMP rep 2 - Chain B"
sasas$label[8009:9009] <- "Open STING V160M rep 1 - Chain A"
sasas$label[9010:10010] <- "Open STING V160M rep 1 - Chain B"
sasas$label[10011:11011] <- "Open STING V160M rep 2 - Chain A"
sasas$label[11012:12012] <- "Open STING V160M rep 2 - Chain B"

sasa_open1 <- ggplot(sasas[1:2002,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("SASA (nm2)")  + 
  scale_color_manual(values = c("palegreen", "black")) + a  +
  ggtitle("SASA of Open STING rep 1") +
  geom_ysidedensity()

sasa_open2 <- ggplot(sasas[2003:4004,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("SASA (nm2)")  + 
  scale_color_manual(values = c("limegreen", "black")) + a  +
  ggtitle("SASA of Open STING rep 2") +
  geom_ysidedensity()

sasa_close1 <- ggplot(sasas[4005:6006,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("SASA (nm2)")  + 
  scale_color_manual(values = c("blue", "black")) + a  +
  ggtitle("SASA of Closed STING rep 1") +
  geom_ysidedensity()
sasa_close2 <- ggplot(sasas[6007:8008,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("SASA (nm2)")  + 
  scale_color_manual(values = c("dodgerblue", "black")) + a  +
  ggtitle("SASA Closed STING rep 2") +
  geom_ysidedensity()

sasa_v160m1 <- ggplot(sasas[8009:10010,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "black")) + a + 
  ggtitle("SASA of Open STING V160M rep 1") +
  geom_ysidedensity()
sasa_v160m2 <- ggplot(sasas[10011:12012,], aes(x = time, y = sasa, color = label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("hotpink1", "black")) + a + 
  ggtitle("SASA of Open STING V160M rep 2") +
  geom_ysidedensity()

png("pics-19Sep/abiertas_sasa.png", height = 300, width = 250, units = "mm", res = 300)
ggplot() + 
  geom_line(data = abierta_rep1_sasa, aes(Time, polsiteA), color = "palegreen", linetype = "solid", size = 1) +
  geom_line(data = abierta_rep1_sasa, aes(Time, polsite), color = "palegreen", linetype = "dashed", size = 1) +
  geom_line(data = abierta_rep2_sasa, aes(Time, polsiteA), color = "limegreen", linetype = "solid", size = 1) +
  geom_line(data = abierta_rep2_sasa, aes(Time, polsiteB), color = "limegreen", linetype = "dashed", size = 1) + a
dev.off()

png("pics-19Sep/cerradas_sasa.png", height = 300, width = 250, units = "mm", res = 300)
ggplot() + 
  geom_line(data = cerrada_rep1_sasa, aes(Time, polsiteA), color = "blue", linetype = "solid", size = 1) +
  geom_line(data = cerrada_rep1_sasa, aes(Time, polsiteB), color = "blue", linetype = "dashed", size = 1) +
  geom_line(data = cerrada_rep2_sasa, aes(Time, polsiteA), color = "dodgerblue", linetype = "solid", size = 1) +
  geom_line(data = cerrada_rep2_sasa, aes(Time, polsiteB), color = "dodgerblue", linetype = "dashed", size = 1) + a
dev.off()

png("pics-19Sep/v160m_sasa.png", height = 300, width = 250, units = "mm", res = 300)
ggplot() + 
  geom_line(data = v160m_rep1_sasa, aes(Time, polsiteA), color = "deeppink", linetype = "solid", size = 1) +
  geom_line(data = v160m_rep1_sasa, aes(Time, polsiteB), color = "deeppink", linetype = "dashed", size = 1) +
  geom_line(data = v160m_rep2_sasa, aes(Time, polsiteA), color = "hotpink1", linetype = "solid", size = 1) +
  geom_line(data = v160m_rep2_sasa, aes(Time, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484), color = "hotpink1", linetype = "dashed", size = 1) + a
dev.off()


##### roll mean #####

sasas <- data.frame("time" = rollmean(abierta_rep1_sasa$Time, 250),
                    "abierta_rep1_A" = rollmean(abierta_rep1_sasa$polsiteA, 250),
                    "abierta_rep1_B" = rollmean(abierta_rep1_sasa$polsite, 250),
                    "abierta_rep2_A" = rollmean(abierta_rep2_sasa$polsiteA, 250),
                    "abierta_rep2_B" = rollmean(abierta_rep2_sasa$polsiteB, 250),
                    "cerrada_rep1_A" = rollmean(cerrada_rep1_sasa$polsiteA, 250),
                    "cerrada_rep1_B" = rollmean(cerrada_rep1_sasa$polsiteB, 250),
                    "cerrada_rep2_A" = rollmean(cerrada_rep2_sasa$polsiteA, 250),
                    "cerrada_rep2_B" = rollmean(cerrada_rep2_sasa$polsiteB, 250),
                    "v160m_rep1_A" = rollmean(v160m_rep1_sasa$polsiteA, 250),
                    "v160m_rep1_B" = rollmean(v160m_rep1_sasa$polsiteB, 250),
                    "v160m_rep2_A" = rollmean(v160m_rep2_sasa$polsiteA, 250),
                    "v160m_rep2_B" = rollmean(v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484, 250))

png("pics-19Sep/SASAs.png", height = 300, width = 250, units = "mm", res = 300)
ggplot(sasas, aes(x = time)) +  xlab("Time (ns)") + ylab("SASA (nm2)") + a +
  geom_line(aes(y = abierta_rep1_A), color = "palegreen", linetype = "solid", size = 1) +
  geom_line(aes(y = abierta_rep1_B), color = "palegreen", linetype = "dashed", size = 1) +
  geom_line(aes(y = abierta_rep2_A), color = "limegreen", linetype = "solid", size = 1) +
  geom_line(aes(y = abierta_rep2_B), color = "limegreen", linetype = "dashed", size = 1) +
  geom_line(aes(y = cerrada_rep1_A), color = "blue", linetype = "solid", size = 1) +
  geom_line(aes(y = cerrada_rep1_B), color = "blue", linetype = "dashed", size = 1) +
  geom_line(aes(y = cerrada_rep2_A), color = "dodgerblue", linetype = "solid", size = 1) +
  geom_line(aes(y = cerrada_rep2_B), color = "dodgerblue", linetype = "dashed", size = 1) + 
  geom_line(aes(y = v160m_rep1_A), color = "deeppink", linetype = "solid", size = 1) +
  geom_line(aes(y = v160m_rep1_B), color = "deeppink", linetype = "dashed", size = 1) +
  geom_line(aes(y = v160m_rep2_A), color = "hotpink1", linetype = "solid", size = 1) +
  geom_line(aes(y = v160m_rep2_B), color = "hotpink1", linetype = "dashed", size = 1) +
  annotate("text", x=-Inf, y=9.7, label="Open Rep 1", hjust=-.2, vjust=2, fontface="bold", color = "palegreen") + 
  annotate("text", x=-Inf, y=9.6, label="Open Rep 2", hjust=-.2, vjust=2, fontface="bold", color = "limegreen") + 
  annotate("text", x=-Inf, y=9.5, label="Closed Rep 1", hjust=-.2, vjust=2, fontface="bold", color = "blue") + 
  annotate("text", x=-Inf, y=9.4, label="Closed Rep 2", hjust=-.2, vjust=2, fontface="bold", color = "dodgerblue") +
  annotate("text", x=-Inf, y=9.3, label="V160M Rep 1", hjust=-.2, vjust=2, fontface="bold", color = "deeppink") + 
  annotate("text", x=-Inf, y=9.2, label="V160M Rep 2", hjust=-.2, vjust=2, fontface="bold", color = "hotpink1")
dev.off()


ggplot(sasas, aes(x = time)) +  xlab("Time (ns)") + ylab("SASA (nm2)") + a +
  geom_line(aes(y = v160m_rep1_A), color = "deeppink", linetype = "solid", size = 1) +
  geom_line(aes(y = v160m_rep1_B), color = "deeppink", linetype = "dashed", size = 1)# +
ggplot(sasas, aes(x = time)) +  xlab("Time (ns)") + ylab("SASA (nm2)") + a +
  geom_line(aes(y = v160m_rep2_A), color = "hotpink1", linetype = "solid", size = 1) +
  geom_line(aes(y = v160m_rep2_B), color = "hotpink1", linetype = "dashed", size = 1) 
 
  
##### box plots ####

boxes <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                               abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                               cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time, 
                               cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time, 
                               v160m_rep1_sasa$Time, v160m_rep1_sasa$Time, 
                               v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                    "sasas" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                                abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB, 
                                cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB, 
                                cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB,
                                v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB,
                                v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                    "label" = "label")

boxes$label[1:1001] <- "Open STING chain A rep 1"
boxes$label[1002:2002] <- "Open STING chain B rep 1"
boxes$label[2003:3003] <- "Open STING chain A rep 2"
boxes$label[3004:4004] <- "Open STING chain B rep 2"
boxes$label[4005:5005] <- "Closed STING + 23cGAMP STING chain A rep 1"
boxes$label[5006:6006] <- "Closed STING + 23cGAMP STING chain B rep 1"
boxes$label[6007:7007] <- "Closed STING + 23cGAMP STING chain A rep 2"
boxes$label[7008:8008] <- "Closed STING + 23cGAMP STING chain B rep 2"
boxes$label[8009:9009] <- "Open STING V160M chain A rep 1"
boxes$label[9010:10010] <- "Open STING V160M chain B rep 1"
boxes$label[10011:11011] <- "Open STING V160M chain A rep 2"
boxes$label[11012:12012] <- "Open STING V160M chain B rep 2"

orden <- c("Open STING chain A rep 1", "Open STING chain B rep 1", 
           "Open STING chain A rep 2", "Open STING chain B rep 2",
           "Closed STING + 23cGAMP STING chain A rep 1", "Closed STING + 23cGAMP STING chain B rep 1",
           "Closed STING + 23cGAMP STING chain A rep 2", "Closed STING + 23cGAMP STING chain B rep 2",
           "Open STING V160M chain A rep 1","Open STING V160M chain B rep 1",
           "Open STING V160M chain A rep 2","Open STING V160M chain B rep 2")

boxes$label <- as.factor(boxes$label)
boxes$label <- ordered(boxes$label, levels=orden)

png("pics-19Sep/sasas_boxplot1.png", height = 300, width = 450, units = "mm", res = 300)
ggplot(boxes, aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen", "limegreen", "black", "slategray4", 
                                "blue", "dodgerblue", "darkorchid2", "magenta4",
                                "red", "firebrick1", "deeppink", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=4.5, ymin= 5, ymax=5.5, alpha=.2, fill="palegreen") +
  annotate("rect", xmin = 4.6, xmax=8.5, ymin= 5, ymax=5.5, alpha=.2, fill="blue") +
  annotate("rect", xmin = 8.6, xmax=12.5, ymin= 5, ymax=5.5, alpha=.2, fill="deeppink") +
  annotate("text", x=2, y=5.25, label="Open STING", fontface="bold") + 
  annotate("text", x=6, y=5.25, label="Closed STING + 23cGAMP",fontface="bold") +
  annotate("text", x=10, y=5.25, label="Open STING V160M", fontface="bold") +
  xlab("Corridas") + a  + ylab("SASA (nm2)") 
dev.off()

boxes2 <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                                abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                                cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time, 
                                cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time, 
                                v160m_rep1_sasa$Time, v160m_rep1_sasa$Time, 
                                v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                     "sasas" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                                 abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB, 
                                 cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB, 
                                 cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB,
                                 v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB,
                                 v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                     "label" = "label")

boxes2$label[c(1:1001, 2003:3003)] <- "Open STING chain A"
boxes2$label[c(1002:2002, 3004:4004)] <- "Open STING chain B"
boxes2$label[c(4005:5005, 6007:7007)] <- "Closed STING + 23cGAMP chain A"
boxes2$label[c(5006:6006, 7008:8008)] <- "Closed STING + 23cGAMP chain B"
boxes2$label[c(8009:9009, 10011:11011)] <- "Open STING V160M chain A"
boxes2$label[c(9010:10010, 11012:12012)] <- "Open STING V160M chain B"

orden2 <- c("Open STING chain A", "Open STING chain B",
           "Closed STING + 23cGAMP chain A", "Closed STING + 23cGAMP chain B",
           "Open STING V160M chain A", "Open STING V160M chain B")

boxes2$label <- as.factor(boxes2$label)
boxes2$label <- ordered(boxes2$label, levels=orden2)

ab <- ggplot(boxes2[1:4004,], aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen", "slategray4")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="palegreen") +
  annotate("text", x=1.5, y=5.25, label="Open STING", fontface="bold") + 
  xlab("Corridas") + a  + ylab("SASA (nm2)") 
cl <- ggplot(boxes2[4005:8008,], aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("blue",  "magenta4")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="blue") +
  annotate("text", x=1.5, y=5.25, label="Closed STING + 23cGAMP",fontface="bold") + 
  xlab("Corridas") + a  + ylab("SASA (nm2)")

ggplot(boxes2[8009:12012,], aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("red", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="hotpink1") +
  annotate("text", x=1.5, y=5.25, label="Open STING V160M",fontface="bold") + 
  xlab("Corridas") + a  + ylab("SASA (nm2)") 


ggarrange(ab, cl, ncol = 2, nrow = 1)

png("pics-19Sep/sasas_boxplot2.png", height = 300, width = 450, units = "mm", res = 300)
ggplot(boxes2, aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen", "slategray4", 
                                "blue",  "magenta4",
                                "red", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="palegreen") +
  annotate("rect", xmin = 2.6, xmax=4.5, ymin= 5, ymax=5.5, alpha=.2, fill="blue") +
  annotate("rect", xmin = 4.6, xmax=6.5, ymin= 5, ymax=5.5, alpha=.2, fill="hotpink1") +
  annotate("text", x=1, y=5.25, label="Open STING", fontface="bold") + 
  annotate("text", x=3, y=5.25, label="Closed STING + 23cGAMP",fontface="bold") +
  annotate("text", x=5, y=5.25, label="Open STING V160M", fontface="bold") +
  xlab("Corridas") + a  + ylab("SASA (nm2)") 
dev.off()

boxes3 <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                                abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                                cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time, 
                                cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time, 
                                v160m_rep1_sasa$Time, v160m_rep1_sasa$Time, 
                                v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                     "sasas" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                                 abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB, 
                                 cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB, 
                                 cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB,
                                 v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB,
                                 v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                     "label" = "label")

boxes3$label[1:2002] <- "Open STING rep 1"
boxes3$label[2003:4004] <- "Open STING rep 2"
boxes3$label[4005:6006] <- "Closed STING + 23cGAMP rep 1"
boxes3$label[6007:8008] <- "Closed STING + 23cGAMP rep 2"
boxes3$label[8009:10010] <- "Open STING V160M rep 1"
boxes3$label[10011:12012] <- "Open STING V160M rep 2"

orden3 <- c("Open STING rep 1", "Open STING rep 2",
            "Closed STING + 23cGAMP rep 1", "Closed STING + 23cGAMP rep 2",
            "Open STING V160M rep 1", "Open STING V160M rep 2")

boxes3$label <- as.factor(boxes3$label)
boxes3$label <- ordered(boxes3$label, levels=orden3)
png("pics-19Sep/sasas_boxplot3.png", height = 300, width = 450, units = "mm", res = 300)
ggplot(boxes3, aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen", "limegreen", 
                                "blue",  "dodgerblue",
                                "deeppink", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="palegreen") +
  annotate("rect", xmin = 2.6, xmax=4.5, ymin= 5, ymax=5.5, alpha=.2, fill="blue") +
  annotate("rect", xmin = 4.6, xmax=6.5, ymin= 5, ymax=5.5, alpha=.2, fill="hotpink1") +
  annotate("text", x=1, y=5.25, label="Open STING", fontface="bold") + 
  annotate("text", x=3, y=5.25, label="Closed STING + 23cGAMP",fontface="bold") +
  annotate("text", x=5, y=5.25, label="Open STING V160M", fontface="bold") +
  xlab("Corridas") + a  + ylab("SASA (nm2)") 
dev.off()
boxes4 <- data.frame("time" = c(abierta_rep1_sasa$Time, abierta_rep1_sasa$Time,
                                abierta_rep2_sasa$Time, abierta_rep2_sasa$Time,
                                cerrada_rep1_sasa$Time, cerrada_rep1_sasa$Time, 
                                cerrada_rep2_sasa$Time, cerrada_rep2_sasa$Time, 
                                v160m_rep1_sasa$Time, v160m_rep1_sasa$Time, 
                                v160m_rep2_sasa$Time, v160m_rep2_sasa$Time),
                     "sasas" = c(abierta_rep1_sasa$polsiteA, abierta_rep1_sasa$polsite,
                                 abierta_rep2_sasa$polsiteA, abierta_rep2_sasa$polsiteB, 
                                 cerrada_rep1_sasa$polsiteA, cerrada_rep1_sasa$polsiteB, 
                                 cerrada_rep2_sasa$polsiteA, cerrada_rep2_sasa$polsiteB,
                                 v160m_rep1_sasa$polsiteA, v160m_rep1_sasa$polsiteB,
                                 v160m_rep2_sasa$polsiteA, v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484),
                     "label" = "label")

boxes4$label[1:4004] <- "Open STING"
boxes4$label[4005:8008] <- "Closed STING + 23cGAMP"
boxes4$label[8009:12012] <- "Open STING V160M"

orden4 <- c("Open STING", "Closed STING + 23cGAMP", "Open STING V160M")

boxes4$label <- as.factor(boxes4$label)
boxes4$label <- ordered(boxes4$label, levels=orden4)
png("pics-19Sep/sasas_boxplot4.png", height = 300, width = 450, units = "mm", res = 300)
ggplot(boxes4, aes(x = label, y = sasas, color = label)) + geom_boxplot() +
  scale_color_manual(values = c("palegreen","blue", "hotpink1")) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  annotate("rect", xmin = 0.5, xmax=1.5, ymin= 5, ymax=5.5, alpha=.2, fill="palegreen") +
  annotate("rect", xmin = 1.6, xmax=2.5, ymin= 5, ymax=5.5, alpha=.2, fill="blue") +
  annotate("rect", xmin = 2.6, xmax=3.5, ymin= 5, ymax=5.5, alpha=.2, fill="hotpink1") +
  annotate("text", x=1, y=5.25, label="Open STING", fontface="bold") + 
  annotate("text", x=2, y=5.25, label="Closed STING + 23cGAMP",fontface="bold") +
  annotate("text", x=3, y=5.25, label="Open STING V160M", fontface="bold") +
  xlab("Corridas") + a  + ylab("SASA (nm2)") 
dev.off()
