##### box plots ####

boxes <- data.frame("time" = c(abierta_rep1_sasa$Time[501:1001], abierta_rep1_sasa$Time[501:1001],
                               abierta_rep2_sasa$Time[501:1001], abierta_rep2_sasa$Time[501:1001],
                               cerrada_rep1_sasa$Time[501:1001], cerrada_rep1_sasa$Time[501:1001], 
                               cerrada_rep2_sasa$Time[501:1001], cerrada_rep2_sasa$Time[501:1001], 
                               v160m_rep1_sasa$Time[501:1001], v160m_rep1_sasa$Time[501:1001], 
                               v160m_rep2_sasa$Time[501:1001], v160m_rep2_sasa$Time[501:1001]),
                    "sasas" = c(abierta_rep1_sasa$polsiteA[501:1001], abierta_rep1_sasa$polsite[501:1001],
                                abierta_rep2_sasa$polsiteA[501:1001], abierta_rep2_sasa$polsiteB[501:1001], 
                                cerrada_rep1_sasa$polsiteA[501:1001], cerrada_rep1_sasa$polsiteB[501:1001], 
                                cerrada_rep2_sasa$polsiteA[501:1001], cerrada_rep2_sasa$polsiteB[501:1001],
                                v160m_rep1_sasa$polsiteA[501:1001], v160m_rep1_sasa$polsiteB[501:1001],
                                v160m_rep2_sasa$polsiteA[501:1001], v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484[501:1001]),
                    "label" = "label")

boxes$label[1:501] <- "Open STING chain A rep 1"
boxes$label[502:1002] <- "Open STING chain B rep 1"
boxes$label[1003:1503] <- "Open STING chain A rep 2"
boxes$label[1504:2004] <- "Open STING chain B rep 2"
boxes$label[2005:2505] <- "Closed STING + 23cGAMP STING chain A rep 1"
boxes$label[2506:3006] <- "Closed STING + 23cGAMP STING chain B rep 1"
boxes$label[3007:3507] <- "Closed STING + 23cGAMP STING chain A rep 2"
boxes$label[3508:4008] <- "Closed STING + 23cGAMP STING chain B rep 2"
boxes$label[4009:4509] <- "Open STING V160M chain A rep 1"
boxes$label[4510:5010] <- "Open STING V160M chain B rep 1"
boxes$label[5011:5511] <- "Open STING V160M chain A rep 2"
boxes$label[5512:6012] <- "Open STING V160M chain B rep 2"

orden <- c("Open STING chain A rep 1", "Open STING chain B rep 1", 
           "Open STING chain A rep 2", "Open STING chain B rep 2",
           "Closed STING + 23cGAMP STING chain A rep 1", "Closed STING + 23cGAMP STING chain B rep 1",
           "Closed STING + 23cGAMP STING chain A rep 2", "Closed STING + 23cGAMP STING chain B rep 2",
           "Open STING V160M chain A rep 1","Open STING V160M chain B rep 1",
           "Open STING V160M chain A rep 2","Open STING V160M chain B rep 2")

boxes$label <- as.factor(boxes$label)
boxes$label <- ordered(boxes$label, levels=orden)

png("pics-19Sep/sasas_100ns_boxplot1.png", height = 300, width = 450, units = "mm", res = 300)
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

boxes2 <- data.frame("time" = c(abierta_rep1_sasa$Time[501:1001], abierta_rep1_sasa$Time[501:1001],
                               abierta_rep2_sasa$Time[501:1001], abierta_rep2_sasa$Time[501:1001],
                               cerrada_rep1_sasa$Time[501:1001], cerrada_rep1_sasa$Time[501:1001], 
                               cerrada_rep2_sasa$Time[501:1001], cerrada_rep2_sasa$Time[501:1001], 
                               v160m_rep1_sasa$Time[501:1001], v160m_rep1_sasa$Time[501:1001], 
                               v160m_rep2_sasa$Time[501:1001], v160m_rep2_sasa$Time[501:1001]),
                    "sasas" = c(abierta_rep1_sasa$polsiteA[501:1001], abierta_rep1_sasa$polsite[501:1001],
                                abierta_rep2_sasa$polsiteA[501:1001], abierta_rep2_sasa$polsiteB[501:1001], 
                                cerrada_rep1_sasa$polsiteA[501:1001], cerrada_rep1_sasa$polsiteB[501:1001], 
                                cerrada_rep2_sasa$polsiteA[501:1001], cerrada_rep2_sasa$polsiteB[501:1001],
                                v160m_rep1_sasa$polsiteA[501:1001], v160m_rep1_sasa$polsiteB[501:1001],
                                v160m_rep2_sasa$polsiteA[501:1001], v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484[501:1001]),
                    "label" = "label")

boxes2$label[c(1:501, 1003:1503)] <- "Open STING chain A"
boxes2$label[c(502:1002, 1504:2004)] <- "Open STING chain B"
boxes2$label[c(2005:2505, 3007:3507)] <- "Closed STING + 23cGAMP chain A"
boxes2$label[c(2506:3006, 3508:4008)] <- "Closed STING + 23cGAMP chain B"
boxes2$label[c(4009:4509, 5011:5511)] <- "Open STING V160M chain A"
boxes2$label[c(4510:5010, 5512:6012)] <- "Open STING V160M chain B"

orden2 <- c("Open STING chain A", "Open STING chain B",
            "Closed STING + 23cGAMP chain A", "Closed STING + 23cGAMP chain B",
            "Open STING V160M chain A", "Open STING V160M chain B")

boxes2$label <- as.factor(boxes2$label)
boxes2$label <- ordered(boxes2$label, levels=orden2)
png("pics-19Sep/sasas_100ns_boxplot2.png", height = 300, width = 450, units = "mm", res = 300)
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

boxes3 <- data.frame("time" = c(abierta_rep1_sasa$Time[501:1001], abierta_rep1_sasa$Time[501:1001],
                               abierta_rep2_sasa$Time[501:1001], abierta_rep2_sasa$Time[501:1001],
                               cerrada_rep1_sasa$Time[501:1001], cerrada_rep1_sasa$Time[501:1001], 
                               cerrada_rep2_sasa$Time[501:1001], cerrada_rep2_sasa$Time[501:1001], 
                               v160m_rep1_sasa$Time[501:1001], v160m_rep1_sasa$Time[501:1001], 
                               v160m_rep2_sasa$Time[501:1001], v160m_rep2_sasa$Time[501:1001]),
                    "sasas" = c(abierta_rep1_sasa$polsiteA[501:1001], abierta_rep1_sasa$polsite[501:1001],
                                abierta_rep2_sasa$polsiteA[501:1001], abierta_rep2_sasa$polsiteB[501:1001], 
                                cerrada_rep1_sasa$polsiteA[501:1001], cerrada_rep1_sasa$polsiteB[501:1001], 
                                cerrada_rep2_sasa$polsiteA[501:1001], cerrada_rep2_sasa$polsiteB[501:1001],
                                v160m_rep1_sasa$polsiteA[501:1001], v160m_rep1_sasa$polsiteB[501:1001],
                                v160m_rep2_sasa$polsiteA[501:1001], v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484[501:1001]),
                    "label" = "label")

boxes3$label[1:1002] <- "Open STING rep 1"
boxes3$label[1003:2004] <- "Open STING rep 2"
boxes3$label[2005:3006] <- "Closed STING + 23cGAMP rep 1"
boxes3$label[3007:4008] <- "Closed STING + 23cGAMP rep 2"
boxes3$label[4009:5010] <- "Open STING V160M rep 1"
boxes3$label[5011:6012] <- "Open STING V160M rep 2"

orden3 <- c("Open STING rep 1", "Open STING rep 2",
            "Closed STING + 23cGAMP rep 1", "Closed STING + 23cGAMP rep 2",
            "Open STING V160M rep 1", "Open STING V160M rep 2")

boxes3$label <- as.factor(boxes3$label)
boxes3$label <- ordered(boxes3$label, levels=orden3)
png("pics-19Sep/sasas_100ns_boxplot3.png", height = 300, width = 450, units = "mm", res = 300)
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

boxes4 <- data.frame("time" = c(abierta_rep1_sasa$Time[501:1001], abierta_rep1_sasa$Time[501:1001],
                               abierta_rep2_sasa$Time[501:1001], abierta_rep2_sasa$Time[501:1001],
                               cerrada_rep1_sasa$Time[501:1001], cerrada_rep1_sasa$Time[501:1001], 
                               cerrada_rep2_sasa$Time[501:1001], cerrada_rep2_sasa$Time[501:1001], 
                               v160m_rep1_sasa$Time[501:1001], v160m_rep1_sasa$Time[501:1001], 
                               v160m_rep2_sasa$Time[501:1001], v160m_rep2_sasa$Time[501:1001]),
                    "sasas" = c(abierta_rep1_sasa$polsiteA[501:1001], abierta_rep1_sasa$polsite[501:1001],
                                abierta_rep2_sasa$polsiteA[501:1001], abierta_rep2_sasa$polsiteB[501:1001], 
                                cerrada_rep1_sasa$polsiteA[501:1001], cerrada_rep1_sasa$polsiteB[501:1001], 
                                cerrada_rep2_sasa$polsiteA[501:1001], cerrada_rep2_sasa$polsiteB[501:1001],
                                v160m_rep1_sasa$polsiteA[501:1001], v160m_rep1_sasa$polsiteB[501:1001],
                                v160m_rep2_sasa$polsiteA[501:1001], v160m_rep2_sasa$r_475_r_476_r_477_r_478_r_479_r_480_r_481_r_482_r_483_r_484[501:1001]),
                    "label" = "label")

boxes4$label[1:2004] <- "Open STING"
boxes4$label[2005:4008] <- "Closed STING + 23cGAMP"
boxes4$label[4009:6012] <- "Open STING V160M"

orden4 <- c("Open STING", "Closed STING + 23cGAMP", "Open STING V160M")

boxes4$label <- as.factor(boxes4$label)
boxes4$label <- ordered(boxes4$label, levels=orden4)
png("pics-19Sep/sasas_100ns_boxplot4.png", height = 300, width = 450, units = "mm", res = 300)
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
