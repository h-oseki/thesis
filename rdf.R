setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 30),
  axis.title = element_text(colour = 'black', face = 'italic', size = 15),
  panel.background = element_rect(fill = "white", colour = "black"),
  panel.border = element_rect(colour="black", fill=NA, size=.6)
)
#
#### Polimerization Site ####

abierta1_all_chainA <- readXVG("rdf/abierta1/abierta1_rdf_polsiteA.xvg")
abierta1_all_chainA <- as.data.frame(sapply(abierta1_all_chainA, as.numeric))
abierta1_all_chainB <- readXVG("rdf/abierta1/abierta1_rdf_polsiteB.xvg")
abierta1_all_chainB <- as.data.frame(sapply(abierta1_all_chainB, as.numeric))
abierta1_last_chainA <- readXVG("rdf/abierta1/abierta1_rdf_last100_polsiteA.xvg")
abierta1_last_chainA <- as.data.frame(sapply(abierta1_last_chainA, as.numeric))
abierta1_last_chainB <- readXVG("rdf/abierta1/abierta1_rdf_last100_polsiteB.xvg")
abierta1_last_chainB <- as.data.frame(sapply(abierta1_last_chainB, as.numeric))

abierta1 <- data.frame("r" = c(abierta1_all_chainA$r, abierta1_all_chainB$r, 
                               abierta1_last_chainA$r, abierta1_last_chainB$r), 
                       "g" = c(abierta1_all_chainA$`Water_&_OW`, abierta1_all_chainB$`Water_&_OW`,
                                  abierta1_last_chainA$`Water_&_OW`, abierta1_last_chainB$`Water_&_OW`),
                       "label" = "label")

abierta1$label[1:2629] <- "Chain A"
abierta1$label[2630:5258] <- "Chain B"
abierta1$label[5259:7879] <- "Chain A Last 100ns"
abierta1$label[7880:10500] <- "Chain B Last 100ns"

chains <- c("Chain A", "Chain B", "Chain A Last 100ns", "Chain B Last 100ns")

abierta1$label <- as.factor(abierta1$label)
abierta1$label <- ordered(abierta1$label, levels=chains)


abierta1_plot1 <- ggplot(abierta1[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 1")
abierta1_plot2 <- ggplot(abierta1[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 1")


abierta2_all_chainA <- readXVG("rdf/abierta2/abierta2_rdf_polsiteA.xvg")
abierta2_all_chainA <- as.data.frame(sapply(abierta2_all_chainA, as.numeric))
abierta2_all_chainB <- readXVG("rdf/abierta2/abierta2_rdf_polsiteB.xvg")
abierta2_all_chainB <- as.data.frame(sapply(abierta2_all_chainB, as.numeric))
abierta2_last_chainA <- readXVG("rdf/abierta2/abierta2_rdf_last100_polsiteA.xvg")
abierta2_last_chainA <- as.data.frame(sapply(abierta2_last_chainA, as.numeric))
abierta2_last_chainB <- readXVG("rdf/abierta2/abierta2_rdf_last100_polsiteB.xvg")
abierta2_last_chainB <- as.data.frame(sapply(abierta2_last_chainB, as.numeric))

abierta2 <- data.frame("r" = c(abierta2_all_chainA$r, abierta2_all_chainB$r, 
                               abierta2_last_chainA$r, abierta2_last_chainB$r), 
                       "g" = c(abierta2_all_chainA$`Water_&_OW`, abierta2_all_chainB$`Water_&_OW`,
                               abierta2_last_chainA$`Water_&_OW`, abierta2_last_chainB$`Water_&_OW`),
                       "label" = "label")

abierta2$label[1:2629] <- "Chain A"
abierta2$label[2630:5258] <- "Chain B"
abierta2$label[5259:7879] <- "Chain A Last 100ns"
abierta2$label[7880:10500] <- "Chain B Last 100ns"

abierta2$label <- as.factor(abierta2$label)
abierta2$label <- ordered(abierta2$label, levels=chains)

abierta2_plot1 <- ggplot(abierta2[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) +  xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 2")
abierta2_plot2 <- ggplot(abierta2[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) +  xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 2")

cerrada1_all_chainA <- readXVG("rdf/cerrada1/cerrada1_rdf_polsiteA.xvg")
cerrada1_all_chainA <- as.data.frame(sapply(cerrada1_all_chainA, as.numeric))
cerrada1_all_chainB <- readXVG("rdf/cerrada1/cerrada1_rdf_polsiteB.xvg")
cerrada1_all_chainB <- as.data.frame(sapply(cerrada1_all_chainB, as.numeric))
cerrada1_last_chainA <- readXVG("rdf/cerrada1/cerrada1_rdf_last100_polsiteA.xvg")
cerrada1_last_chainA <- as.data.frame(sapply(cerrada1_last_chainA, as.numeric))
cerrada1_last_chainB <- readXVG("rdf/cerrada1/cerrada1_rdf_last100_polsiteB.xvg")
cerrada1_last_chainB <- as.data.frame(sapply(cerrada1_last_chainB, as.numeric))

cerrada1 <- data.frame("r" = c(cerrada1_all_chainA$r, cerrada1_all_chainB$r, 
                               cerrada1_last_chainA$r, cerrada1_last_chainB$r), 
                       "g" = c(cerrada1_all_chainA$`Water_&_OW`, cerrada1_all_chainB$`Water_&_OW`,
                               cerrada1_last_chainA$`Water_&_OW`, cerrada1_last_chainB$`Water_&_OW`),
                       "label" = "label")

cerrada1$label[1:2479] <- "Chain A"
cerrada1$label[2480:4958] <- "Chain B"
cerrada1$label[4959:7433] <- "Chain A Last 100ns"
cerrada1$label[7434:9908] <- "Chain B Last 100ns"

cerrada1$label <- as.factor(cerrada1$label)
cerrada1$label <- ordered(cerrada1$label, levels=chains)

cerrada1_plot1 <- ggplot(cerrada1[1:4958,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 1")
cerrada1_plot2 <- ggplot(cerrada1[4959:9908,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 1")

cerrada2_all_chainA <- readXVG("rdf/cerrada2/cerrada2_rdf_polsiteA.xvg")
cerrada2_all_chainA <- as.data.frame(sapply(cerrada2_all_chainA, as.numeric))
cerrada2_all_chainB <- readXVG("rdf/cerrada2/cerrada2_rdf_polsiteB.xvg")
cerrada2_all_chainB <- as.data.frame(sapply(cerrada2_all_chainB, as.numeric))
cerrada2_last_chainA <- readXVG("rdf/cerrada2/cerrada2_rdf_last100_polsiteA.xvg")
cerrada2_last_chainA <- as.data.frame(sapply(cerrada2_last_chainA, as.numeric))
cerrada2_last_chainB <- readXVG("rdf/cerrada2/cerrada2_rdf_last100_polsiteB.xvg")
cerrada2_last_chainB <- as.data.frame(sapply(cerrada2_last_chainB, as.numeric))

cerrada2 <- data.frame("r" = c(cerrada2_all_chainA$r, cerrada2_all_chainB$r, 
                               cerrada2_last_chainA$r, cerrada2_last_chainB$r), 
                       "g" = c(cerrada2_all_chainA$`Water_&_OW`, cerrada2_all_chainB$`Water_&_OW`,
                               cerrada2_last_chainA$`Water_&_OW`, cerrada2_last_chainB$`Water_&_OW`),
                       "label" = "label")

cerrada2$label[1:2479] <- "Chain A"
cerrada2$label[2480:4958] <- "Chain B"
cerrada2$label[4959:7433] <- "Chain A Last 100ns"
cerrada2$label[7434:9908] <- "Chain B Last 100ns"

cerrada2$label <- as.factor(cerrada2$label)
cerrada2$label <- ordered(cerrada2$label, levels=chains)

cerrada2_plot1 <- ggplot(cerrada2[1:4958,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 2")

cerrada2_plot2 <- ggplot(cerrada2[4959:9908,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 2")

v160m1_all_chainA <- readXVG("rdf/v160m1/v160m1_rdf_polsiteA.xvg")
v160m1_all_chainA <- as.data.frame(sapply(v160m1_all_chainA, as.numeric))
v160m1_all_chainB <- readXVG("rdf/v160m1/v160m1_rdf_polsiteB.xvg")
v160m1_all_chainB <- as.data.frame(sapply(v160m1_all_chainB, as.numeric))
v160m1_last_chainA <- readXVG("rdf/v160m1/v160m1_rdf_last100_polsiteA.xvg")
v160m1_last_chainA <- as.data.frame(sapply(v160m1_last_chainA, as.numeric))
v160m1_last_chainB <- readXVG("rdf/v160m1/v160m1_rdf_last100_polsiteB.xvg")
v160m1_last_chainB <- as.data.frame(sapply(v160m1_last_chainB, as.numeric))

v160m1 <- data.frame("r" = c(v160m1_all_chainA$r, v160m1_all_chainB$r, 
                             v160m1_last_chainA$r, v160m1_last_chainB$r), 
                     "g" = c(v160m1_all_chainA$`Water_&_OW`, v160m1_all_chainB$`Water_&_OW`,
                             v160m1_last_chainA$`Water_&_OW`, v160m1_last_chainB$`Water_&_OW`),
                     "label" = "label")

v160m1$label[1:2629] <- "Chain A"
v160m1$label[2630:5258] <- "Chain B"
v160m1$label[5259:7879] <- "Chain A Last 100ns"
v160m1$label[7880:10500] <- "Chain B Last 100ns"

v160m1$label <- as.factor(v160m1$label)
v160m1$label <- ordered(v160m1$label, levels=chains)

v160m1_plot1 <- ggplot(v160m1[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 1")
v160m1_plot2 <- ggplot(v160m1[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 1")


v160m2_all_chainA <- readXVG("rdf/v160m2/v160m2_rdf_polsiteA.xvg")
v160m2_all_chainA <- as.data.frame(sapply(v160m2_all_chainA, as.numeric))
v160m2_all_chainB <- readXVG("rdf/v160m2/v160m2_rdf_polsiteB.xvg")
v160m2_all_chainB <- as.data.frame(sapply(v160m2_all_chainB, as.numeric))
v160m2_last_chainA <- readXVG("rdf/v160m2/v160m2_rdf_last100_polsiteA.xvg")
v160m2_last_chainA <- as.data.frame(sapply(v160m2_last_chainA, as.numeric))
v160m2_last_chainB <- readXVG("rdf/v160m2/v160m2_rdf_last100_polsiteB.xvg")
v160m2_last_chainB <- as.data.frame(sapply(v160m2_last_chainB, as.numeric))

v160m2 <- data.frame("r" = c(v160m2_all_chainA$r, v160m2_all_chainB$r, 
                             v160m2_last_chainA$r, v160m2_last_chainB$r), 
                     "g" = c(v160m2_all_chainA$`Water_&_OW`, v160m2_all_chainB$`Water_&_OW`,
                             v160m2_last_chainA$`Water_&_OW`, v160m2_last_chainB$`Water_&_OW`),
                     "label" = "label")

v160m2$label[1:2629] <- "Chain A"
v160m2$label[2630:5258] <- "Chain B"
v160m2$label[5259:7879] <- "Chain A Last 100ns"
v160m2$label[7880:10500] <- "Chain B Last 100ns"

v160m2$label <- as.factor(v160m2$label)
v160m2$label <- ordered(v160m2$label, levels=chains)

v160m2_plot1 <- ggplot(v160m2[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 2")
v160m2_plot2 <- ggplot(v160m2[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 2")

ggarrange(abierta1_plot2, abierta2_plot2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(cerrada1_plot2, cerrada2_plot2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(v160m1_plot2, v160m2_plot2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")


#### Residue 286 #### 

abierta1_all_chainA <- readXVG("rdf/abierta1/abierta1_rdf_r286A.xvg")
abierta1_all_chainA <- as.data.frame(sapply(abierta1_all_chainA, as.numeric))
abierta1_all_chainB <- readXVG("rdf/abierta1/abierta1_rdf_r286B.xvg")
abierta1_all_chainB <- as.data.frame(sapply(abierta1_all_chainB, as.numeric))
abierta1_last_chainA <- readXVG("rdf/abierta1/abierta1_rdf_last100_r286A.xvg")
abierta1_last_chainA <- as.data.frame(sapply(abierta1_last_chainA, as.numeric))
abierta1_last_chainB <- readXVG("rdf/abierta1/abierta1_rdf_last100_r286B.xvg")
abierta1_last_chainB <- as.data.frame(sapply(abierta1_last_chainB, as.numeric))

abierta1 <- data.frame("r" = c(abierta1_all_chainA$r, abierta1_all_chainB$r, 
                               abierta1_last_chainA$r, abierta1_last_chainB$r), 
                       "g" = c(abierta1_all_chainA$`Water_&_OW`, abierta1_all_chainB$`Water_&_OW`,
                               abierta1_last_chainA$`Water_&_OW`, abierta1_last_chainB$`Water_&_OW`),
                       "label" = "label")

abierta1$label[1:2629] <- "Chain A"
abierta1$label[2630:5258] <- "Chain B"
abierta1$label[5259:7879] <- "Chain A Last 100ns"
abierta1$label[7880:10500] <- "Chain B Last 100ns"

chains <- c("Chain A", "Chain B", "Chain A Last 100ns", "Chain B Last 100ns")

abierta1$label <- as.factor(abierta1$label)
abierta1$label <- ordered(abierta1$label, levels=chains)


abierta1_plot1 <- ggplot(abierta1[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 1")
abierta1_plot2 <- ggplot(abierta1[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 1")


abierta2_all_chainA <- readXVG("rdf/abierta2/abierta2_rdf_r286A.xvg")
abierta2_all_chainA <- as.data.frame(sapply(abierta2_all_chainA, as.numeric))
abierta2_all_chainB <- readXVG("rdf/abierta2/abierta2_rdf_r286B.xvg")
abierta2_all_chainB <- as.data.frame(sapply(abierta2_all_chainB, as.numeric))
abierta2_last_chainA <- readXVG("rdf/abierta2/abierta2_rdf_last100_r286A.xvg")
abierta2_last_chainA <- as.data.frame(sapply(abierta2_last_chainA, as.numeric))
abierta2_last_chainB <- readXVG("rdf/abierta2/abierta2_rdf_last100_r286B.xvg")
abierta2_last_chainB <- as.data.frame(sapply(abierta2_last_chainB, as.numeric))

abierta2 <- data.frame("r" = c(abierta2_all_chainA$r, abierta2_all_chainB$r, 
                               abierta2_last_chainA$r, abierta2_last_chainB$r), 
                       "g" = c(abierta2_all_chainA$`Water_&_OW`, abierta2_all_chainB$`Water_&_OW`,
                               abierta2_last_chainA$`Water_&_OW`, abierta2_last_chainB$`Water_&_OW`),
                       "label" = "label")

abierta2$label[1:2629] <- "Chain A"
abierta2$label[2630:5258] <- "Chain B"
abierta2$label[5259:7879] <- "Chain A Last 100ns"
abierta2$label[7880:10500] <- "Chain B Last 100ns"

abierta2$label <- as.factor(abierta2$label)
abierta2$label <- ordered(abierta2$label, levels=chains)

abierta2_plot1 <- ggplot(abierta2[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) +  xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 2")
abierta2_plot2 <- ggplot(abierta2[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) +  xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING Replica 2")

cerrada1_all_chainA <- readXVG("rdf/cerrada1/cerrada1_rdf_r286A.xvg")
cerrada1_all_chainA <- as.data.frame(sapply(cerrada1_all_chainA, as.numeric))
cerrada1_all_chainB <- readXVG("rdf/cerrada1/cerrada1_rdf_r286B.xvg")
cerrada1_all_chainB <- as.data.frame(sapply(cerrada1_all_chainB, as.numeric))
cerrada1_last_chainA <- readXVG("rdf/cerrada1/cerrada1_rdf_last100_r286A.xvg")
cerrada1_last_chainA <- as.data.frame(sapply(cerrada1_last_chainA, as.numeric))
cerrada1_last_chainB <- readXVG("rdf/cerrada1/cerrada1_rdf_last100_r286B.xvg")
cerrada1_last_chainB <- as.data.frame(sapply(cerrada1_last_chainB, as.numeric))

cerrada1 <- data.frame("r" = c(cerrada1_all_chainA$r, cerrada1_all_chainB$r, 
                               cerrada1_last_chainA$r, cerrada1_last_chainB$r), 
                       "g" = c(cerrada1_all_chainA$`Water_&_OW`, cerrada1_all_chainB$`Water_&_OW`,
                               cerrada1_last_chainA$`Water_&_OW`, cerrada1_last_chainB$`Water_&_OW`),
                       "label" = "label")

cerrada1$label[1:2479] <- "Chain A"
cerrada1$label[2480:4958] <- "Chain B"
cerrada1$label[4959:7433] <- "Chain A Last 100ns"
cerrada1$label[7434:9908] <- "Chain B Last 100ns"

cerrada1$label <- as.factor(cerrada1$label)
cerrada1$label <- ordered(cerrada1$label, levels=chains)

cerrada1_plot1 <- ggplot(cerrada1[1:4958,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 1")
cerrada1_plot2 <- ggplot(cerrada1[4959:9908,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 1")


cerrada2_all_chainA <- readXVG("rdf/cerrada2/cerrada2_rdf_r286A.xvg")
cerrada2_all_chainA <- as.data.frame(sapply(cerrada2_all_chainA, as.numeric))
cerrada2_all_chainB <- readXVG("rdf/cerrada2/cerrada2_rdf_r286B.xvg")
cerrada2_all_chainB <- as.data.frame(sapply(cerrada2_all_chainB, as.numeric))
cerrada2_last_chainA <- readXVG("rdf/cerrada2/cerrada2_rdf_last100_r286A.xvg")
cerrada2_last_chainA <- as.data.frame(sapply(cerrada2_last_chainA, as.numeric))
cerrada2_last_chainB <- readXVG("rdf/cerrada2/cerrada2_rdf_last100_r286B.xvg")
cerrada2_last_chainB <- as.data.frame(sapply(cerrada2_last_chainB, as.numeric))

cerrada2 <- data.frame("r" = c(cerrada2_all_chainA$r, cerrada2_all_chainB$r, 
                               cerrada2_last_chainA$r, cerrada2_last_chainB$r), 
                       "g" = c(cerrada2_all_chainA$`Water_&_OW`, cerrada2_all_chainB$`Water_&_OW`,
                               cerrada2_last_chainA$`Water_&_OW`, cerrada2_last_chainB$`Water_&_OW`),
                       "label" = "label")

cerrada2$label[1:2479] <- "Chain A"
cerrada2$label[2480:4958] <- "Chain B"
cerrada2$label[4959:7433] <- "Chain A Last 100ns"
cerrada2$label[7434:9908] <- "Chain B Last 100ns"

cerrada2$label <- as.factor(cerrada2$label)
cerrada2$label <- ordered(cerrada2$label, levels=chains)

cerrada2_plot1 <- ggplot(cerrada2[1:4958,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 2")
cerrada2_plot2 <- ggplot(cerrada2[4959:9908,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Closed STING Replica 2")


v160m1_all_chainA <- readXVG("rdf/v160m1/v160m1_rdf_r286A.xvg")
v160m1_all_chainA <- as.data.frame(sapply(v160m1_all_chainA, as.numeric))
v160m1_all_chainB <- readXVG("rdf/v160m1/v160m1_rdf_r286B.xvg")
v160m1_all_chainB <- as.data.frame(sapply(v160m1_all_chainB, as.numeric))
v160m1_last_chainA <- readXVG("rdf/v160m1/v160m1_rdf_last100_r286A.xvg")
v160m1_last_chainA <- as.data.frame(sapply(v160m1_last_chainA, as.numeric))
v160m1_last_chainB <- readXVG("rdf/v160m1/v160m1_rdf_last100_r286B.xvg")
v160m1_last_chainB <- as.data.frame(sapply(v160m1_last_chainB, as.numeric))

v160m1 <- data.frame("r" = c(v160m1_all_chainA$r, v160m1_all_chainB$r, 
                             v160m1_last_chainA$r, v160m1_last_chainB$r), 
                     "g" = c(v160m1_all_chainA$`Water_&_OW`, v160m1_all_chainB$`Water_&_OW`,
                             v160m1_last_chainA$`Water_&_OW`, v160m1_last_chainB$`Water_&_OW`),
                     "label" = "label")

v160m1$label[1:2629] <- "Chain A"
v160m1$label[2630:5258] <- "Chain B"
v160m1$label[5259:7879] <- "Chain A Last 100ns"
v160m1$label[7880:10500] <- "Chain B Last 100ns"

v160m1$label <- as.factor(v160m1$label)
v160m1$label <- ordered(v160m1$label, levels=chains)

v160m1_plot1 <- ggplot(v160m1[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 1")
v160m1_plot2 <- ggplot(v160m1[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 1")


v160m2_all_chainA <- readXVG("rdf/v160m2/v160m2_rdf_r286A.xvg")
v160m2_all_chainA <- as.data.frame(sapply(v160m2_all_chainA, as.numeric))
v160m2_all_chainB <- readXVG("rdf/v160m2/v160m2_rdf_r286B.xvg")
v160m2_all_chainB <- as.data.frame(sapply(v160m2_all_chainB, as.numeric))
v160m2_last_chainA <- readXVG("rdf/v160m2/v160m2_rdf_last100_r286A.xvg")
v160m2_last_chainA <- as.data.frame(sapply(v160m2_last_chainA, as.numeric))
v160m2_last_chainB <- readXVG("rdf/v160m2/v160m2_rdf_last100_r286B.xvg")
v160m2_last_chainB <- as.data.frame(sapply(v160m2_last_chainB, as.numeric))

v160m2 <- data.frame("r" = c(v160m2_all_chainA$r, v160m2_all_chainB$r, 
                             v160m2_last_chainA$r, v160m2_last_chainB$r), 
                     "g" = c(v160m2_all_chainA$`Water_&_OW`, v160m2_all_chainB$`Water_&_OW`,
                             v160m2_last_chainA$`Water_&_OW`, v160m2_last_chainB$`Water_&_OW`),
                     "label" = "label")

v160m2$label[1:2629] <- "Chain A"
v160m2$label[2630:5258] <- "Chain B"
v160m2$label[5259:7879] <- "Chain A Last 100ns"
v160m2$label[7880:10500] <- "Chain B Last 100ns"

v160m2$label <- as.factor(v160m2$label)
v160m2$label <- ordered(v160m2$label, levels=chains)

v160m2_plot1 <- ggplot(v160m2[1:5258,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 2")
v160m2_plot2 <- ggplot(v160m2[5259:10500,], aes(x = r, y = g)) + 
  geom_line(aes(color = label), size = 0.5) + xlim(0,0.5) + 
  scale_color_manual(values = c("darkgreen", "red")) +
  xlab("r (nm)") + ylab("g(r)") + a + ggtitle("Open STING V160M Replica 2")

ggarrange(abierta1_plot1, abierta2_plot1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(cerrada1_plot1, cerrada2_plot1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggarrange(v160m1_plot1, v160m2_plot1, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
