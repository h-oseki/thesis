hbonds <- read.csv("~/Desktop/Sysbio/cSTING/hbonds-from-ligand2.csv")
ggplot(hbonds[c(1,13),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) + ylim(0,100)
hbonds[c(1,13),]
ggplot(hbonds[c(2,14),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(2,14),]
ggplot(hbonds[c(3,15),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(3,15),]
ggplot(hbonds[c(4,16),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(4,16),]
ggplot(hbonds[c(5,17),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(5,17),]
ggplot(hbonds[c(6,18),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(6,18),]
ggplot(hbonds[c(7,19),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(7,19),]
ggplot(hbonds[c(8,20),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(8,20),]
ggplot(hbonds[c(9,21),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(9,21),]
ggplot(hbonds[c(10,22),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(10,22),]
ggplot(hbonds[c(11,23),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(11,23),]
ggplot(hbonds[c(12,24),], aes(hb, occ, fill = str)) + geom_col(position = "dodge2") + 
  geom_hline(yintercept = 0) + scale_fill_manual(values = c("mediumseagreen", "mediumorchid2")) +
  theme(panel.grid.major = element_line(color =NA), 
        panel.grid.minor = element_line(color =NA),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")  + ylim(0,100)
hbonds[c(12,24),]
