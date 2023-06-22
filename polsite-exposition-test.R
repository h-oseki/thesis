holo_chainA <- read.delim("6nt7_HOLO_NO_BORRAR/distances/chainA-holo.xvg")
holo_chainB <- read.delim("6nt7_HOLO_NO_BORRAR/distances/chainB-holo.xvg")
holo_protein <- read.delim("6nt7_HOLO_NO_BORRAR/distances/protein-holo.xvg")
rot_chainA <- read.delim("rot_cgamp/production/distances/chainA-rot.xvg")
rot_chainB <- read.delim("rot_cgamp/production/distances/chainB-rot.xvg")
rot_protein <- read.delim("rot_cgamp/production/distances/protein-rot.xvg")

holo <- data.frame("time" = holo_chainA$time/1000, 
                   "chainA_x" = holo_chainA$x,"chainA_y" = holo_chainA$y, "chainA_z" = holo_chainA$z,
                   "chainB_x" = holo_chainB$x, "chainB_y" = holo_chainB$y,"chainB_z" = holo_chainB$z, 
                   "protein_x" = holo_protein$x,"protein_y" = holo_protein$y, "protein_z" = holo_protein$z)

rot <- data.frame("time" = holo_chainA$time/1000, 
                   "chainA_x" = rot_chainA$x,"chainA_y" = rot_chainA$y, "chainA_z" = rot_chainA$z,
                   "chainB_x" = rot_chainB$x, "chainB_y" = rot_chainB$y,"chainB_z" = rot_chainB$z, 
                   "protein_x" = rot_protein$x,"protein_y" = rot_protein$y, "protein_z" = rot_protein$z)
holo$chainA_dist <- sqrt((holo$chainA_x - holo$protein_x)^2 + 
                           (holo$chainA_y - holo$protein_y)^2 + 
                           (holo$chainA_z - holo$protein_z)^2)
holo$chainB_dist <- sqrt((holo$chainB_x - holo$protein_x)^2 + 
                           (holo$chainB_y - holo$protein_y)^2 + 
                           (holo$chainB_z - holo$protein_z)^2)
rot$chainA_dist <- sqrt((rot$chainA_x - rot$protein_x)^2 + 
                           (rot$chainA_y - rot$protein_y)^2 + 
                           (rot$chainA_z - rot$protein_z)^2)
rot$chainB_dist <- sqrt((rot$chainB_x - rot$protein_x)^2 + 
                           (rot$chainB_y - rot$protein_y)^2 + 
                           (rot$chainB_z - rot$protein_z)^2)

holo$normdistA <- holo$chainA_dist - holo$chainA_dist[1]
holo$normdistB <- holo$chainB_dist - holo$chainB_dist[1]
rot$normdistA <- rot$chainA_dist - rot$chainA_dist[1]
rot$normdistB <- rot$chainB_dist - rot$chainB_dist[1]

distances2 <- data.frame("Time" = c(rot$time, rot$time, rot$time, rot$time),
                        "Distance" = c(holo$chainA_dist, holo$chainB_dist, 
                                       rot$chainA_dist, rot$chainB_dist),
                        "Norm_distance" = c(holo$normdistA, holo$normdistB, 
                                          rot$normdistA, rot$normdistB),
                        "Label" = "label")
distances2$Label[1:2501] <- "HOLO - chain A"
distances2$Label[2502:5002] <- "HOLO - chain B"
distances2$Label[5003:7503] <- "ROT - chain A"
distances2$Label[7504:10004] <- "ROT - chain B"

ggplot(distances2[1:5002,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("deeppink", "chartreuse3", "dodgerblue3", "purple")) + a + 
  ggtitle("Distance between the polymerization site and the center of mass of STING")

ggplot(distances2[5003:10004,], aes(x = Time, y = Distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("dodgerblue3", "purple")) + a + 
  ggtitle("Distance between the polymerization site and the center of mass of STING")

control_holo <- ggplot(distances2[1:5002,], aes(x = Time, y = Norm_distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("dodgerblue3", "purple")) + a + 
  ggtitle("Normalized distance between the polymerization site and the center of mass of STING") +
  annotate("rect", xmin = 25, xmax=50, ymin= -0.2, ymax=0.45, alpha=.1, fill="deeppink") +
  annotate("rect", xmin = 115, xmax=200, ymin= -0.2, ymax=0.45, alpha=.1, fill="deeppink") +
  annotate("rect", xmin = 360, xmax=425, ymin= -0.2, ymax=0.45, alpha=.1, fill="deeppink") 

control_rot <- ggplot(distances2[5003:10004,], aes(x = Time, y = Norm_distance, color = Label)) +
  geom_line(size = .8) + xlab("Time (ns)") + ylab("Distance (nm)") + 
  scale_color_manual(values = c("dodgerblue3", "purple")) + a + 
  ggtitle("Normalized distance between the polymerization site and the center of mass of STING") +
  annotate("rect", xmin = 90, xmax=275, ymin= -0.2, ymax=0.45, alpha=.1, fill="deeppink") +
  annotate("rect", xmin = 325, xmax=390, ymin= -0.2, ymax=0.45, alpha=.1, fill="deeppink") 




