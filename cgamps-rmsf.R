cgamp <- read.csv("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/cgamp-rmsf.xvg", sep="")
rot_cgamp <- read.csv("~/Desktop/Sysbio/cSTING/rot_cgamp/production/cgamp-rot-rmsf.xvg", sep="")

cgamps <- merge(x = cgamp[,-1], y = rot_cgamp[,-1], by = c("atom"), all = TRUE)
names(cgamps)[names(cgamps)=="rmsf.x"] <- "cgamp"
names(cgamps)[names(cgamps)=="rmsf.y"] <- "rot"
cgamps$num <- c(1:67)
  
cgamps2 <- data.frame("num" = c(cgamps$num, cgamps$num),
                      "rmsf" = c(cgamps$cgamp, cgamps$rot),
                      "label" = "label")
cgamps2$label[1:67] <- "cgamp"
cgamps2$label[68:134] <- "rot-cgamp"

ggplot(cgamps2, aes(x = num, y = rmsf)) + 
  geom_line(aes(color = label, linetype = label), size = 0.8) + 
  scale_color_manual(values = c("mediumorchid2", "dodgerblue3")) +
  scale_linetype_manual(values=c("twodash", "solid")) + a +
  ylab("B-factor") +  
  scale_x_continuous(breaks = c(1:67),
                     labels = paste(c("C1", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", 
                                      "C19", "C2", "C20", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "H1", 
                                      "H10", "H11", "H12", "H13", "H14", "H15", "H16", "H17", "H18", "H19", 
                                      "H2", "H20", "H21", "H22", "H3", "H4", "H5", "H6", "H7","H8", "H9", 
                                      "N1", "N10", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "O1", 
                                      "O10", "O11", "O12", "O13", "O2", "O3", "O4", "O5", "O6", "O7", "O8", 
                                      "O9", "P1", "P2"))) +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", size = 7),
    panel.grid.major = element_line(color = NA), 
    panel.grid.minor = element_line(color = NA),
    plot.title = element_text(colour = 'black', size = 15),
    axis.title = element_text(colour = 'black', face = 'italic'),
    panel.background = element_rect(fill = "white", colour = "black",),
    axis.title.x = element_blank()) 
