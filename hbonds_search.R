abierta1 <- read.delim("~/Desktop/Sysbio/cSTING/abierta-rep1/crop-200/hbonds/hbonds-details-abierta1.dat")
abierta2 <- read.delim("~/Desktop/Sysbio/cSTING/abierta-rep2/hbond/hbonds-details_abierta2.dat")
close1 <- read.delim("~/Desktop/Sysbio/cSTING/cerrada-rep1/prod/hbonds/hbonds-details-close1.dat")
close2 <- read.delim("~/Desktop/Sysbio/cSTING/cerrada-rep2/hbonds/hbonds-details-close2.dat")
mutant1 <- read.delim("~/Desktop/Sysbio/cSTING/v160m-rep1/prod_200ns/hbonds/hbonds-details.dat")
mutant2 <- read.delim("~/Desktop/Sysbio/cSTING/v160m-rep2/hbonds/hbonds-details.dat")

abierta1[] <- lapply(abierta1, gsub, pattern = " ", replacement = "", fixed = TRUE)
abierta2[] <- lapply(abierta2, gsub, pattern = " ", replacement = "", fixed = TRUE)
close1[] <- lapply(close1, gsub, pattern = " ", replacement = "", fixed = TRUE)
close2[] <- lapply(close2, gsub, pattern = " ", replacement = "", fixed = TRUE)
mutant1[] <- lapply(mutant1, gsub, pattern = " ", replacement = "", fixed = TRUE)
mutant2[] <- lapply(mutant2, gsub, pattern = " ", replacement = "", fixed = TRUE)

abierta1[] <- lapply(abierta1, gsub, pattern = "%", replacement = "", fixed = TRUE)
abierta2[] <- lapply(abierta2, gsub, pattern = "%", replacement = "", fixed = TRUE)
close1[] <- lapply(close1, gsub, pattern = "%", replacement = "", fixed = TRUE)
close2[] <- lapply(close2, gsub, pattern = "%", replacement = "", fixed = TRUE)
mutant1[] <- lapply(mutant1, gsub, pattern = "%", replacement = "", fixed = TRUE)
mutant2[] <- lapply(mutant2, gsub, pattern = "%", replacement = "", fixed = TRUE)

abierta1$occupancy <- as.numeric(abierta1$occupancy)
abierta2$occupancy <- as.numeric(abierta2$occupancy)
close1$occupancy <- as.numeric(close1$occupancy)
close2$occupancy <- as.numeric(close2$occupancy)
mutant1$occupancy <- as.numeric(mutant1$occupancy)
mutant2$occupancy <- as.numeric(mutant2$occupancy)

abierta1 <- subset(abierta1, occupancy > 25)
abierta2 <- subset(abierta2, occupancy > 25)
close1 <- subset(close1, occupancy > 25)
close2 <- subset(close2, occupancy > 25)
mutant1 <- subset(mutant1, occupancy > 25)
mutant2 <- subset(mutant2, occupancy > 25)

abierta1 <- separate(abierta1, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
abierta1 <- separate(abierta1, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
abierta2 <- separate(abierta2, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
abierta2 <- separate(abierta2, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
close1 <- separate(close1, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
close1 <- separate(close1, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
close2 <- separate(close2, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
close2 <- separate(close2, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
mutant1 <- separate(mutant1, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
mutant1 <- separate(mutant1, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
mutant2 <- separate(mutant2, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
mutant2 <- separate(mutant2, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)

show_hbonds <- function(db, residue){
  donor <- which(db$Donor_res == residue)
  aceptor <- which(db$Aceptor_res == residue)
  return(db[c(donor,aceptor),c(1,3,4,6,7)])
}

residues <- c("ALA146", "VAL147", "GLU148", "VAL149", "SER150", "GLU151", "LEU152", "THR153", 
              "GLU154", "SER155", "SER156", "LYS157", "LYS158", "ASN159", "VAL160", "ALA161", 
              "HIS162", "GLY163", "LEU164", "ALA165", "TRP166", "SER167", "TYR168", "TYR169", 
              "ILE170", "GLY171", "TYR172", "LEU173", "LYS174", "VAL175", "VAL176", "LEU177", 
              "PRO178", "ARG179", "LEU180", "ALA343", "VAL344", "GLU345", "VAL346", "SER347", 
              "GLU348", "LEU349", "THR350", "GLU351", "SER352", "SER353", "LYS354", "LYS355", 
              "ASN356", "VAL357", "ALA358", "HIS359", "GLY360", "LEU361", "ALA362", "TRP363", 
              "SER364", "TYR365", "TYR366", "ILE367", "GLY368", "TYR369", "LEU370", "LYS371", 
              "VAL372", "VAL373", "LEU374", "PRO375", "ARG376", "LEU377")
residues <- as.data.frame(residues)
residues[2,1]
show_hbonds(mutant1, "GLN292")
show_hbonds(mutant2, "GLN292")

show_hbonds(mutant1, "GLN489")
show_hbonds(mutant2, "GLN489")


closed_hbonds$occupancy.x[which(is.na(closed_hbonds$occupancy.x))] = 0
closed_hbonds$occupancy.y[which(is.na(closed_hbonds$occupancy.y))] = 0



for(i in 1:70){
  if(i == 1){
    abierta1_HB <- show_hbonds(abierta1, as.character(residues[i,1]))
  }
  else{
    aux <- show_hbonds(abierta1, as.character(residues[i,1]))
    abierta1_HB <- rbind(abierta1_HB, aux)
    rm(aux)
  }
}
  
for(i in 1:70){
  if(i == 1){
    abierta2_HB <- show_hbonds(abierta2, as.character(residues[i,1]))
  }
  else{
    aux <- show_hbonds(abierta2, as.character(residues[i,1]))
    abierta2_HB <- rbind(abierta2_HB, aux)
    rm(aux)
  }
}

for(i in 1:70){
  if(i == 1){
    cerrada1_HB <- show_hbonds(close1, as.character(residues[i,1]))
  }
  else{
    aux <- show_hbonds(close1, as.character(residues[i,1]))
    cerrada1_HB <- rbind(cerrada1_HB, aux)
    rm(aux)
  }
}

for(i in 1:70){
  if(i == 1){
    cerrada2_HB <- show_hbonds(close2, as.character(residues[i,1]))
  }
  else{
    aux <- show_hbonds(close2, as.character(residues[i,1]))
    cerrada2_HB <- rbind(cerrada2_HB, aux)
    rm(aux)
  }
}
  

hbonds <- merge(x = abierta1_HB, y = abierta2_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy.x"] <- "Abierta1"
names(hbonds)[names(hbonds)=="occupancy.y"] <- "Abierta2"
hbonds <- merge(x = hbonds, y = cerrada1_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy"] <- "Cerrada1"
hbonds <- merge(x = hbonds, y = cerrada2_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy"] <- "Cerrada2"

hbonds$Abierta1[which(is.na(hbonds$Abierta1))] = 0
hbonds$Abierta2[which(is.na(hbonds$Abierta2))] = 0
hbonds$Cerrada1[which(is.na(hbonds$Cerrada1))] = 0
hbonds$Cerrada2[which(is.na(hbonds$Cerrada2))] = 0

hbonds <- distinct(hbonds)

hbonds$open_abs = abs(hbonds$Abierta1 - hbonds$Abierta2)
hbonds$closed_abs = abs(hbonds$Cerrada1 - hbonds$Cerrada2)

hbonds <- subset(hbonds, open_abs < 10)
hbonds <- subset(hbonds, closed_abs < 10)

hbonds$abierta_av <- 0
hbonds$cerrada_av <- 0

for(i in 1:1055){
  hbonds$abierta_av[i] <- mean(hbonds$Abierta1[i], hbonds$Abierta2[i])
  hbonds$cerrada_av[i] <- mean(hbonds$Cerrada1[i], hbonds$Cerrada2[i])
}

hbonds$ab_cl <- hbonds$abierta_av - hbonds$cerrada_av
hbonds$resta_abs <- abs(hbonds$ab_cl)

hbonds2 <- subset(hbonds, resta_abs > 25)

hbonds3 <- merge(x = abierta1, y = abierta2, by = c("Donor_res","Donor_part",  "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds3)[names(hbonds3)=="occupancy.x"] <- "Abierta1"
names(hbonds3)[names(hbonds3)=="occupancy.y"] <- "Abierta2"
hbonds3 <- merge(x = hbonds3, y = close1, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds3)[names(hbonds3)=="occupancy"] <- "Cerrada1"
hbonds3 <- merge(x = hbonds3, y = close2, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part","Aceptor_atom"), all = TRUE)
names(hbonds3)[names(hbonds3)=="occupancy"] <- "Cerrada2"

hbonds3$Abierta1[which(is.na(hbonds3$Abierta1))] = 0
hbonds3$Abierta2[which(is.na(hbonds3$Abierta2))] = 0
hbonds3$Cerrada1[which(is.na(hbonds3$Cerrada1))] = 0
hbonds3$Cerrada2[which(is.na(hbonds3$Cerrada2))] = 0

hbonds3 <- distinct(hbonds3)

hbonds3$open_abs = abs(hbonds3$Abierta1 - hbonds3$Abierta2)
hbonds3$closed_abs = abs(hbonds3$Cerrada1 - hbonds3$Cerrada2)

hbonds3 <- subset(hbonds3, open_abs < 10)
hbonds3 <- subset(hbonds3, closed_abs < 10)

hbonds3$abierta_av <- 0
hbonds3$cerrada_av <- 0

for(i in 1:4821){
  hbonds3$abierta_av[i] <- mean(hbonds3$Abierta1[i], hbonds3$Abierta2[i])
  hbonds3$cerrada_av[i] <- mean(hbonds3$Cerrada1[i], hbonds3$Cerrada2[i])
}

hbonds3$ab_cl <- hbonds3$abierta_av - hbonds3$cerrada_av
hbonds3$resta_abs <- abs(hbonds3$ab_cl)

hbonds3 <- subset(hbonds3, resta_abs > 25)

hbonds3[which(hbonds3$ab_cl > 0), c(1,4, 13:15)]

hbonds3[which(hbonds3$ab_cl < 0), c(1,4, 13:15)]


##########

hbonds4 <- merge(x = abierta1, y = abierta2, by = c("Donor_res","Donor_part",  "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds4)[names(hbonds4)=="occupancy.x"] <- "Abierta1"
names(hbonds4)[names(hbonds4)=="occupancy.y"] <- "Abierta2"
hbonds4 <- merge(x = hbonds4, y = close1, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds4)[names(hbonds4)=="occupancy"] <- "Cerrada1"
hbonds4 <- merge(x = hbonds4, y = close2, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part","Aceptor_atom"), all = TRUE)
names(hbonds4)[names(hbonds4)=="occupancy"] <- "Cerrada2"

hbonds4$Abierta1[which(is.na(hbonds4$Abierta1))] = 0
hbonds4$Abierta2[which(is.na(hbonds4$Abierta2))] = 0
hbonds4$Cerrada1[which(is.na(hbonds4$Cerrada1))] = 0
hbonds4$Cerrada2[which(is.na(hbonds4$Cerrada2))] = 0

hbonds4 <- distinct(hbonds4)

hbonds4$open_abs = abs(hbonds4$Abierta1 - hbonds4$Abierta2)
hbonds4$closed_abs = abs(hbonds4$Cerrada1 - hbonds4$Cerrada2)

hbonds4 <- subset(hbonds4, open_abs < 10)
hbonds4 <- subset(hbonds4, closed_abs < 10)

hbonds4$abierta_av <- 0
hbonds4$cerrada_av <- 0

for(i in 1:4821){
  hbonds4$abierta_av[i] <- mean(hbonds4$Abierta1[i], hbonds4$Abierta2[i])
  hbonds4$cerrada_av[i] <- mean(hbonds4$Cerrada1[i], hbonds4$Cerrada2[i])
}

hbonds4 <- subset(hbonds4, abierta_av > 10)
hbonds4 <- subset(hbonds4, cerrada_av > 10)

hbonds4$ab_cl <- hbonds4$abierta_av - hbonds4$cerrada_av
hbonds4$resta_abs <- abs(hbonds4$ab_cl)

hbonds4 <- subset(hbonds4, resta_abs < 5)


hbonds4[, c(1,4, 13:15)]


########

residues2 <- c("LEU304", "GLY305", "SER306", "SER307", "LYS308", "GLU309", "CYS310", "ALA311",
               "LEU501", "GLY502", "SER503", "SER504", "LYS505", "GLU506", "CYS507", "ALA508")
residues2 <- as.data.frame(residues2)
residues[2,1]
show_hbonds(close2, "GLU309")


for(i in 1:16){
  if(i == 1){
    abierta1_HB <- show_hbonds(abierta1, as.character(residues2[i,1]))
  }
  else{
    aux <- show_hbonds(abierta1, as.character(residues2[i,1]))
    abierta1_HB <- rbind(abierta1_HB, aux)
    rm(aux)
  }
}

for(i in 1:16){
  if(i == 1){
    abierta2_HB <- show_hbonds(abierta2, as.character(residues2[i,1]))
  }
  else{
    aux <- show_hbonds(abierta2, as.character(residues2[i,1]))
    abierta2_HB <- rbind(abierta2_HB, aux)
    rm(aux)
  }
}

for(i in 1:16){
  if(i == 1){
    cerrada1_HB <- show_hbonds(close1, as.character(residues2[i,1]))
  }
  else{
    aux <- show_hbonds(close1, as.character(residues2[i,1]))
    cerrada1_HB <- rbind(cerrada1_HB, aux)
    rm(aux)
  }
}

for(i in 1:16){
  if(i == 1){
    cerrada2_HB <- show_hbonds(close2, as.character(residues2[i,1]))
  }
  else{
    aux <- show_hbonds(close2, as.character(residues2[i,1]))
    cerrada2_HB <- rbind(cerrada2_HB, aux)
    rm(aux)
  }
}


hbonds <- merge(x = abierta1_HB, y = abierta2_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy.x"] <- "Abierta1"
names(hbonds)[names(hbonds)=="occupancy.y"] <- "Abierta2"
hbonds <- merge(x = hbonds, y = cerrada1_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy"] <- "Cerrada1"
hbonds <- merge(x = hbonds, y = cerrada2_HB, by = c("Donor_res", "Donor_atom", "Aceptor_res", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy"] <- "Cerrada2"

hbonds$Abierta1[which(is.na(hbonds$Abierta1))] = 0
hbonds$Abierta2[which(is.na(hbonds$Abierta2))] = 0
hbonds$Cerrada1[which(is.na(hbonds$Cerrada1))] = 0
hbonds$Cerrada2[which(is.na(hbonds$Cerrada2))] = 0

