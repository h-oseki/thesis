setwd("~/Desktop/Sysbio/cSTING")
index <- read.delim("~/Desktop/Sysbio/cSTING/corridas restringidas/index.txt")
index$index <- as.numeric(index$index)
index$X6nt6 <- as.character(index$X6nt6)
index$X6nt7_holo <- as.character(index$X6nt7_holo)

############################## RES 1 ##############################
res1_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res1.csv")
res1_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res1.csv")

res1_abierta <- res1_abierta[-1]
names(res1_abierta)[names(res1_abierta)=="Ocupación...."] <- "Ocurrencia"
res1_abierta[] <- lapply(res1_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res1_abierta[] <- lapply(res1_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res1_abierta$Ocurrencia <- as.numeric(res1_abierta$Ocurrencia)
res1_abierta <- separate(res1_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res1_abierta <- separate(res1_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res1_holo <- res1_holo[-1]
names(res1_holo)[names(res1_holo)=="Ocupación...."] <- "Ocurrencia"
res1_holo[] <- lapply(res1_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res1_holo[] <- lapply(res1_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res1_holo$Ocurrencia <- as.numeric(res1_holo$Ocurrencia)
res1_holo <- separate(res1_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res1_holo <- separate(res1_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:58){
  res1_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res1_abierta$Donor_res[i])] 
  res1_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res1_abierta$Aceptor_res[i])] 
}
for(i in 1:37){
  res1_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res1_holo$Donor_res[i])]
  res1_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res1_holo$Aceptor_res[i])]
}

res1 <- merge(x = res1_abierta, y = res1_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res1)[names(res1)=="Ocurrencia.x"] <- "Abierta"
names(res1)[names(res1)=="Ocurrencia.y"] <- "Holo"

res1$Abierta[which(is.na(res1$Abierta))] = 0
res1$Holo[which(is.na(res1$Holo))] = 0

res1$abMholo = res1$Abierta - res1$Holo
res1 <- res1[order(abs(res1$abMholo), decreasing = TRUE),]

res1ab <- subset(res1, abMholo > 0)
subset(res1ab, abMholo > 25)

res1ab_holo <- res1ab
res1ab_holo$Donor_res <- as.numeric(res1ab_holo$Donor_res)
res1ab_holo$Aceptor_res <- as.numeric(res1ab_holo$Aceptor_res)

for(i in 1:42){
  res1ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res1ab_holo$Donor_res[i])]
  res1ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res1ab_holo$Aceptor_res[i])]
}
subset(res1ab_holo, abMholo > 25)
residues2_holo <- c("271", "218", "439")
residues2_holo <- as.numeric(residues2_holo)
residues2_holo <- as.data.frame(residues2_holo)
write.csv(residues2_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)

res1ab_abierta <- res1ab
res1ab_abierta$Donor_res <- as.numeric(res1ab_abierta$Donor_res)
res1ab_abierta$Aceptor_res <- as.numeric(res1ab_abierta$Aceptor_res)

for(i in 1:42){
  res1ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res1ab_abierta$Donor_res[i])]
  res1ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res1ab_abierta$Aceptor_res[i])]
}
subset(res1ab_abierta, abMholo > 25)

residues2_abierta <- c("218", "271", "439")
residues2_abierta <- as.numeric(residues2_abierta)
residues2_abierta <- as.data.frame(residues2_abierta)
write.csv(residues2_abierta, "6nt6_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)

############################## RES 2 ##############################
res2_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res2.csv")
res2_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res2.csv")

res2_abierta <- res2_abierta[-1]
names(res2_abierta)[names(res2_abierta)=="Ocupación...."] <- "Ocurrencia"
res2_abierta[] <- lapply(res2_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res2_abierta[] <- lapply(res2_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res2_abierta$Ocurrencia <- as.numeric(res2_abierta$Ocurrencia)
res2_abierta <- separate(res2_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res2_abierta <- separate(res2_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res2_holo <- res2_holo[-1]
names(res2_holo)[names(res2_holo)=="Ocupación...."] <- "Ocurrencia"
res2_holo[] <- lapply(res2_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res2_holo[] <- lapply(res2_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res2_holo$Ocurrencia <- as.numeric(res2_holo$Ocurrencia)
res2_holo <- separate(res2_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res2_holo <- separate(res2_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:52){
  res2_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res2_abierta$Donor_res[i])] 
  res2_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res2_abierta$Aceptor_res[i])] 
}
for(i in 1:25){
  res2_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res2_holo$Donor_res[i])]
  res2_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res2_holo$Aceptor_res[i])]
}

res2 <- merge(x = res2_abierta, y = res2_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res2)[names(res2)=="Ocurrencia.x"] <- "Abierta"
names(res2)[names(res2)=="Ocurrencia.y"] <- "Holo"

res2$Abierta[which(is.na(res2$Abierta))] = 0
res2$Holo[which(is.na(res2$Holo))] = 0

res2$abMholo = res2$Abierta - res2$Holo
res2 <- res2[order(abs(res2$abMholo), decreasing = TRUE),]

res2ab <- subset(res2, abMholo > 0)
subset(res2ab, abMholo > 25)

res2ab_holo <- res2ab
res2ab_holo$Donor_res <- as.numeric(res2ab_holo$Donor_res)
res2ab_holo$Aceptor_res <- as.numeric(res2ab_holo$Aceptor_res)

for(i in 1:47){
  res2ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res2ab_holo$Donor_res[i])]
  res2ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res2ab_holo$Aceptor_res[i])]
}
subset(res2ab_holo, abMholo > 25)

residues3_holo <- c("213", "215", "250")
residues3_holo <- as.numeric(residues3_holo)
residues3_holo <- as.data.frame(residues3_holo)
write.csv(residues3_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)

res2ab_abierta <- res2ab
res2ab_abierta$Donor_res <- as.numeric(res2ab_abierta$Donor_res)
res2ab_abierta$Aceptor_res <- as.numeric(res2ab_abierta$Aceptor_res)

for(i in 1:47){
  res2ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res2ab_abierta$Donor_res[i])]
  res2ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res2ab_abierta$Aceptor_res[i])]
}
subset(res2ab_abierta, abMholo > 25)

residues3_abierta <- c("213", "215", "250")
residues3_abierta <- as.numeric(residues3_abierta)
residues3_abierta <- as.data.frame(residues3_abierta)
write.csv(residues3_abierta, "6nt6_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)

############################## RES 3 ##############################
res3_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res3.csv")
res3_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res3.csv")

res3_abierta <- res3_abierta[-1]
names(res3_abierta)[names(res3_abierta)=="Ocupación...."] <- "Ocurrencia"
res3_abierta[] <- lapply(res3_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res3_abierta[] <- lapply(res3_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res3_abierta$Ocurrencia <- as.numeric(res3_abierta$Ocurrencia)
res3_abierta <- separate(res3_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res3_abierta <- separate(res3_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res3_holo <- res3_holo[-1]
names(res3_holo)[names(res3_holo)=="Ocupación...."] <- "Ocurrencia"
res3_holo[] <- lapply(res3_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res3_holo[] <- lapply(res3_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res3_holo$Ocurrencia <- as.numeric(res3_holo$Ocurrencia)
res3_holo <- separate(res3_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res3_holo <- separate(res3_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:23){
  res3_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res3_abierta$Donor_res[i])] 
  res3_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res3_abierta$Aceptor_res[i])] 
}
for(i in 1:16){
  res3_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res3_holo$Donor_res[i])]
  res3_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res3_holo$Aceptor_res[i])]
}

res3 <- merge(x = res3_abierta, y = res3_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res3)[names(res3)=="Ocurrencia.x"] <- "Abierta"
names(res3)[names(res3)=="Ocurrencia.y"] <- "Holo"

res3$Abierta[which(is.na(res3$Abierta))] = 0
res3$Holo[which(is.na(res3$Holo))] = 0

res3$abMholo = res3$Abierta - res3$Holo
res3 <- res3[order(abs(res3$abMholo), decreasing = TRUE),]

res3ab <- subset(res3, abMholo > 0)
subset(res3ab, abMholo > 25)

res3ab_holo <- res3ab
res3ab_holo$Donor_res <- as.numeric(res3ab_holo$Donor_res)
res3ab_holo$Aceptor_res <- as.numeric(res3ab_holo$Aceptor_res)

for(i in 1:19){
  res3ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res3ab_holo$Donor_res[i])]
  res3ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res3ab_holo$Aceptor_res[i])]
}
subset(res3ab_holo, abMholo > 25)
#residues4_holo <- c("161")
#residues4_holo <- as.numeric(residues4_holo)
#residues4_holo <- as.data.frame(residues4_holo)
#write.csv(residues4_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)

res3ab_abierta <- res3ab
res3ab_abierta$Donor_res <- as.numeric(res3ab_abierta$Donor_res)
res3ab_abierta$Aceptor_res <- as.numeric(res3ab_abierta$Aceptor_res)

for(i in 1:19){
  res3ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res3ab_abierta$Donor_res[i])]
  res3ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res3ab_abierta$Aceptor_res[i])]
}
subset(res3ab_abierta, abMholo > 25)

#residues4_abierta <- c("358")
#residues4_abierta <- as.numeric(residues4_abierta)
#residues4_abierta <- as.data.frame(residues4_abierta)
#write.csv(residues4_abierta, "6nt6_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)
