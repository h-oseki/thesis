setwd("~/Desktop/Sysbio/cSTING")
index <- read.delim("~/Desktop/Sysbio/cSTING/corridas restringidas/index.txt")
index$index <- as.numeric(index$index)
index$X6nt6 <- as.character(index$X6nt6)
index$X6nt7_holo <- as.character(index$X6nt7_holo)

############################## RES 1 ##############################
res1_abierta <- read.csv("6nt6_NO_BORRAR/hbonds-polsite-positive/part1/hbonds_res1.csv")
res1_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds-polsite-positive/part1/hbonds_res1.csv")

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

for(i in 1:112){
  res1_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res1_abierta$Donor_res[i])] 
  res1_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res1_abierta$Aceptor_res[i])] 
}
for(i in 1:112){
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

for(i in 1:80){
  res1ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res1ab_holo$Donor_res[i])]
  res1ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res1ab_holo$Aceptor_res[i])]
}
subset(res1ab_holo, abMholo > 25)
residues2_holo <- c("308", "366", "371", "472", "473")
residues2_holo <- as.numeric(residues2_holo)
residues2_holo <- as.data.frame(residues2_holo)
write.csv(residues2_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)

res1ab_abierta <- res1ab
res1ab_abierta$Donor_res <- as.numeric(res1ab_abierta$Donor_res)
res1ab_abierta$Aceptor_res <- as.numeric(res1ab_abierta$Aceptor_res)

for(i in 1:80){
  res1ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res1ab_abierta$Donor_res[i])]
  res1ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res1ab_abierta$Aceptor_res[i])]
}
subset(res1ab_abierta, abMholo > 30)

residues2_abierta <- c("308", "366", "371", "472", "473")
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

for(i in 1:36){
  res2_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res2_abierta$Donor_res[i])] 
  res2_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res2_abierta$Aceptor_res[i])] 
}
for(i in 1:34){
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

for(i in 1:24){
  res2ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res2ab_holo$Donor_res[i])]
  res2ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res2ab_holo$Aceptor_res[i])]
}
subset(res2ab_holo, abMholo > 25)

residues3_holo <- c("362", "468")
residues3_holo <- as.numeric(residues3_holo)
residues3_holo <- as.data.frame(residues3_holo)
write.csv(residues3_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)

res2ab_abierta <- res2ab
res2ab_abierta$Donor_res <- as.numeric(res2ab_abierta$Donor_res)
res2ab_abierta$Aceptor_res <- as.numeric(res2ab_abierta$Aceptor_res)

for(i in 1:24){
  res2ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res2ab_abierta$Donor_res[i])]
  res2ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res2ab_abierta$Aceptor_res[i])]
}
subset(res2ab_abierta, abMholo > 25)

residues3_abierta <- c("362", "468")
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

for(i in 1:21){
  res3_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res3_abierta$Donor_res[i])] 
  res3_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res3_abierta$Aceptor_res[i])] 
}
for(i in 1:17){
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

for(i in 1:17){
  res3ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res3ab_holo$Donor_res[i])]
  res3ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res3ab_holo$Aceptor_res[i])]
}
subset(res3ab_holo, abMholo > 25)
residues4_holo <- c("161")
residues4_holo <- as.numeric(residues4_holo)
residues4_holo <- as.data.frame(residues4_holo)
write.csv(residues4_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)

res3ab_abierta <- res3ab
res3ab_abierta$Donor_res <- as.numeric(res3ab_abierta$Donor_res)
res3ab_abierta$Aceptor_res <- as.numeric(res3ab_abierta$Aceptor_res)

for(i in 1:17){
  res3ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res3ab_abierta$Donor_res[i])]
  res3ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res3ab_abierta$Aceptor_res[i])]
}
subset(res3ab_abierta, abMholo > 25)

residues4_abierta <- c("358")
residues4_abierta <- as.numeric(residues4_abierta)
residues4_abierta <- as.data.frame(residues4_abierta)
write.csv(residues4_abierta, "6nt6_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)

############################## RES 4 ##############################
res4_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res4.csv")
res4_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res4.csv")

res4_abierta <- res4_abierta[-1]
names(res4_abierta)[names(res4_abierta)=="Ocupación...."] <- "Ocurrencia"
res4_abierta[] <- lapply(res4_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res4_abierta[] <- lapply(res4_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res4_abierta$Ocurrencia <- as.numeric(res4_abierta$Ocurrencia)
res4_abierta <- separate(res4_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res4_abierta <- separate(res4_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res4_holo <- res4_holo[-1]
names(res4_holo)[names(res4_holo)=="Ocupación...."] <- "Ocurrencia"
res4_holo[] <- lapply(res4_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res4_holo[] <- lapply(res4_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res4_holo$Ocurrencia <- as.numeric(res4_holo$Ocurrencia)
res4_holo <- separate(res4_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res4_holo <- separate(res4_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:4){
  res4_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res4_abierta$Donor_res[i])] 
  res4_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res4_abierta$Aceptor_res[i])] 
}
for(i in 1:5){
  res4_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res4_holo$Donor_res[i])]
  res4_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res4_holo$Aceptor_res[i])]
}

res4 <- merge(x = res4_abierta, y = res4_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res4)[names(res4)=="Ocurrencia.x"] <- "Abierta"
names(res4)[names(res4)=="Ocurrencia.y"] <- "Holo"

res4$Abierta[which(is.na(res4$Abierta))] = 0
res4$Holo[which(is.na(res4$Holo))] = 0

res4$abMholo = res4$Abierta - res4$Holo
res4 <- res4[order(abs(res4$abMholo), decreasing = TRUE),]

res4ab <- subset(res4, abMholo > 0)
subset(res4ab, abMholo > 25)

res4ab_holo <- res4ab
res4ab_holo$Donor_res <- as.numeric(res4ab_holo$Donor_res)
res4ab_holo$Aceptor_res <- as.numeric(res4ab_holo$Aceptor_res)

for(i in 1:3){
  res4ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res4ab_holo$Donor_res[i])]
  res4ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res4ab_holo$Aceptor_res[i])]
}
subset(res4ab_holo, abMholo > 25)
residues5_holo <- c("489")
residues5_holo <- as.numeric(residues5_holo)
residues5_holo <- as.data.frame(residues5_holo)
write.csv(residues5_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res5.txt", row.names = FALSE)

res4ab_abierta <- res4ab
res4ab_abierta$Donor_res <- as.numeric(res4ab_abierta$Donor_res)
res4ab_abierta$Aceptor_res <- as.numeric(res4ab_abierta$Aceptor_res)

for(i in 1:3){
  res4ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res4ab_abierta$Donor_res[i])]
  res4ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res4ab_abierta$Aceptor_res[i])]
}
subset(res4ab_abierta, abMholo > 25)

residues5_abierta <- c("489")
residues5_abierta <- as.numeric(residues5_abierta)
residues5_abierta <- as.data.frame(residues5_abierta)
write.csv(residues5_abierta, "6nt6_NO_BORRAR/hbonds/res5.txt", row.names = FALSE)
############################## RES 5 ##############################
res5_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res5.csv")
res5_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res5.csv")

res5_abierta <- res5_abierta[-1]
names(res5_abierta)[names(res5_abierta)=="Ocupación...."] <- "Ocurrencia"
res5_abierta[] <- lapply(res5_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res5_abierta[] <- lapply(res5_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res5_abierta$Ocurrencia <- as.numeric(res5_abierta$Ocurrencia)
res5_abierta <- separate(res5_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res5_abierta <- separate(res5_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res5_holo <- res5_holo[-1]
names(res5_holo)[names(res5_holo)=="Ocupación...."] <- "Ocurrencia"
res5_holo[] <- lapply(res5_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res5_holo[] <- lapply(res5_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res5_holo$Ocurrencia <- as.numeric(res5_holo$Ocurrencia)
res5_holo <- separate(res5_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res5_holo <- separate(res5_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:9){
  res5_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res5_abierta$Donor_res[i])] 
  res5_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res5_abierta$Aceptor_res[i])] 
}
for(i in 1:10){
  res5_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res5_holo$Donor_res[i])]
  res5_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res5_holo$Aceptor_res[i])]
}

res5 <- merge(x = res5_abierta, y = res5_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res5)[names(res5)=="Ocurrencia.x"] <- "Abierta"
names(res5)[names(res5)=="Ocurrencia.y"] <- "Holo"

res5$Abierta[which(is.na(res5$Abierta))] = 0
res5$Holo[which(is.na(res5$Holo))] = 0

res5$abMholo = res5$Abierta - res5$Holo
res5 <- res5[order(abs(res5$abMholo), decreasing = TRUE),]

res5ab <- subset(res5, abMholo > 0)
subset(res5ab, abMholo > 25)

res5ab_holo <- res5ab
res5ab_holo$Donor_res <- as.numeric(res5ab_holo$Donor_res)
res5ab_holo$Aceptor_res <- as.numeric(res5ab_holo$Aceptor_res)

for(i in 1:8){
  res5ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res5ab_holo$Donor_res[i])]
  res5ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res5ab_holo$Aceptor_res[i])]
}
subset(res5ab_holo, abMholo > 25)
residues6_holo <- c("159", "160")
residues6_holo <- as.numeric(residues6_holo)
residues6_holo <- as.data.frame(residues6_holo)
write.csv(residues6_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res6.txt", row.names = FALSE)

res5ab_abierta <- res5ab
res5ab_abierta$Donor_res <- as.numeric(res5ab_abierta$Donor_res)
res5ab_abierta$Aceptor_res <- as.numeric(res5ab_abierta$Aceptor_res)

for(i in 1:8){
  res5ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res5ab_abierta$Donor_res[i])]
  res5ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res5ab_abierta$Aceptor_res[i])]
}
subset(res5ab_abierta, abMholo > 25)

residues6_abierta <- c("356", "357")
residues6_abierta <- as.numeric(residues6_abierta)
residues6_abierta <- as.data.frame(residues6_abierta)
write.csv(residues6_abierta, "6nt6_NO_BORRAR/hbonds/res6.txt", row.names = FALSE)


############################## RES 6 ##############################
res6_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res6.csv")
res6_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res6.csv")

res6_abierta <- res6_abierta[-1]
names(res6_abierta)[names(res6_abierta)=="Ocupación...."] <- "Ocurrencia"
res6_abierta[] <- lapply(res6_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res6_abierta[] <- lapply(res6_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res6_abierta$Ocurrencia <- as.numeric(res6_abierta$Ocurrencia)
res6_abierta <- separate(res6_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res6_abierta <- separate(res6_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res6_holo <- res6_holo[-1]
names(res6_holo)[names(res6_holo)=="Ocupación...."] <- "Ocurrencia"
res6_holo[] <- lapply(res6_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res6_holo[] <- lapply(res6_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res6_holo$Ocurrencia <- as.numeric(res6_holo$Ocurrencia)
res6_holo <- separate(res6_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res6_holo <- separate(res6_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:22){
  res6_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res6_abierta$Donor_res[i])] 
  res6_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res6_abierta$Aceptor_res[i])] 
}
for(i in 1:12){
  res6_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res6_holo$Donor_res[i])]
  res6_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res6_holo$Aceptor_res[i])]
}

res6 <- merge(x = res6_abierta, y = res6_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res6)[names(res6)=="Ocurrencia.x"] <- "Abierta"
names(res6)[names(res6)=="Ocurrencia.y"] <- "Holo"

res6$Abierta[which(is.na(res6$Abierta))] = 0
res6$Holo[which(is.na(res6$Holo))] = 0

res6$abMholo = res6$Abierta - res6$Holo
res6 <- res6[order(abs(res6$abMholo), decreasing = TRUE),]

res6ab <- subset(res6, abMholo > 0)
subset(res6ab, abMholo > 25)

res6ab_holo <- res6ab
res6ab_holo$Donor_res <- as.numeric(res6ab_holo$Donor_res)
res6ab_holo$Aceptor_res <- as.numeric(res6ab_holo$Aceptor_res)

for(i in 1:21){
  res6ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res6ab_holo$Donor_res[i])]
  res6ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res6ab_holo$Aceptor_res[i])]
}
subset(res6ab_holo, abMholo > 25)
residues7_holo <- c("353", "361")
residues7_holo <- as.numeric(residues7_holo)
residues7_holo <- as.data.frame(residues7_holo)
write.csv(residues7_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res7.txt", row.names = FALSE)

res6ab_abierta <- res6ab
res6ab_abierta$Donor_res <- as.numeric(res6ab_abierta$Donor_res)
res6ab_abierta$Aceptor_res <- as.numeric(res6ab_abierta$Aceptor_res)

for(i in 1:21){
  res6ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res6ab_abierta$Donor_res[i])]
  res6ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res6ab_abierta$Aceptor_res[i])]
}
subset(res6ab_abierta, abMholo > 25)

residues7_abierta <- c("156", "361")
residues7_abierta <- as.numeric(residues7_abierta)
residues7_abierta <- as.data.frame(residues7_abierta)
write.csv(residues7_abierta, "6nt6_NO_BORRAR/hbonds/res7.txt", row.names = FALSE)


############################## RES 7 ##############################
res7_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res7.csv")
res7_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res7.csv")

res7_abierta <- res7_abierta[-1]
names(res7_abierta)[names(res7_abierta)=="Ocupación...."] <- "Ocurrencia"
res7_abierta[] <- lapply(res7_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res7_abierta[] <- lapply(res7_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res7_abierta$Ocurrencia <- as.numeric(res7_abierta$Ocurrencia)
res7_abierta <- separate(res7_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res7_abierta <- separate(res7_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res7_holo <- res7_holo[-1]
names(res7_holo)[names(res7_holo)=="Ocupación...."] <- "Ocurrencia"
res7_holo[] <- lapply(res7_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res7_holo[] <- lapply(res7_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res7_holo$Ocurrencia <- as.numeric(res7_holo$Ocurrencia)
res7_holo <- separate(res7_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res7_holo <- separate(res7_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:12){
  res7_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res7_abierta$Donor_res[i])] 
  res7_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res7_abierta$Aceptor_res[i])] 
}
for(i in 1:19){
  res7_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res7_holo$Donor_res[i])]
  res7_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res7_holo$Aceptor_res[i])]
}

res7 <- merge(x = res7_abierta, y = res7_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res7)[names(res7)=="Ocurrencia.x"] <- "Abierta"
names(res7)[names(res7)=="Ocurrencia.y"] <- "Holo"

res7$Abierta[which(is.na(res7$Abierta))] = 0
res7$Holo[which(is.na(res7$Holo))] = 0

res7$abMholo = res7$Abierta - res7$Holo
res7 <- res7[order(abs(res7$abMholo), decreasing = TRUE),]

res7ab <- subset(res7, abMholo > 0)
subset(res7ab, abMholo > 25)

res7ab_holo <- res7ab
res7ab_holo$Donor_res <- as.numeric(res7ab_holo$Donor_res)
res7ab_holo$Aceptor_res <- as.numeric(res7ab_holo$Aceptor_res)

for(i in 1:9){
  res7ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res7ab_holo$Donor_res[i])]
  res7ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res7ab_holo$Aceptor_res[i])]
}
subset(res7ab_holo, abMholo > 25)
residues8_holo <- c("349", "356")
residues8_holo <- as.numeric(residues8_holo)
residues8_holo <- as.data.frame(residues8_holo)
write.csv(residues8_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res8.txt", row.names = FALSE)

res7ab_abierta <- res7ab
res7ab_abierta$Donor_res <- as.numeric(res7ab_abierta$Donor_res)
res7ab_abierta$Aceptor_res <- as.numeric(res7ab_abierta$Aceptor_res)

for(i in 1:9){
  res7ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res7ab_abierta$Donor_res[i])]
  res7ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res7ab_abierta$Aceptor_res[i])]
}
subset(res7ab_abierta, abMholo > 25)

residues8_abierta <- c("152", "159")
residues8_abierta <- as.numeric(residues8_abierta)
residues8_abierta <- as.data.frame(residues8_abierta)
write.csv(residues8_abierta, "6nt6_NO_BORRAR/hbonds/res8.txt", row.names = FALSE)


############################## RES 8 ##############################
res8_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res8.csv")
res8_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res8.csv")

res8_abierta <- res8_abierta[-1]
names(res8_abierta)[names(res8_abierta)=="Ocupación...."] <- "Ocurrencia"
res8_abierta[] <- lapply(res8_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res8_abierta[] <- lapply(res8_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res8_abierta$Ocurrencia <- as.numeric(res8_abierta$Ocurrencia)
res8_abierta <- separate(res8_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res8_abierta <- separate(res8_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res8_holo <- res8_holo[-1]
names(res8_holo)[names(res8_holo)=="Ocupación...."] <- "Ocurrencia"
res8_holo[] <- lapply(res8_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res8_holo[] <- lapply(res8_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res8_holo$Ocurrencia <- as.numeric(res8_holo$Ocurrencia)
res8_holo <- separate(res8_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res8_holo <- separate(res8_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:26){
  res8_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res8_abierta$Donor_res[i])] 
  res8_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res8_abierta$Aceptor_res[i])] 
}
for(i in 1:19){
  res8_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res8_holo$Donor_res[i])]
  res8_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res8_holo$Aceptor_res[i])]
}

res8 <- merge(x = res8_abierta, y = res8_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res8)[names(res8)=="Ocurrencia.x"] <- "Abierta"
names(res8)[names(res8)=="Ocurrencia.y"] <- "Holo"

res8$Abierta[which(is.na(res8$Abierta))] = 0
res8$Holo[which(is.na(res8$Holo))] = 0

res8$abMholo = res8$Abierta - res8$Holo
res8 <- res8[order(abs(res8$abMholo), decreasing = TRUE),]

res8ab <- subset(res8, abMholo > 0)
subset(res8ab, abMholo > 25)

res8ab_holo <- res8ab
res8ab_holo$Donor_res <- as.numeric(res8ab_holo$Donor_res)
res8ab_holo$Aceptor_res <- as.numeric(res8ab_holo$Aceptor_res)

for(i in 1:23){
  res8ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res8ab_holo$Donor_res[i])]
  res8ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res8ab_holo$Aceptor_res[i])]
}
subset(res8ab_holo, abMholo > 25)
residues9_holo <- c("292", "345")
residues9_holo <- as.numeric(residues9_holo)
residues9_holo <- as.data.frame(residues9_holo)
write.csv(residues9_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res9.txt", row.names = FALSE)

res8ab_abierta <- res8ab
res8ab_abierta$Donor_res <- as.numeric(res8ab_abierta$Donor_res)
res8ab_abierta$Aceptor_res <- as.numeric(res8ab_abierta$Aceptor_res)

for(i in 1:23){
  res8ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res8ab_abierta$Donor_res[i])]
  res8ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res8ab_abierta$Aceptor_res[i])]
}
subset(res8ab_abierta, abMholo > 25)

residues9_abierta <- c("292", "148")
residues9_abierta <- as.numeric(residues9_abierta)
residues9_abierta <- as.data.frame(residues9_abierta)
write.csv(residues9_abierta, "6nt6_NO_BORRAR/hbonds/res9.txt", row.names = FALSE)


############################## RES 9 ##############################
res9_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res9.csv")
res9_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res9.csv")

res9_abierta <- res9_abierta[-1]
names(res9_abierta)[names(res9_abierta)=="Ocupación...."] <- "Ocurrencia"
res9_abierta[] <- lapply(res9_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res9_abierta[] <- lapply(res9_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res9_abierta$Ocurrencia <- as.numeric(res9_abierta$Ocurrencia)
res9_abierta <- separate(res9_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res9_abierta <- separate(res9_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res9_holo <- res9_holo[-1]
names(res9_holo)[names(res9_holo)=="Ocupación...."] <- "Ocurrencia"
res9_holo[] <- lapply(res9_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res9_holo[] <- lapply(res9_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res9_holo$Ocurrencia <- as.numeric(res9_holo$Ocurrencia)
res9_holo <- separate(res9_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res9_holo <- separate(res9_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:15){
  res9_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res9_abierta$Donor_res[i])] 
  res9_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res9_abierta$Aceptor_res[i])] 
}
for(i in 1:12){
  res9_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res9_holo$Donor_res[i])]
  res9_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res9_holo$Aceptor_res[i])]
}

res9 <- merge(x = res9_abierta, y = res9_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res9)[names(res9)=="Ocurrencia.x"] <- "Abierta"
names(res9)[names(res9)=="Ocurrencia.y"] <- "Holo"

res9$Abierta[which(is.na(res9$Abierta))] = 0
res9$Holo[which(is.na(res9$Holo))] = 0

res9$abMholo = res9$Abierta - res9$Holo
res9 <- res9[order(abs(res9$abMholo), decreasing = TRUE),]

res9ab <- subset(res9, abMholo > 0)
subset(res9ab, abMholo > 25)

res9ab_holo <- res9ab
res9ab_holo$Donor_res <- as.numeric(res9ab_holo$Donor_res)
res9ab_holo$Aceptor_res <- as.numeric(res9ab_holo$Aceptor_res)

for(i in 1:11){
  res9ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res9ab_holo$Donor_res[i])]
  res9ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res9ab_holo$Aceptor_res[i])]
}
subset(res9ab_holo, abMholo > 25)
residues10_holo <- c("357", "358")
residues10_holo <- as.numeric(residues10_holo)
residues10_holo <- as.data.frame(residues10_holo)
write.csv(residues10_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res10.txt", row.names = FALSE)

res9ab_abierta <- res9ab
res9ab_abierta$Donor_res <- as.numeric(res9ab_abierta$Donor_res)
res9ab_abierta$Aceptor_res <- as.numeric(res9ab_abierta$Aceptor_res)

for(i in 1:11){
  res9ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res9ab_abierta$Donor_res[i])]
  res9ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res9ab_abierta$Aceptor_res[i])]
}
subset(res9ab_abierta, abMholo > 25)

residues10_abierta <- c("160", "161")
residues10_abierta <- as.numeric(residues10_abierta)
residues10_abierta <- as.data.frame(residues10_abierta)
write.csv(residues10_abierta, "6nt6_NO_BORRAR/hbonds/res10.txt", row.names = FALSE)

############################## RES 10 #############################
res10_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res10.csv")
res10_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res10.csv")

res10_abierta <- res10_abierta[-1]
names(res10_abierta)[names(res10_abierta)=="Ocupación...."] <- "Ocurrencia"
res10_abierta[] <- lapply(res10_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res10_abierta[] <- lapply(res10_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res10_abierta$Ocurrencia <- as.numeric(res10_abierta$Ocurrencia)
res10_abierta <- separate(res10_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res10_abierta <- separate(res10_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res10_holo <- res10_holo[-1]
names(res10_holo)[names(res10_holo)=="Ocupación...."] <- "Ocurrencia"
res10_holo[] <- lapply(res10_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res10_holo[] <- lapply(res10_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res10_holo$Ocurrencia <- as.numeric(res10_holo$Ocurrencia)
res10_holo <- separate(res10_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res10_holo <- separate(res10_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:9){
  res10_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res10_abierta$Donor_res[i])] 
  res10_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res10_abierta$Aceptor_res[i])] 
}
for(i in 1:5){
  res10_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res10_holo$Donor_res[i])]
  res10_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res10_holo$Aceptor_res[i])]
}

res10 <- merge(x = res10_abierta, y = res10_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res10)[names(res10)=="Ocurrencia.x"] <- "Abierta"
names(res10)[names(res10)=="Ocurrencia.y"] <- "Holo"

res10$Abierta[which(is.na(res10$Abierta))] = 0
res10$Holo[which(is.na(res10$Holo))] = 0

res10$abMholo = res10$Abierta - res10$Holo
res10 <- res10[order(abs(res10$abMholo), decreasing = TRUE),]

res10ab <- subset(res10, abMholo > 0)
subset(res10ab, abMholo > 25)

res10ab_holo <- res10ab
res10ab_holo$Donor_res <- as.numeric(res10ab_holo$Donor_res)
res10ab_holo$Aceptor_res <- as.numeric(res10ab_holo$Aceptor_res)

for(i in 1:8){
  res10ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res10ab_holo$Donor_res[i])]
  res10ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res10ab_holo$Aceptor_res[i])]
}
subset(res10ab_holo, abMholo > 25)
residues11_holo <- c("164", "165")
residues11_holo <- as.numeric(residues11_holo)
residues11_holo <- as.data.frame(residues11_holo)
write.csv(residues11_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res11.txt", row.names = FALSE)

res10ab_abierta <- res10ab
res10ab_abierta$Donor_res <- as.numeric(res10ab_abierta$Donor_res)
res10ab_abierta$Aceptor_res <- as.numeric(res10ab_abierta$Aceptor_res)

for(i in 1:8){
  res10ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res10ab_abierta$Donor_res[i])]
  res10ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res10ab_abierta$Aceptor_res[i])]
}
subset(res10ab_abierta, abMholo > 25)

residues11_abierta <- c("164", "165")
residues11_abierta <- as.numeric(residues11_abierta)
residues11_abierta <- as.data.frame(residues11_abierta)
write.csv(residues11_abierta, "6nt6_NO_BORRAR/hbonds/res11.txt", row.names = FALSE)

############################## RES 11 #############################
res11_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res11.csv")
res11_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res11.csv")

res11_abierta <- res11_abierta[-1]
names(res11_abierta)[names(res11_abierta)=="Ocupación...."] <- "Ocurrencia"
res11_abierta[] <- lapply(res11_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res11_abierta[] <- lapply(res11_abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
res11_abierta$Ocurrencia <- as.numeric(res11_abierta$Ocurrencia)
res11_abierta <- separate(res11_abierta, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res11_abierta <- separate(res11_abierta, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)
res11_holo <- res11_holo[-1]
names(res11_holo)[names(res11_holo)=="Ocupación...."] <- "Ocurrencia"
res11_holo[] <- lapply(res11_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res11_holo[] <- lapply(res11_holo, gsub, pattern = " ", replacement = "", fixed = TRUE)
res11_holo$Ocurrencia <- as.numeric(res11_holo$Ocurrencia)
res11_holo <- separate(res11_holo, Aceptor, c("Aceptor_res", "Aceptor_part"), sep = "-", remove = TRUE)
res11_holo <- separate(res11_holo, Donor, c("Donor_res", "Donor_part"), sep = "-", remove = TRUE)

for(i in 1:8){
  res11_abierta$Donor_res[i] <- index$index[which(index$X6nt6 == res11_abierta$Donor_res[i])] 
  res11_abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == res11_abierta$Aceptor_res[i])] 
}
for(i in 1:8){
  res11_holo$Donor_res[i] <- index$index[which(index$X6nt7_holo == res11_holo$Donor_res[i])]
  res11_holo$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == res11_holo$Aceptor_res[i])]
}

res11 <- merge(x = res11_abierta, y = res11_holo, by = c("Donor_res", "Donor_part", "Aceptor_res", "Aceptor_part"), all = TRUE)
names(res11)[names(res11)=="Ocurrencia.x"] <- "Abierta"
names(res11)[names(res11)=="Ocurrencia.y"] <- "Holo"

res11$Abierta[which(is.na(res11$Abierta))] = 0
res11$Holo[which(is.na(res11$Holo))] = 0

res11$abMholo = res11$Abierta - res11$Holo
res11 <- res11[order(abs(res11$abMholo), decreasing = TRUE),]

res11ab <- subset(res11, abMholo > 0)
subset(res11ab, abMholo > 25)

res11ab_holo <- res11ab
res11ab_holo$Donor_res <- as.numeric(res11ab_holo$Donor_res)
res11ab_holo$Aceptor_res <- as.numeric(res11ab_holo$Aceptor_res)

for(i in 1:5){
  res11ab_holo$Donor_res[i] <- index$X6nt7_holo[as.numeric(res11ab_holo$Donor_res[i])]
  res11ab_holo$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(res11ab_holo$Aceptor_res[i])]
}
subset(res11ab_holo, abMholo > 25)
#residues12_holo <- c("308", "366", "371", "472", "473")
#residues12_holo <- as.numeric(residues12_holo)
#residues12_holo <- as.data.frame(residues12_holo)
#write.csv(residues12_holo, "6nt7_HOLO_NO_BORRAR/hbonds/res12.txt", row.names = FALSE)

res11ab_abierta <- res11ab
res11ab_abierta$Donor_res <- as.numeric(res11ab_abierta$Donor_res)
res11ab_abierta$Aceptor_res <- as.numeric(res11ab_abierta$Aceptor_res)

for(i in 1:5){
  res11ab_abierta$Donor_res[i] <- index$X6nt6[as.numeric(res11ab_abierta$Donor_res[i])]
  res11ab_abierta$Aceptor_res[i] <- index$X6nt6[as.numeric(res11ab_abierta$Aceptor_res[i])]
}
subset(res11ab_abierta, abMholo > 25)

#residues12_abierta <- c("308", "366", "371", "472", "473")
#residues12_abierta <- as.numeric(residues12_abierta)
#residues12_abierta <- as.data.frame(residues12_abierta)
#write.csv(residues12_abierta, "6nt6_NO_BORRAR/hbonds/res12.txt", row.names = FALSE)