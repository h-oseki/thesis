setwd("~/Desktop/Sysbio/cSTING")

############################## RES 1 ##############################
###                                                             ###
### Residuos:                                                   ###
###   - ASP280A                                                 ###
###   - CYS281A                                                 ###
###   - ALA282A                                                 ###
###   - ALA283A                                                 ###
###   - PHE284A                                                 ###
###   - SER285A                                                 ###
###   - ARG286A                                                 ###
###   - GLU287A                                                 ###
###   - ASP280B                                                 ###
###   - CYS271B                                                 ###
###   - ALA282B                                                 ###
###   - ALA283B                                                 ###
###   - PHE284B                                                 ###
###   - SER285B                                                 ###
###   - ARG286B                                                 ###
###   - GLU287B                                                 ###
###################################################################

res1_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res1.csv")
res1_apo <- read.csv("6nt7_APO_NO_BORRAR/hbonds/hbonds_res1.csv")
res1_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res1.csv")

res1_abierta <- res1_abierta[-1]
names(res1_abierta)[names(res1_abierta)=="Ocupación...."] <- "Ocurrencia"
res1_abierta[] <- lapply(res1_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res1_abierta$Ocurrencia <- as.numeric(res1_abierta$Ocurrencia)
res1_apo <- res1_apo[-1]
names(res1_apo)[names(res1_apo)=="Ocupación...."] <- "Ocurrencia"
res1_apo[] <- lapply(res1_apo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res1_apo$Ocurrencia <- as.numeric(res1_apo$Ocurrencia)
res1_holo <- res1_holo[-1]
names(res1_holo)[names(res1_holo)=="Ocupación...."] <- "Ocurrencia"
res1_holo[] <- lapply(res1_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res1_holo$Ocurrencia <- as.numeric(res1_holo$Ocurrencia)

res1 <- merge(x = res1_abierta, y = res1_apo, by = c("Donor", "Aceptor"), all = TRUE)
names(res1)[names(res1)=="Ocurrencia.x"] <- "Abierta"
names(res1)[names(res1)=="Ocurrencia.y"] <- "Apo"
res1 <- merge(x = res1, y = res1_holo, by = c("Donor", "Aceptor"), all = TRUE)
names(res1)[names(res1)=="Ocurrencia"] <- "Holo"

res1$Abierta[which(is.na(res1$Abierta))] = 0
res1$Apo[which(is.na(res1$Apo))] = 0
res1$Holo[which(is.na(res1$Holo))] = 0

res1$abMholo = res1$Abierta - res1$Holo
res1 <- res1[order(abs(res1$abMholo), decreasing = TRUE),]

res1ab <- subset(res1, abMholo > 0)
subset(res1ab, abMholo > 30)

residues1 <- read.csv("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds/res1.txt", sep="")
residues2 <- c("279", "473")
residues2 <- as.numeric(residues2)
residues2 <- as.data.frame(residues2)
write.csv(residues2, "6nt6_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)
write.csv(residues2, "6nt7_HOLO_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)
write.csv(residues2, "6nt7_APO_NO_BORRAR/hbonds/res2.txt", row.names = FALSE)

############################## RES 2 ##############################
###                                                             ###
### Residuos:                                                   ###
###   - CYS281A:ASP279A                                         ###
###   - ALA282A:ASP279A                                         ###
###   - ALA282B:MET276B                                         ###
###################################################################

res2_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res2.csv")
res2_apo <- read.csv("6nt7_APO_NO_BORRAR/hbonds/hbonds_res2.csv")
res2_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res2.csv")

res2_abierta <- res2_abierta[-1]
names(res2_abierta)[names(res2_abierta)=="Ocupación...."] <- "Ocurrencia"
res2_abierta[] <- lapply(res2_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res2_abierta$Ocurrencia <- as.numeric(res2_abierta$Ocurrencia)
res2_apo <- res2_apo[-1]
names(res2_apo)[names(res2_apo)=="Ocupación...."] <- "Ocurrencia"
res2_apo[] <- lapply(res2_apo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res2_apo$Ocurrencia <- as.numeric(res2_apo$Ocurrencia)
res2_holo <- res2_holo[-1]
names(res2_holo)[names(res2_holo)=="Ocupación...."] <- "Ocurrencia"
res2_holo[] <- lapply(res2_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res2_holo$Ocurrencia <- as.numeric(res2_holo$Ocurrencia)

res2 <- merge(x = res2_abierta, y = res2_apo, by = c("Donor", "Aceptor"), all = TRUE)
names(res2)[names(res2)=="Ocurrencia.x"] <- "Abierta"
names(res2)[names(res2)=="Ocurrencia.y"] <- "Apo"
res2 <- merge(x = res2, y = res2_holo, by = c("Donor", "Aceptor"), all = TRUE)
names(res2)[names(res2)=="Ocurrencia"] <- "Holo"

res2$Abierta[which(is.na(res2$Abierta))] = 0
res2$Apo[which(is.na(res2$Apo))] = 0
res2$Holo[which(is.na(res2$Holo))] = 0

res2$abMholo = res2$Abierta - res2$Holo
res2 <- res2[order(abs(res2$abMholo), decreasing = TRUE),]

res2ab <- subset(res2, abMholo > 0)
subset(res2ab, abMholo > 30)

residues2 <- read.csv("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds/res2.txt", sep="")
residues2 <- c("366", "371")
residues2 <- as.numeric(residues2)
residues2 <- as.data.frame(residues2)
write.csv(residues2, "6nt6_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)
write.csv(residues2, "6nt7_HOLO_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)
write.csv(residues2, "6nt7_APO_NO_BORRAR/hbonds/res3.txt", row.names = FALSE)

############################## RES 3 ##############################
###                                                             ###
### Residuos:                                                   ###
###   - CYS281A:ASP279A:TYR169B                                 ###
###   - ALA282A:ASP279A:LYS174B                                 ###
###################################################################

res3_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res3.csv")
res3_apo <- read.csv("6nt7_APO_NO_BORRAR/hbonds/hbonds_res3.csv")
res3_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res3.csv")

res3_abierta <- res3_abierta[-1]
names(res3_abierta)[names(res3_abierta)=="Ocupación...."] <- "Ocurrencia"
res3_abierta[] <- lapply(res3_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res3_abierta$Ocurrencia <- as.numeric(res3_abierta$Ocurrencia)
res3_apo <- res3_apo[-1]
names(res3_apo)[names(res3_apo)=="Ocupación...."] <- "Ocurrencia"
res3_apo[] <- lapply(res3_apo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res3_apo$Ocurrencia <- as.numeric(res3_apo$Ocurrencia)
res3_holo <- res3_holo[-1]
names(res3_holo)[names(res3_holo)=="Ocupación...."] <- "Ocurrencia"
res3_holo[] <- lapply(res3_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res3_holo$Ocurrencia <- as.numeric(res3_holo$Ocurrencia)

res3 <- merge(x = res3_abierta, y = res3_apo, by = c("Donor", "Aceptor"), all = TRUE)
names(res3)[names(res3)=="Ocurrencia.x"] <- "Abierta"
names(res3)[names(res3)=="Ocurrencia.y"] <- "Apo"
res3 <- merge(x = res3, y = res3_holo, by = c("Donor", "Aceptor"), all = TRUE)
names(res3)[names(res3)=="Ocurrencia"] <- "Holo"

res3$Abierta[which(is.na(res3$Abierta))] = 0
res3$Apo[which(is.na(res3$Apo))] = 0
res3$Holo[which(is.na(res3$Holo))] = 0

res3$abMholo = res3$Abierta - res3$Holo
res3 <- res3[order(abs(res3$abMholo), decreasing = TRUE),]

res3ab <- subset(res3, abMholo > 0)
subset(res3ab, abMholo > 30)

residues3 <- read.csv("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/hbonds/res3.txt", sep="")
residues4 <- c("362", "504")
residues4 <- as.numeric(residues4)
residues4 <- as.data.frame(residues4)
write.csv(residues4, "6nt6_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)
write.csv(residues4, "6nt7_HOLO_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)
write.csv(residues4, "6nt7_APO_NO_BORRAR/hbonds/res4.txt", row.names = FALSE)

############################## RES 4 ##############################
###                                                             ###
### Residuos:                                                   ###
###   - CYS281A:ASP279A:TYR169B:ALA165B                         ###
###   - ALA282A:ASP279A:LYS174B:SER307B                         ###
###################################################################

res4_abierta <- read.csv("6nt6_NO_BORRAR/hbonds/hbonds_res4.csv")
res4_apo <- read.csv("6nt7_APO_NO_BORRAR/hbonds/hbonds_res4.csv")
res4_holo <- read.csv("6nt7_HOLO_NO_BORRAR/hbonds/hbonds_res4.csv")

res4_abierta <- res4_abierta[-1]
names(res4_abierta)[names(res4_abierta)=="Ocupación...."] <- "Ocurrencia"
res4_abierta[] <- lapply(res4_abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
res4_abierta$Ocurrencia <- as.numeric(res4_abierta$Ocurrencia)
res4_apo <- res4_apo[-1]
names(res4_apo)[names(res4_apo)=="Ocupación...."] <- "Ocurrencia"
res4_apo[] <- lapply(res4_apo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res4_apo$Ocurrencia <- as.numeric(res4_apo$Ocurrencia)
res4_holo <- res4_holo[-1]
names(res4_holo)[names(res4_holo)=="Ocupación...."] <- "Ocurrencia"
res4_holo[] <- lapply(res4_holo, gsub, pattern = "%", replacement = "", fixed = TRUE)
res4_holo$Ocurrencia <- as.numeric(res4_holo$Ocurrencia)

res4 <- merge(x = res4_abierta, y = res4_apo, by = c("Donor", "Aceptor"), all = TRUE)
names(res4)[names(res4)=="Ocurrencia.x"] <- "Abierta"
names(res4)[names(res4)=="Ocurrencia.y"] <- "Apo"
res4 <- merge(x = res4, y = res4_holo, by = c("Donor", "Aceptor"), all = TRUE)
names(res4)[names(res4)=="Ocurrencia"] <- "Holo"

res4$Abierta[which(is.na(res4$Abierta))] = 0
res4$Apo[which(is.na(res4$Apo))] = 0
res4$Holo[which(is.na(res4$Holo))] = 0

res4$abMholo = res4$Abierta - res4$Holo
res4 <- res4[order(abs(res4$abMholo), decreasing = TRUE),]

res4ab <- subset(res4, abMholo > 0)
subset(res4ab, abMholo > 30)

residues5 <- c("362", "504")
residues5 <- as.numeric(residues5)
residues5 <- as.data.frame(residues5)
write.csv(residues5, "6nt6_NO_BORRAR/hbonds/res5.txt", row.names = FALSE)
write.csv(residues5, "6nt7_HOLO_NO_BORRAR/hbonds/res5.txt", row.names = FALSE)
write.csv(residues5, "6nt7_APO_NO_BORRAR/hbonds/res5.txt", row.names = FALSE)
