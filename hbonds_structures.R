index <- read.delim("~/Desktop/Sysbio/cSTING/corridas restringidas/index.txt")
index$index <- as.numeric(index$index)
index$X6nt6 <- as.character(index$X6nt6)
index$X6nt7_holo <- as.character(index$X6nt7_holo)

index2 <- data.frame("num" = 1:394,
                     "res" = c("ALA146","VAL147","GLU148","VAL149","SER150","GLU151","LEU152","THR153",
                               "GLU154","SER155","SER156","LYS157","LYS158","ASN159","VAL160","ALA161",
                               "HIS162","GLY163","LEU164","ALA165","TRP166","SER167","TYR168","TYR169",
                               "ILE170","GLY171","TYR172","LEU173","LYS174","VAL175","VAL176","LEU177",
                               "PRO178","ARG179","LEU180","LYS181","GLU182","CYS183","MET184","GLU185",
                               "GLU186","LEU187","SER188","ARG189","THR190","ASN191","PRO192","MET193",
                               "LEU194","ARG195","ALA196","HIS197","ARG198","ASP199","THR200","TRP201",
                               "LYS202","LEU203","HIS204","ILE205","LEU206","VAL207","PRO208","LEU209",
                               "GLY210","CYS211","ASP212","ILE213","TRP214","ASP215","ASP216","LEU217",
                               "GLU218","LYS219","ALA220","ASP221","SER222","ASN223","ILE224","GLN225",
                               "TYR226","LEU227","ALA228","ASP229","LEU230","PRO231","GLU232","THR233",
                               "ILE234","LEU235","THR236","ARG237","ALA238","GLY239","ILE240","LYS241",
                               "ARG242","ARG243","VAL244","TYR245","LYS246","HIS247","SER248","LEU249",
                               "TYR250","VAL251","ILE252","ARG253","ASP254","LYS255","ASP256","ASN257",
                               "LYS258","LEU259","ARG260","PRO261","CYS262","VAL263","LEU264","GLU265",
                               "PHE266","ALA267","SER268","PRO269","LEU270","GLN271","THR272","LEU273",
                               "CYS274","ALA275","MET276","SER277","GLN278","ASP279","ASP280","CYS281",
                               "ALA282","ALA283","PHE284","SER285","ARG286","GLU287","GLN288","ARG289",
                               "LEU290","GLU291","GLN292","ALA293","ARG294","LEU295","PHE296","TYR297",
                               "ARG298","SER299","LEU300","ARG301","ASP302","ILE303","LEU304","GLY305",
                               "SER306","SER307","LYS308","GLU309","CYS310","ALA311","GLY312","LEU313",
                               "TYR314","ARG315","LEU316","ILE317","ALA318","TYR319","GLU320","GLU321",
                               "PRO322","ALA323","GLU324","PRO325","GLU326","SER327","HIS328","PHE329",
                               "LEU330","SER331","GLY332","LEU333","ILE334","LEU335","TRP336","HIS337",
                               "LEU338","GLN339","GLN340","GLN341","GLN342","ALA343","VAL344","GLU345",
                               "VAL346","SER347","GLU348","LEU349","THR350","GLU351","SER352","SER353",
                               "LYS354","LYS355","ASN356","VAL357","ALA358","HIS359","GLY360","LEU361",
                               "ALA362","TRP363","SER364","TYR365","TYR366","ILE367","GLY368","TYR369",
                               "LEU370","LYS371","VAL372","VAL373","LEU374","PRO375","ARG376","LEU377",
                               "LYS378","GLU379","CYS380","MET381","GLU382","GLU383","LEU384","SER385",
                               "ARG386","THR387","ASN388","PRO389","MET390","LEU391","ARG392","ALA393",
                               "HIS394","ARG395","ASP396","THR397","TRP398","LYS399","LEU400","HIS401",
                               "ILE402","LEU403","VAL404","PRO405","LEU406","GLY407","CYS408","ASP409",
                               "ILE410","TRP411","ASP412","ASP413","LEU414","GLU415","LYS416","ALA417",
                               "ASP418","SER419","ASN420","ILE421","GLN422","TYR423","LEU424","ALA425",
                               "ASP426","LEU427","PRO428","GLU429","THR430","ILE431","LEU432","THR433",
                               "ARG434","ALA435","GLY436","ILE437","LYS438","ARG439","ARG440","VAL441",
                               "TYR442","LYS443","HIS444","SER445","LEU446","TYR447","VAL448","ILE449",
                               "ARG450","ASP451","LYS452","ASP453","ASN454","LYS455","LEU456","ARG457",
                               "PRO458","CYS459","VAL460","LEU461","GLU462","PHE463","ALA464","SER465",
                               "PRO466","LEU467","GLN468","THR469","LEU470","CYS471","ALA472","MET473",
                               "SER474","GLN475","ASP476","ASP477","CYS478","ALA479","ALA480","PHE481",
                               "SER482","ARG483","GLU484","GLN485","ARG486","LEU487","GLU488","GLN489",
                               "ALA490","ARG491","LEU492","PHE493","TYR494","ARG495","SER496","LEU497",
                               "ARG498","ASP499","ILE500","LEU501","GLY502","SER503","SER504","LYS505",
                               "GLU506","CYS507","ALA508","GLY509","LEU510","TYR511","ARG512","LEU513",
                               "ILE514","ALA515","TYR516","GLU517","GLU518","PRO519","ALA520","GLU521",
                               "PRO522","GLU523","SER524","HIS525","PHE526","LEU527","SER528","GLY529",
                               "LEU530","ILE531","LEU532","TRP533","HIS534","LEU535","GLN536","GLN537",
                               "GLN538","GLN539"))

cerrada <- read.delim("~/Desktop/Sysbio/cSTING/hbonds/hbonds-cerrada-details.dat")
abierta <- read.delim("~/Desktop/Sysbio/cSTING/hbonds/hbonds-details-abierta.dat")
rotada <- read.delim("~/Desktop/Sysbio/cSTING/hbonds/hbonds-details-rotada.dat")

abierta[] <- lapply(abierta, gsub, pattern = " ", replacement = "", fixed = TRUE)
cerrada[] <- lapply(cerrada, gsub, pattern = " ", replacement = "", fixed = TRUE)
rotada[] <- lapply(rotada, gsub, pattern = " ", replacement = "", fixed = TRUE)

abierta[] <- lapply(abierta, gsub, pattern = "%", replacement = "", fixed = TRUE)
cerrada[] <- lapply(cerrada, gsub, pattern = "%", replacement = "", fixed = TRUE)
rotada[] <- lapply(rotada, gsub, pattern = "%", replacement = "", fixed = TRUE)

abierta$occupancy <- as.numeric(abierta$occupancy)
cerrada$occupancy <- as.numeric(cerrada$occupancy)
rotada$occupancy <- as.numeric(rotada$occupancy)

abierta <- separate(abierta, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
abierta <- separate(abierta, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
cerrada <- separate(cerrada, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
cerrada <- separate(cerrada, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)
rotada <- separate(rotada, acceptor, c("Aceptor_res", "Aceptor_part", "Aceptor_atom"), sep = "-", remove = TRUE)
rotada <- separate(rotada, donor, c("Donor_res", "Donor_part", "Donor_atom"), sep = "-", remove = TRUE)

for(i in 1:2301){
  abierta$Donor_res[i] <- index$index[which(index$X6nt6 == abierta$Donor_res[i])] 
  abierta$Aceptor_res[i] <- index$index[which(index$X6nt6 == abierta$Aceptor_res[i])] 
}
for(i in 1:2197){
  cerrada$Donor_res[i] <- index$index[which(index$X6nt7_holo == cerrada$Donor_res[i])]
  cerrada$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == cerrada$Aceptor_res[i])]
}
for(i in 1:7094){
  rotada$Donor_res[i] <- index$index[which(index$X6nt7_holo == rotada$Donor_res[i])]
  rotada$Aceptor_res[i] <- index$index[which(index$X6nt7_holo == rotada$Aceptor_res[i])]
}

hbonds <- merge(x = abierta, y = cerrada, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy.x"] <- "Abierta"
names(hbonds)[names(hbonds)=="occupancy.y"] <- "Cerrada"
hbonds <- merge(x = hbonds, y = rotada, by = c("Donor_res", "Donor_part", "Donor_atom", "Aceptor_res", "Aceptor_part", "Aceptor_atom"), all = TRUE)
names(hbonds)[names(hbonds)=="occupancy"] <- "Rotada"


hbonds$Abierta[which(is.na(hbonds$Abierta))] = 0
hbonds$Cerrada[which(is.na(hbonds$Cerrada))] = 0
hbonds$Rotada[which(is.na(hbonds$Rotada))] = 0

hbonds$CmA <- hbonds$Cerrada - hbonds$Abierta
hbonds$RmA <- hbonds$Rotada - hbonds$Abierta
hbonds$CmR <- hbonds$Cerrada - hbonds$Rotada


#### Diferencia entre Cerrada y Abierta, con ocurrencia diferencial en CERRADA ####
CmA <- hbonds[which(hbonds$CmA > 20),c(1:8,10)]
CmA$Donor_res <- as.numeric(CmA$Donor_res)
CmA$Aceptor_res <- as.numeric(CmA$Aceptor_res)
for(i in 1:167){
  CmA$Donor_res[i] <- index$X6nt7_holo[as.numeric(CmA$Donor_res[i])]
  CmA$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(CmA$Aceptor_res[i])]
}
names(CmA)[names(CmA)=="CmA"] <- "Dif"
write.csv(CmA, "hbonds/CmA2.csv", row.names = FALSE)

#
#### Diferencia entre Cerrada y Abierta, con ocurrencia diferencial en ABIERTA #### 
AmC <- hbonds[which(hbonds$CmA < -20),c(1:8,10)]
AmC$Donor_res <- as.numeric(AmC$Donor_res)
AmC$Aceptor_res <- as.numeric(AmC$Aceptor_res)
for(i in 1:136){
  AmC$Donor_res[i] <- index$X6nt6[as.numeric(AmC$Donor_res[i])]
  AmC$Aceptor_res[i] <- index$X6nt6[as.numeric(AmC$Aceptor_res[i])]
}
names(AmC)[names(AmC)=="CmA"] <- "Dif"
AmC$Dif <- abs(AmC$Dif)
write.csv(AmC, "hbonds/AmC2.csv", row.names = FALSE)
#
#### Diferencia entre Cerrada y Rotada, con ocurrencia diferencial en CERRADA #### 
CmR <- hbonds[which(hbonds$CmR > 20),c(1:6,8,9,12)]
CmR$Donor_res <- as.numeric(CmR$Donor_res)
CmR$Aceptor_res <- as.numeric(CmR$Aceptor_res)
for(i in 1:123){
  CmR$Donor_res[i] <- index$X6nt7_holo[as.numeric(CmR$Donor_res[i])]
  CmR$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(CmR$Aceptor_res[i])]
}
names(CmR)[names(CmR)=="CmR"] <- "Dif"
write.csv(CmR, "hbonds/CmR2.csv", row.names = FALSE)
#
#### Diferencia entre Cerrada y Rotada, con ocurrencia diferencial en ROTADA #### 
RmC <- hbonds[which(hbonds$CmR < -20),c(1:6,8,9,12)]
RmC$Donor_res <- as.numeric(RmC$Donor_res)
RmC$Aceptor_res <- as.numeric(RmC$Aceptor_res)
for(i in 1:115){
  RmC$Donor_res[i] <- index$X6nt7_holo[as.numeric(RmC$Donor_res[i])]
  RmC$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(RmC$Aceptor_res[i])]
}
names(RmC)[names(RmC)=="CmR"] <- "Dif"
RmC$Dif <- abs(RmC$Dif)
write.csv(RmC, "hbonds/RmC2.csv", row.names = FALSE)
#### Diferencia entre Rotada y Abierta, con ocurrencia diferencial en ROTADA #### 
RmA <- hbonds[which(hbonds$RmA > 20),c(1:7,9,11)]
RmA$Donor_res <- as.numeric(RmA$Donor_res)
RmA$Aceptor_res <- as.numeric(RmA$Aceptor_res)
for(i in 1:161){
  RmA$Donor_res[i] <- index$X6nt7_holo[as.numeric(RmA$Donor_res[i])]
  RmA$Aceptor_res[i] <- index$X6nt7_holo[as.numeric(RmA$Aceptor_res[i])]
}
names(RmA)[names(RmA)=="RmA"] <- "Dif"
write.csv(RmA, "hbonds/RmA2.csv", row.names = FALSE)
#### Diferencia entre Rotada y Abierta, con ocurrencia diferencial en ABIERTA #### 
AmR <- hbonds[which(hbonds$RmA < -20),c(1:7,9,11)]
AmR$Donor_res <- as.numeric(AmR$Donor_res)
AmR$Aceptor_res <- as.numeric(AmR$Aceptor_res)
for(i in 1:127){
  AmR$Donor_res[i] <- index$X6nt6[as.numeric(AmR$Donor_res[i])]
  AmR$Aceptor_res[i] <- index$X6nt6[as.numeric(AmR$Aceptor_res[i])]
}
names(AmR)[names(AmR)=="RmA"] <- "Dif"
AmR$Dif <- abs(AmR$Dif)
write.csv(AmR, "hbonds/AmR2.csv", row.names = FALSE)
