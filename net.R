setwd("~/Desktop/Sysbio/cSTING")
a <- theme(
  panel.grid.major = element_line(color ='gray90'), 
  panel.grid.minor = element_line(color ='gray95'),
  plot.title = element_text(colour = 'black', size = 15),
  axis.title = element_text(colour = 'black', face = 'italic'),
  panel.background = element_rect(fill = "white", colour = "black")
) 
#### data sets ####
abierta_bet <- read.table("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/network_analysis/abierta-bet.dat", quote="\"", comment.char="")
abierta_net <- read.table("~/Desktop/Sysbio/cSTING/6nt6_NO_BORRAR/network_analysis/abierta-net.dat", quote="\"", comment.char="")
cerrada_bet <- read.table("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/network_analysis/cerrada-bet.dat", quote="\"", comment.char="")
cerrada_net <- read.table("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/network_analysis/cerrada-net.dat", quote="\"", comment.char="")
abierta_res <- read.table("~/Desktop/Sysbio/cSTING/6nt7_HOLO_NO_BORRAR/6nt7_holo_aa.txt", quote="\"", comment.char="")
abierta_res <- c("ALA146A", "VAL147A", "GLU148A", "VAL149A", "SER150A", "GLU151A",
                "LEU152A", "THR153A", "GLU154A", "SER155A", "SER156A", "LYS157A",
                "LYS158A", "ASN159A", "VAL160A", "ALA161A", "HIS162A", "GLY163A",
                "LEU164A", "ALA165A", "TRP166A", "SER167A", "TYR168A", "TYR169A",
                "ILE170A", "GLY171A", "TYR172A", "LEU173A", "LYS174A", "VAL175A",
                "VAL176A", "LEU177A", "PRO178A", "ARG179A", "LEU180A", "LYS181A",
                "GLU182A", "CYS183A", "MET184A", "GLU185A", "GLU186A", "LEU187A",
                "SER188A", "ARG189A", "THR190A", "ASN191A", "PRO192A", "MET193A",
                "LEU194A", "ARG195A", "ALA196A", "HIS197A", "ARG198A", "ASP199A",
                "THR200A", "TRP201A", "LYS202A", "LEU203A", "HIS204A", "ILE205A",
                "LEU206A", "VAL207A", "PRO208A", "LEU209A", "GLY210A", "CYS211A",
                "ASP212A", "ILE213A", "TRP214A", "ASP215A", "ASP216A", "LEU217A",
                "GLU218A", "LYS219A", "ALA220A", "ASP221A", "SER222A", "ASN223A",
                "ILE224A", "GLN225A", "TYR226A", "LEU227A", "ALA228A", "ASP229A",
                "LEU230A", "PRO231A", "GLU232A", "THR233A", "ILE234A", "LEU235A",
                "THR236A", "ARG237A", "ALA238A", "GLY239A", "ILE240A", "LYS241A",
                "ARG242A", "ARG243A", "VAL244A", "TYR245A", "LYS246A", "HIS247A",
                "SER248A", "LEU249A", "TYR250A", "VAL251A", "ILE252A", "ARG253A",
                "ASP254A", "LYS255A", "ASP256A", "ASN257A", "LYS258A", "LEU259A",
                "ARG260A", "PRO261A", "CYS262A", "VAL263A", "LEU264A", "GLU265A",
                "PHE266A", "ALA267A", "SER268A", "PRO269A", "LEU270A", "GLN271A",
                "THR272A", "LEU273A", "CYS274A", "ALA275A", "MET276A", "SER277A", 
                "GLN278A", "ASP279A", "ASP280A", "CYS281A", "ALA282A", "ALA283A",
                "PHE284A", "SER285A", "ARG286A", "GLU287A", "GLN288A", "ARG289A",
                "LEU290A", "GLU291A", "GLN292A", "ALA293A", "ARG294A", "LEU295A",
                "PHE296A", "TYR297A", "ARG298A", "SER299A", "LEU300A", "ARG301A",
                "ASP302A", "ILE303A", "LEU304A", "GLY305A", "SER306A", "SER307A",
                "LYS308A", "GLU309A", "CYS310A", "ALA311A", "GLY312A", "LEU313A",
                "TYR314A", "ARG315A", "LEU316A", "ILE317A", "ALA318A", "TYR319A",
                "GLU320A", "GLU321A", "PRO322A", "ALA323A", "GLU324A", "PRO325A",
                "GLU326A", "SER327A", "HIS328A", "PHE329A", "LEU330A", "SER331A",
                "GLY332A", "LEU333A", "ILE334A", "LEU335A", "TRP336A", "HIS337A",
                "LEU338A", "GLN339A", "GLN340A", "GLN341A", "GLN342A", "ALA146B",
                "VAL147B", "GLU148B", "VAL149B", "SER150B", "GLU151B", "LEU152B",
                "THR153B", "GLU154B", "SER155B", "SER156B", "LYS157B", "LYS158B",
                "ASN159B", "VAL160B", "ALA161B", "HIS162B", "GLY163B", "LEU164B",
                "ALA165B", "TRP166B", "SER167B", "TYR168B", "TYR169B", "ILE170B",
                "GLY171B", "TYR172B", "LEU173B", "LYS174B", "VAL175B", "VAL176B",
                "LEU177B", "PRO178B", "ARG179B", "LEU180B", "LYS181B", "GLU182B",
                "CYS183B", "MET184B", "GLU185B", "GLU186B", "LEU187B", "SER188B",
                "ARG189B", "THR190B", "ASN191B", "PRO192B", "MET193B", "LEU194B",
                "ARG195B", "ALA196B", "HIS197B", "ARG198B", "ASP199B", "THR200B",
                "TRP201B", "LYS202B", "LEU203B", "HIS204B", "ILE205B", "LEU206B",
                "VAL207B", "PRO208B", "LEU209B", "GLY210B", "CYS211B", "ASP212B",
                "ILE213B", "TRP214B", "ASP215B", "ASP216B", "LEU217B", "GLU218B",
                "LYS219B", "ALA220B", "ASP221B", "SER222B", "ASN223B", "ILE224B",
                "GLN225B", "TYR226B", "LEU227B", "ALA228B", "ASP229B", "LEU230B",
                "PRO231B", "GLU232B", "THR233B", "ILE234B", "LEU235B", "THR236B",
                "ARG237B", "ALA238B", "GLY239B", "ILE240B", "LYS241B", "ARG242B",
                "ARG243B", "VAL244B", "TYR245B", "LYS246B", "HIS247B", "SER248B",
                "LEU249B", "TYR250B", "VAL251B", "ILE252B", "ARG253B", "ASP254B",
                "LYS255B", "ASP256B", "ASN257B", "LYS258B", "LEU259B", "ARG260B",
                "PRO261B", "CYS262B", "VAL263B", "LEU264B", "GLU265B", "PHE266B",
                "ALA267B", "SER268B", "PRO269B", "LEU270B", "GLN271B", "THR272B",
                "LEU273B", "CYS274B", "ALA275B", "MET276B", "SER277B", "GLN278B",
                "ASP279B", "ASP280B", "CYS281B", "ALA282B", "ALA283B", "PHE284B",
                "SER285B", "ARG286B", "GLU287B", "GLN288B", "ARG289B", "LEU290B",
                "GLU291B", "GLN292B", "ALA293B", "ARG294B", "LEU295B", "PHE296B",
                "TYR297B", "ARG298B", "SER299B", "LEU300B", "ARG301B", "ASP302B",
                "ILE303B", "LEU304B", "GLY305B", "SER306B", "SER307B", "LYS308B",
                "GLU309B", "CYS310B", "ALA311B", "GLY312B", "LEU313B", "TYR314B",
                "ARG315B", "LEU316B", "ILE317B", "ALA318B", "TYR319B", "GLU320B",
                "GLU321B", "PRO322B", "ALA323B", "GLU324B", "PRO325B", "GLU326B",
                "SER327B", "HIS328B", "PHE329B", "LEU330B", "SER331B", "GLY332B",
                "LEU333B", "ILE334B", "LEU335B", "TRP336B", "HIS337B", "LEU338B",
                "GLN339B", "GLN340B", "GLN341B", "GLN342B")
cerrada_res <- c("ALA146A", "VAL147A", "GLU148A", "VAL149A", "SER150A", "GLU151A",
                "LEU152A", "THR153A", "GLU154A", "SER155A", "SER156A", "LYS157A",
                "LYS158A", "ASN159A", "VAL160A", "ALA161A", "HIS162A", "GLY163A",
                "LEU164A", "ALA165A", "TRP166A", "SER167A", "TYR168A", "TYR169A",
                "ILE170A", "GLY171A", "TYR172A", "LEU173A", "LYS174A", "VAL175A",
                "VAL176A", "LEU177A", "PRO178A", "ARG179A", "LEU180A", "LYS181A",
                "GLU182A", "CYS183A", "MET184A", "GLU185A", "GLU186A", "LEU187A",
                "SER188A", "ARG189A", "THR190A", "ASN191A", "PRO192A", "MET193A",
                "LEU194A", "ARG195A", "ALA196A", "HIS197A", "ARG198A", "ASP199A",
                "THR200A", "TRP201A", "LYS202A", "LEU203A", "HIS204A", "ILE205A",
                "LEU206A", "VAL207A", "PRO208A", "LEU209A", "GLY210A", "CYS211A",
                "ASP212A", "ILE213A", "TRP214A", "ASP215A", "ASP216A", "LEU217A",
                "GLU218A", "LYS219A", "ALA220A", "ASP221A", "SER222A", "ASN223A",
                "ILE224A", "GLN225A", "TYR226A", "LEU227A", "ALA228A", "ASP229A",
                "LEU230A", "PRO231A", "GLU232A", "THR233A", "ILE234A", "LEU235A",
                "THR236A", "ARG237A", "ALA238A", "GLY239A", "ILE240A", "LYS241A",
                "ARG242A", "ARG243A", "VAL244A", "TYR245A", "LYS246A", "HIS247A",
                "SER248A", "LEU249A", "TYR250A", "VAL251A", "ILE252A", "ARG253A",
                "ASP254A", "LYS255A", "ASP256A", "ASN257A", "LYS258A", "LEU259A",
                "ARG260A", "PRO261A", "CYS262A", "VAL263A", "LEU264A", "GLU265A",
                "PHE266A", "ALA267A", "SER268A", "PRO269A", "LEU270A", "GLN271A",
                "THR272A", "LEU273A", "CYS274A", "ALA275A", "MET276A", "SER277A", 
                "GLN278A", "ASP279A", "ASP280A", "CYS281A", "ALA282A", "ALA283A",
                "PHE284A", "SER285A", "ARG286A", "GLU287A", "GLN288A", "ARG289A",
                "LEU290A", "GLU291A", "GLN292A", "ALA293A", "ARG294A", "LEU295A",
                "PHE296A", "TYR297A", "ARG298A", "SER299A", "LEU300A", "ARG301A",
                "ASP302A", "ILE303A", "LEU304A", "GLY305A", "SER306A", "SER307A",
                "LYS308A", "GLU309A", "CYS310A", "ALA311A", "GLY312A", "LEU313A",
                "TYR314A", "ARG315A", "LEU316A", "ILE317A", "ALA318A", "TYR319A",
                "GLU320A", "GLU321A", "PRO322A", "ALA323A", "GLU324A", "PRO325A",
                "GLU326A", "SER327A", "HIS328A", "PHE329A", "LEU330A", "SER331A",
                "GLY332A", "LEU333A", "ILE334A", "LEU335A", "TRP336A", "HIS337A",
                "LEU338A", "GLN339A", "GLN340A", "GLN341A", "ALA146B",
                "VAL147B", "GLU148B", "VAL149B", "SER150B", "GLU151B", "LEU152B",
                "THR153B", "GLU154B", "SER155B", "SER156B", "LYS157B", "LYS158B",
                "ASN159B", "VAL160B", "ALA161B", "HIS162B", "GLY163B", "LEU164B",
                "ALA165B", "TRP166B", "SER167B", "TYR168B", "TYR169B", "ILE170B",
                "GLY171B", "TYR172B", "LEU173B", "LYS174B", "VAL175B", "VAL176B",
                "LEU177B", "PRO178B", "ARG179B", "LEU180B", "LYS181B", "GLU182B",
                "CYS183B", "MET184B", "GLU185B", "GLU186B", "LEU187B", "SER188B",
                "ARG189B", "THR190B", "ASN191B", "PRO192B", "MET193B", "LEU194B",
                "ARG195B", "ALA196B", "HIS197B", "ARG198B", "ASP199B", "THR200B",
                "TRP201B", "LYS202B", "LEU203B", "HIS204B", "ILE205B", "LEU206B",
                "VAL207B", "PRO208B", "LEU209B", "GLY210B", "CYS211B", "ASP212B",
                "ILE213B", "TRP214B", "ASP215B", "ASP216B", "LEU217B", "GLU218B",
                "LYS219B", "ALA220B", "ASP221B", "SER222B", "ASN223B", "ILE224B",
                "GLN225B", "TYR226B", "LEU227B", "ALA228B", "ASP229B", "LEU230B",
                "PRO231B", "GLU232B", "THR233B", "ILE234B", "LEU235B", "THR236B",
                "ARG237B", "ALA238B", "GLY239B", "ILE240B", "LYS241B", "ARG242B",
                "ARG243B", "VAL244B", "TYR245B", "LYS246B", "HIS247B", "SER248B",
                "LEU249B", "TYR250B", "VAL251B", "ILE252B", "ARG253B", "ASP254B",
                "LYS255B", "ASP256B", "ASN257B", "LYS258B", "LEU259B", "ARG260B",
                "PRO261B", "CYS262B", "VAL263B", "LEU264B", "GLU265B", "PHE266B",
                "ALA267B", "SER268B", "PRO269B", "LEU270B", "GLN271B", "THR272B",
                "LEU273B", "CYS274B", "ALA275B", "MET276B", "SER277B", "GLN278B",
                "ASP279B", "ASP280B", "CYS281B", "ALA282B", "ALA283B", "PHE284B",
                "SER285B", "ARG286B", "GLU287B", "GLN288B", "ARG289B", "LEU290B",
                "GLU291B", "GLN292B", "ALA293B", "ARG294B", "LEU295B", "PHE296B",
                "TYR297B", "ARG298B", "SER299B", "LEU300B", "ARG301B", "ASP302B",
                "ILE303B", "LEU304B", "GLY305B", "SER306B", "SER307B", "LYS308B",
                "GLU309B", "CYS310B", "ALA311B", "GLY312B", "LEU313B", "TYR314B",
                "ARG315B", "LEU316B", "ILE317B", "ALA318B", "TYR319B", "GLU320B",
                "GLU321B", "PRO322B", "ALA323B", "GLU324B", "PRO325B", "GLU326B",
                "SER327B", "HIS328B", "PHE329B", "LEU330B", "SER331B", "GLY332B",
                "LEU333B", "ILE334B", "LEU335B", "TRP336B", "HIS337B", "LEU338B",
                "GLN339B", "GLN340B", "GLN341B")
names(abierta_bet)[] <- abierta_res
names(abierta_net)[] <- abierta_res
names(cerrada_bet)[] <- cerrada_res
names(cerrada_net)[] <- cerrada_res

c_net <- network(cerrada_net)
a_net <- network(abierta_net)
plot(c_net, vertex.cex = 1)
plot(a_net, vertex.cex = 1)
ggraph(c_net) + geom_edge_link() + geom_residues_point() + theme_graph()
ggraph(a_net) + geom_edge_link() + geom_residues_point() + theme_graph()


abierta_values <- data.frame("residues" = abierta_res,
                             "degree" = c(1:394), 
                             "betweeness" = c(1:394),
                             "abs_degree" = c(1:394),
                             "P" = 0, 
                             "z" = 0)
for(i in 1:394){
  abierta_values$degree[i] = sum(abierta_net[,i])
  abierta_values$betweeness[i] = sum(abierta_bet[,i])
}
for(i in 1:394){
  abierta_values$abs_degree[i] = 0
  for(j in 1:394){
    if(abierta_net[j,i] > 0){
      abierta_values$abs_degree[i] = abierta_values$abs_degree[i] + 1
    }
  }
}

cerrada_values <- data.frame("residues" = cerrada_res,
                             "degree" = c(1:392), 
                             "betweeness" = c(1:392),
                             "abs_degree" = c(1:392),
                             "P" = 0,
                             "z" = 0)
for(i in 1:392){
  cerrada_values$degree[i] = sum(cerrada_net[,i])
  cerrada_values$betweeness[i] = sum(cerrada_bet[,i])
}
for(i in 1:392){
  cerrada_values$abs_degree[i] = 0
  for(j in 1:392){
    if(cerrada_net[j,i] > 0){
      cerrada_values$abs_degree[i] = cerrada_values$abs_degree[i] + 1
    }
  }
}
residues <- merge(x = abierta_values, y = cerrada_values, by = c("residues"), all = TRUE)
names(residues)[] <- c("residues", "ab_degree", "ab_betweeness", "ab_abs", "ab_P", "ab_z", "ce_degree", "ce_betweeness", "ce_abs", "ce_P", "ce_z")
ggpairs(residues[,c(2:7)])
residues <- residues[-which(is.na(residues$ce_abs)),] 
residues$ab_cluster <- 0
residues$ce_cluster <- 0

#### COMUNIDADES ESTRUCTURA ABIERTA ####
abierta_index <- data.frame("num" = (0:393),
                            "res" = as.character(abierta_res))
abierta_index$res <- as.character(abierta_index$res)
abierta_comm1 <- data.frame("residues" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 210, 212, 213, 216, 258, 346, 347, 348, 349, 350, 351, 352, 353, 367, 368, 369),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
  
for(i in 1:28){
  abierta_comm1$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm1$residues[i])]
  abierta_comm1$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm1$residues[i])]
  abierta_comm1$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm1$residues[i])]
}
for(i in 1:28){
  abierta_comm1$connectivity[i] = (abierta_comm1$degree[i] - mean(abierta_comm1$degree))/std(abierta_comm1$degree)
}
aux <- abierta_net[match(abierta_comm1$residues, names(abierta_net)),match(abierta_comm1$residues, names(abierta_net))]
for(i in 1:28){
  abierta_comm1$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm1$partition = 1 - ((abierta_comm1$ksi/abierta_comm1$degree)^2)
for(i in 1:28){
  residues$ab_P[which(residues$residues == abierta_comm1$residues[i])] = abierta_comm1$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm1$residues[i])] = abierta_comm1$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm1$residues[i])] = 0
}
nrow(subset(residues, ab_cluster == 0))

abierta_comm2 <- data.frame("residues" = c(292),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
  
abierta_comm2$residues[1] = abierta_index$res[which(abierta_index$num == abierta_comm2$residues[1])]
abierta_comm2$degree[1] = abierta_values$degree[which(abierta_values$residues == abierta_comm2$residues[1])]
abierta_comm2$abs_degree[1] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm2$residues[1])]
abierta_comm2$connectivity[1] = (abierta_comm2$degree[1] - mean(abierta_comm2$degree))/std(abierta_comm2$degree)
aux <- abierta_net[match(abierta_comm2$residues, names(abierta_net)),match(abierta_comm2$residues, names(abierta_net))]
abierta_comm2$ksi[] = sum(aux)
rm(aux)
abierta_comm2$partition = 1 - ((abierta_comm2$ksi/abierta_comm2$degree)^2)
residues$ab_P[which(residues$residues == abierta_comm2$residues[1])] = abierta_comm2$partition[1]
residues$ab_z[which(residues$residues == abierta_comm2$residues[1])] = abierta_comm2$connectivity[1]
residues$ab_cluster[which(residues$residues == abierta_comm2$residues[1])] = 1
nrow(subset(residues, ab_cluster == 1))


abierta_comm3 <- data.frame("residues" = c(12, 13, 14, 15, 18, 19, 22, 23, 25, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 134, 136, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 159, 160, 161, 162, 170, 172, 174),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:65){
  abierta_comm3$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm3$residues[i])]
  abierta_comm3$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm3$residues[i])]
  abierta_comm3$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm3$residues[i])]
}
for(i in 1:65){
  abierta_comm3$connectivity[i] = (abierta_comm3$degree[i] - mean(abierta_comm3$degree))/std(abierta_comm3$degree)
}
aux <- abierta_net[match(abierta_comm3$residues, names(abierta_net)),match(abierta_comm3$residues, names(abierta_net))]
for(i in 1:65){
  abierta_comm3$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm3$partition = 1 - ((abierta_comm3$ksi/abierta_comm3$degree)^2)
for(i in 1:65){
  residues$ab_P[which(residues$residues == abierta_comm3$residues[i])] = abierta_comm3$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm3$residues[i])] = abierta_comm3$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm3$residues[i])] = 2
}
nrow(subset(residues, ab_cluster == 2))

abierta_comm4 <- data.frame("residues" = c(197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:12){
  abierta_comm4$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm4$residues[i])]
  abierta_comm4$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm4$residues[i])]
  abierta_comm4$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm4$residues[i])]
}
for(i in 1:12){
  abierta_comm4$connectivity[i] = (abierta_comm4$degree[i] - mean(abierta_comm4$degree))/std(abierta_comm4$degree)
}
aux <- abierta_net[match(abierta_comm4$residues, names(abierta_net)),match(abierta_comm4$residues, names(abierta_net))]
for(i in 1:12){
  abierta_comm4$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm4$partition = 1 - ((abierta_comm4$ksi/abierta_comm4$degree)^2)
for(i in 1:12){
  residues$ab_P[which(residues$residues == abierta_comm4$residues[i])] = abierta_comm4$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm4$residues[i])] = abierta_comm4$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm4$residues[i])] = 3
}
nrow(subset(residues, ab_cluster == 3))

abierta_comm5 <- data.frame("residues" = c(266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 278, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 315, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:49){
  abierta_comm5$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm5$residues[i])]
  abierta_comm5$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm5$residues[i])]
  abierta_comm5$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm5$residues[i])]
}
for(i in 1:49){
  abierta_comm5$connectivity[i] = (abierta_comm5$degree[i] - mean(abierta_comm5$degree))/std(abierta_comm5$degree)
}
aux <- abierta_net[match(abierta_comm5$residues, names(abierta_net)),match(abierta_comm5$residues, names(abierta_net))]
for(i in 1:49){
  abierta_comm5$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm5$partition = 1 - ((abierta_comm5$ksi/abierta_comm5$degree)^2)
for(i in 1:49){
  residues$ab_P[which(residues$residues == abierta_comm5$residues[i])] = abierta_comm5$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm5$residues[i])] = abierta_comm5$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm5$residues[i])] = 4
}
nrow(subset(residues, ab_cluster == 4))

abierta_comm6 <- data.frame("residues" = c(95, 96, 97, 230, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 277, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288, 289, 290, 291, 293, 295, 297, 299, 300),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:44){
  abierta_comm6$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm6$residues[i])]
  abierta_comm6$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm6$residues[i])]
  abierta_comm6$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm6$residues[i])]
}
for(i in 1:44){
  abierta_comm6$connectivity[i] = (abierta_comm6$degree[i] - mean(abierta_comm6$degree))/std(abierta_comm6$degree)
}
aux <- abierta_net[match(abierta_comm6$residues, names(abierta_net)),match(abierta_comm6$residues, names(abierta_net))]
for(i in 1:44){
  abierta_comm6$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm6$partition = 1 - ((abierta_comm6$ksi/abierta_comm6$degree)^2)
for(i in 1:44){
  residues$ab_P[which(residues$residues == abierta_comm6$residues[i])] = abierta_comm6$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm6$residues[i])] = abierta_comm6$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm6$residues[i])] = 5
}
nrow(subset(residues, ab_cluster == 5))

abierta_comm7 <- data.frame("residues" = c(26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46, 50, 52, 54, 56, 57, 58, 59, 60, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 98, 99, 100, 101, 117, 119, 158, 163, 164, 165, 166, 167, 168, 169, 171),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:55){
  abierta_comm7$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm7$residues[i])]
  abierta_comm7$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm7$residues[i])]
  abierta_comm7$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm7$residues[i])]
}
for(i in 1:55){
  abierta_comm7$connectivity[i] = (abierta_comm7$degree[i] - mean(abierta_comm7$degree))/std(abierta_comm7$degree)
}
aux <- abierta_net[match(abierta_comm7$residues, names(abierta_net)),match(abierta_comm7$residues, names(abierta_net))]
for(i in 1:55){
  abierta_comm7$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm7$partition = 1 - ((abierta_comm7$ksi/abierta_comm7$degree)^2)
for(i in 1:55){
  residues$ab_P[which(residues$residues == abierta_comm7$residues[i])] = abierta_comm7$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm7$residues[i])] = abierta_comm7$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm7$residues[i])] = 6
}
nrow(subset(residues, ab_cluster == 6))

abierta_comm8 <- data.frame("residues" = c(45, 47, 48, 49, 51, 53, 55, 75, 76, 77, 78, 79, 80, 81, 82, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:45){
  abierta_comm8$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm8$residues[i])]
  abierta_comm8$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm8$residues[i])]
  abierta_comm8$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm8$residues[i])]
}
for(i in 1:45){
  abierta_comm8$connectivity[i] = (abierta_comm8$degree[i] - mean(abierta_comm8$degree))/std(abierta_comm8$degree)
}
aux <- abierta_net[match(abierta_comm8$residues, names(abierta_net)),match(abierta_comm8$residues, names(abierta_net))]
for(i in 1:45){
  abierta_comm8$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm8$partition = 1 - ((abierta_comm8$ksi/abierta_comm8$degree)^2)
for(i in 1:45){
  residues$ab_P[which(residues$residues == abierta_comm8$residues[i])] = abierta_comm8$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm8$residues[i])] = abierta_comm8$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm8$residues[i])] = 7
}
nrow(subset(residues, ab_cluster == 7))

abierta_comm9 <- data.frame("residues" = c(133, 135, 137, 214, 215, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 231, 253, 254, 255, 256, 294, 296, 298, 314, 316, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363, 364, 365, 366),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:41){
  abierta_comm9$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm9$residues[i])]
  abierta_comm9$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm9$residues[i])]
  abierta_comm9$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm9$residues[i])]
}
for(i in 1:41){
  abierta_comm9$connectivity[i] = (abierta_comm9$degree[i] - mean(abierta_comm9$degree))/std(abierta_comm9$degree)
}
aux <- abierta_net[match(abierta_comm9$residues, names(abierta_net)),match(abierta_comm9$residues, names(abierta_net))]
for(i in 1:41){
  abierta_comm9$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm9$partition = 1 - ((abierta_comm9$ksi/abierta_comm9$degree)^2)
for(i in 1:41){
  residues$ab_P[which(residues$residues == abierta_comm9$residues[i])] = abierta_comm9$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm9$residues[i])] = abierta_comm9$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm9$residues[i])] = 8
}
nrow(subset(residues, ab_cluster == 8))

abierta_comm10 <- data.frame("residues" = c(16, 17, 20, 21, 24, 209, 211, 257, 259, 260, 261, 262, 263, 264, 265, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 370),
                             "degree" = 0,
                             "partition" = 0,
                             "connectivity" = 0,
                             "abs_degree" = 0,
                             "ksi" = 0)
for(i in 1:45){
  abierta_comm10$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm10$residues[i])]
  abierta_comm10$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm10$residues[i])]
  abierta_comm10$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm10$residues[i])]
}
for(i in 1:45){
  abierta_comm10$connectivity[i] = (abierta_comm10$degree[i] - mean(abierta_comm10$degree))/std(abierta_comm10$degree)
}
aux <- abierta_net[match(abierta_comm10$residues, names(abierta_net)),match(abierta_comm10$residues, names(abierta_net))]
for(i in 1:45){
  abierta_comm10$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm10$partition = 1 - ((abierta_comm10$ksi/abierta_comm10$degree)^2)
for(i in 1:45){
  residues$ab_P[which(residues$residues == abierta_comm10$residues[i])] = abierta_comm10$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm10$residues[i])] = abierta_comm10$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm10$residues[i])] = 9
}
nrow(subset(residues, ab_cluster == 9))

abierta_comm11 <- data.frame("residues" = c(173, 175, 176, 177, 178, 179, 180, 181, 182),
                             "degree" = 0,
                             "partition" = 0,
                             "connectivity" = 0,
                             "abs_degree" = 0,
                             "ksi" = 0)
for(i in 1:9){
  abierta_comm11$residues[i] = abierta_index$res[which(abierta_index$num == abierta_comm11$residues[i])]
  abierta_comm11$degree[i] = abierta_values$degree[which(abierta_values$residues == abierta_comm11$residues[i])]
  abierta_comm11$abs_degree[i] = abierta_values$abs_degree[which(abierta_values$residues == abierta_comm11$residues[i])]
}
for(i in 1:9){
  abierta_comm11$connectivity[i] = (abierta_comm11$degree[i] - mean(abierta_comm11$degree))/std(abierta_comm11$degree)
}
aux <- abierta_net[match(abierta_comm11$residues, names(abierta_net)),match(abierta_comm11$residues, names(abierta_net))]
for(i in 1:9){
  abierta_comm11$ksi[i] = sum(aux[,i])
}
rm(aux)
abierta_comm11$partition = 1 - ((abierta_comm11$ksi/abierta_comm11$degree)^2)
for(i in 1:9){
  residues$ab_P[which(residues$residues == abierta_comm11$residues[i])] = abierta_comm11$partition[i]
  residues$ab_z[which(residues$residues == abierta_comm11$residues[i])] = abierta_comm11$connectivity[i]
  residues$ab_cluster[which(residues$residues == abierta_comm11$residues[i])] = 10
}
nrow(subset(residues, ab_cluster == 10))

#### COMUNIDADES STING CERRADA ####
cerrada_index <- data.frame("num" = (0:391),
                            "res" = as.character(cerrada_res))
cerrada_index$res <- as.character(cerrada_index$res)
cerrada_comm1 <- data.frame("residues" = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 63, 65, 122, 124, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 149, 174, 206, 207, 208),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)

for(i in 1:44){
  cerrada_comm1$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm1$residues[i])]
  cerrada_comm1$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm1$residues[i])]
  cerrada_comm1$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm1$residues[i])]
}
for(i in 1:44){
  cerrada_comm1$connectivity[i] = (cerrada_comm1$degree[i] - mean(cerrada_comm1$degree))/std(cerrada_comm1$degree)
}
aux <- cerrada_net[match(cerrada_comm1$residues, names(cerrada_net)),match(cerrada_comm1$residues, names(cerrada_net))]
for(i in 1:44){
  cerrada_comm1$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm1$partition = 1 - ((cerrada_comm1$ksi/cerrada_comm1$degree)^2)
for(i in 1:44){
  residues$ce_P[which(residues$residues == cerrada_comm1$residues[i])] = cerrada_comm1$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm1$residues[i])] = cerrada_comm1$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm1$residues[i])] = 0
}
nrow(subset(residues, ce_cluster == 0))

cerrada_comm2 <- data.frame("residues" = c(179),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)

cerrada_comm2$residues[1] = cerrada_index$res[which(cerrada_index$num == cerrada_comm2$residues[1])]
cerrada_comm2$degree[1] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm2$residues[1])]
cerrada_comm2$abs_degree[1] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm2$residues[1])]
cerrada_comm2$connectivity[1] = (cerrada_comm2$degree[1] - mean(cerrada_comm2$degree))/std(cerrada_comm2$degree)
aux <- cerrada_net[match(cerrada_comm2$residues, names(cerrada_net)),match(cerrada_comm2$residues, names(cerrada_net))]
cerrada_comm2$ksi[] = sum(aux)
rm(aux)
cerrada_comm2$partition = 1 - ((cerrada_comm2$ksi/cerrada_comm2$degree)^2)
residues$ce_P[which(residues$residues == cerrada_comm2$residues[1])] = cerrada_comm2$partition[1]
residues$ce_z[which(residues$residues == cerrada_comm2$residues[1])] = cerrada_comm2$connectivity[1]
residues$ce_cluster[which(residues$residues == cerrada_comm2$residues[1])] = 1
nrow(subset(residues, ce_cluster == 1))


cerrada_comm3 <- data.frame("residues" = c(305, 307),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:2){
  cerrada_comm3$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm3$residues[i])]
  cerrada_comm3$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm3$residues[i])]
  cerrada_comm3$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm3$residues[i])]
}
for(i in 1:2){
  cerrada_comm3$connectivity[i] = (cerrada_comm3$degree[i] - mean(cerrada_comm3$degree))/std(cerrada_comm3$degree)
}
aux <- cerrada_net[match(cerrada_comm3$residues, names(cerrada_net)),match(cerrada_comm3$residues, names(cerrada_net))]
for(i in 1:2){
  cerrada_comm3$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm3$partition = 1 - ((cerrada_comm3$ksi/cerrada_comm3$degree)^2)
for(i in 1:2){
  residues$ce_P[which(residues$residues == cerrada_comm3$residues[i])] = cerrada_comm3$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm3$residues[i])] = cerrada_comm3$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm3$residues[i])] = 2
}
nrow(subset(residues, ce_cluster == 2))


cerrada_comm4 <- data.frame("residues" = c(92, 93, 94, 95, 254, 256, 258, 260, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 277, 300, 301, 302, 303, 304, 306, 308, 309, 310, 312, 314, 316, 365, 367, 369, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:59){
  cerrada_comm4$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm4$residues[i])]
  cerrada_comm4$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm4$residues[i])]
  cerrada_comm4$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm4$residues[i])]
}
for(i in 1:59){
  cerrada_comm4$connectivity[i] = (cerrada_comm4$degree[i] - mean(cerrada_comm4$degree))/std(cerrada_comm4$degree)
}
aux <- cerrada_net[match(cerrada_comm4$residues, names(cerrada_net)),match(cerrada_comm4$residues, names(cerrada_net))]
for(i in 1:59){
  cerrada_comm4$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm4$partition = 1 - ((cerrada_comm4$ksi/cerrada_comm4$degree)^2)
for(i in 1:59){
  residues$ce_P[which(residues$residues == cerrada_comm4$residues[i])] = cerrada_comm4$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm4$residues[i])] = cerrada_comm4$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm4$residues[i])] = 3
}
nrow(subset(residues, ce_cluster == 3))

cerrada_comm5 <- data.frame("residues" = c(30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 48, 49, 50, 52, 55, 57, 83, 84, 85, 86, 87, 88, 89, 90, 98, 99, 100, 101, 117, 119, 285, 287, 290, 292),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:41){
  cerrada_comm5$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm5$residues[i])]
  cerrada_comm5$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm5$residues[i])]
  cerrada_comm5$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm5$residues[i])]
}
for(i in 1:41){
  cerrada_comm5$connectivity[i] = (cerrada_comm5$degree[i] - mean(cerrada_comm5$degree))/std(cerrada_comm5$degree)
}
aux <- cerrada_net[match(cerrada_comm5$residues, names(cerrada_net)),match(cerrada_comm5$residues, names(cerrada_net))]
for(i in 1:41){
  cerrada_comm5$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm5$partition = 1 - ((cerrada_comm5$ksi/cerrada_comm5$degree)^2)
for(i in 1:41){
  residues$ce_P[which(residues$residues == cerrada_comm5$residues[i])] = cerrada_comm5$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm5$residues[i])] = cerrada_comm5$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm5$residues[i])] = 4
}
nrow(subset(residues, ce_cluster == 4))

cerrada_comm6 <- data.frame("residues" = c(322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:12){
  cerrada_comm6$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm6$residues[i])]
  cerrada_comm6$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm6$residues[i])]
  cerrada_comm6$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm6$residues[i])]
}
for(i in 1:12){
  cerrada_comm6$connectivity[i] = (cerrada_comm6$degree[i] - mean(cerrada_comm6$degree))/std(cerrada_comm6$degree)
}
aux <- cerrada_net[match(cerrada_comm6$residues, names(cerrada_net)),match(cerrada_comm6$residues, names(cerrada_net))]
for(i in 1:12){
  cerrada_comm6$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm6$partition = 1 - ((cerrada_comm6$ksi/cerrada_comm6$degree)^2)
for(i in 1:12){
  residues$ce_P[which(residues$residues == cerrada_comm6$residues[i])] = cerrada_comm6$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm6$residues[i])] = cerrada_comm6$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm6$residues[i])] = 5
}
nrow(subset(residues, ce_cluster == 5))

cerrada_comm7 <- data.frame("residues" = c(209, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 223, 257, 259, 261, 317, 318, 319, 320, 321, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 361, 366, 368, 370),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:50){
  cerrada_comm7$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm7$residues[i])]
  cerrada_comm7$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm7$residues[i])]
  cerrada_comm7$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm7$residues[i])]
}
for(i in 1:50){
  cerrada_comm7$connectivity[i] = (cerrada_comm7$degree[i] - mean(cerrada_comm7$degree))/std(cerrada_comm7$degree)
}
aux <- cerrada_net[match(cerrada_comm7$residues, names(cerrada_net)),match(cerrada_comm7$residues, names(cerrada_net))]
for(i in 1:50){
  cerrada_comm7$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm7$partition = 1 - ((cerrada_comm7$ksi/cerrada_comm7$degree)^2)
for(i in 1:50){
  residues$ce_P[which(residues$residues == cerrada_comm7$residues[i])] = cerrada_comm7$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm7$residues[i])] = cerrada_comm7$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm7$residues[i])] = 6
}
nrow(subset(residues, ce_cluster == 6))

cerrada_comm8 <- data.frame("residues" = c(47, 58, 60, 62, 64, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 118, 120, 125, 169, 171, 173, 175, 176, 177, 178, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 286, 288, 289, 291),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:67){
  cerrada_comm8$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm8$residues[i])]
  cerrada_comm8$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm8$residues[i])]
  cerrada_comm8$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm8$residues[i])]
}
for(i in 1:67){
  cerrada_comm8$connectivity[i] = (cerrada_comm8$degree[i] - mean(cerrada_comm8$degree))/std(cerrada_comm8$degree)
}
aux <- cerrada_net[match(cerrada_comm8$residues, names(cerrada_net)),match(cerrada_comm8$residues, names(cerrada_net))]
for(i in 1:67){
  cerrada_comm8$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm8$partition = 1 - ((cerrada_comm8$ksi/cerrada_comm8$degree)^2)
for(i in 1:67){
  residues$ce_P[which(residues$residues == cerrada_comm8$residues[i])] = cerrada_comm8$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm8$residues[i])] = cerrada_comm8$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm8$residues[i])] = 7
}
nrow(subset(residues, ce_cluster == 7))

cerrada_comm9 <- data.frame("residues" = c(96, 97, 222, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 255, 276, 278, 279, 280, 281, 282, 283, 284, 293, 294, 295, 296, 297, 298, 299, 311, 313, 315, 359, 360, 362, 363, 364),
                            "degree" = 0,
                            "partition" = 0,
                            "connectivity" = 0,
                            "abs_degree" = 0,
                            "ksi" = 0)
for(i in 1:57){
  cerrada_comm9$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm9$residues[i])]
  cerrada_comm9$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm9$residues[i])]
  cerrada_comm9$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm9$residues[i])]
}
for(i in 1:57){
  cerrada_comm9$connectivity[i] = (cerrada_comm9$degree[i] - mean(cerrada_comm9$degree))/std(cerrada_comm9$degree)
}
aux <- cerrada_net[match(cerrada_comm9$residues, names(cerrada_net)),match(cerrada_comm9$residues, names(cerrada_net))]
for(i in 1:57){
  cerrada_comm9$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm9$partition = 1 - ((cerrada_comm9$ksi/cerrada_comm9$degree)^2)
for(i in 1:57){
  residues$ce_P[which(residues$residues == cerrada_comm9$residues[i])] = cerrada_comm9$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm9$residues[i])] = cerrada_comm9$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm9$residues[i])] = 8
}
nrow(subset(residues, ce_cluster == 8))

cerrada_comm10 <- data.frame("residues" = c(196, 197, 198, 199, 200, 201, 202, 203, 204, 205),
                             "degree" = 0,
                             "partition" = 0,
                             "connectivity" = 0,
                             "abs_degree" = 0,
                             "ksi" = 0)
for(i in 1:10){
  cerrada_comm10$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm10$residues[i])]
  cerrada_comm10$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm10$residues[i])]
  cerrada_comm10$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm10$residues[i])]
}
for(i in 1:10){
  cerrada_comm10$connectivity[i] = (cerrada_comm10$degree[i] - mean(cerrada_comm10$degree))/std(cerrada_comm10$degree)
}
aux <- cerrada_net[match(cerrada_comm10$residues, names(cerrada_net)),match(cerrada_comm10$residues, names(cerrada_net))]
for(i in 1:10){
  cerrada_comm10$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm10$partition = 1 - ((cerrada_comm10$ksi/cerrada_comm10$degree)^2)
for(i in 1:10){
  residues$ce_P[which(residues$residues == cerrada_comm10$residues[i])] = cerrada_comm10$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm10$residues[i])] = cerrada_comm10$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm10$residues[i])] = 9
}
nrow(subset(residues, ce_cluster == 9))

cerrada_comm11 <- data.frame("residues" = c(12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 51, 53, 54, 56, 59, 61, 91, 121, 123, 148, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 170, 172, 210),
                             "degree" = 0,
                             "partition" = 0,
                             "connectivity" = 0,
                             "abs_degree" = 0,
                             "ksi" = 0)
for(i in 1:49){
  cerrada_comm11$residues[i] = cerrada_index$res[which(cerrada_index$num == cerrada_comm11$residues[i])]
  cerrada_comm11$degree[i] = cerrada_values$degree[which(cerrada_values$residues == cerrada_comm11$residues[i])]
  cerrada_comm11$abs_degree[i] = cerrada_values$abs_degree[which(cerrada_values$residues == cerrada_comm11$residues[i])]
}
for(i in 1:49){
  cerrada_comm11$connectivity[i] = (cerrada_comm11$degree[i] - mean(cerrada_comm11$degree))/std(cerrada_comm11$degree)
}
aux <- cerrada_net[match(cerrada_comm11$residues, names(cerrada_net)),match(cerrada_comm11$residues, names(cerrada_net))]
for(i in 1:49){
  cerrada_comm11$ksi[i] = sum(aux[,i])
}
rm(aux)
cerrada_comm11$partition = 1 - ((cerrada_comm11$ksi/cerrada_comm11$degree)^2)
for(i in 1:49){
  residues$ce_P[which(residues$residues == cerrada_comm11$residues[i])] = cerrada_comm11$partition[i]
  residues$ce_z[which(residues$residues == cerrada_comm11$residues[i])] = cerrada_comm11$connectivity[i]
  residues$ce_cluster[which(residues$residues == cerrada_comm11$residues[i])] = 10
}
nrow(subset(residues, ce_cluster == 10))

#### GRAFICOS ####
ggplot(residues, aes(x = ab_P, y = ab_z, color = as.factor(ab_cluster))) + 
  geom_point() + a + xlab("P (partition coefficient)") + 
  ylab("z (intramolecule connectivity)")
ggplot(residues, aes(x = ce_P, y = ce_z, color = as.factor(ce_cluster))) + 
  geom_point() + a + xlab("P (partition coefficient)") + 
  ylab("z (intramolecule connectivity)")
residues$atom = c(146:341, 343:538)
ggplot(residues, aes(x = atom)) + 
  geom_line(aes(y = ab_P), color = "darkgreen") +
  geom_point(aes(y = ab_P), color = "darkgreen") +
  geom_line(aes(y = ce_P), color = "magenta4") + 
  geom_point(aes(y = ce_P), color = "magenta4") +
  xlab("residue") + ylab("P")
residues$AP = residues$ab_P - residues$ce_P
ggplot(residues, aes(x = atom)) + 
  geom_line(aes(y = AP), color = "dodgerblue3") +
  geom_point(aes(y = AP), color = "dodgerblue3") +
  xlab("residue") + ylab("P")
residues[which(residues$AP > 0.5),]
residues[which(residues$AP > 0.5),]

ggplot(residues, aes(x = atom)) + 
  geom_point(aes(y = ab_betweeness), color = "darkgreen") + 
  geom_point(aes(y = ce_betweeness), color = "magenta4")

ggplot(residues, aes(x = atom, y = ab_betweeness-ce_betweeness)) + geom_point(color = "dodgerblue4")
      