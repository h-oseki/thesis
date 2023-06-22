#### ABIERTA REP 1 ####

abierta_rep1_GLN278A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-GLN278A.xvg")
abierta_rep1_GLN278A <- as.data.frame(sapply(abierta_rep1_GLN278A, as.numeric))
abierta_rep1_GLN278B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-GLN278B.xvg")
abierta_rep1_GLN278B <- as.data.frame(sapply(abierta_rep1_GLN278B, as.numeric))
abierta_rep1_ASP279A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ASP279A.xvg")
abierta_rep1_ASP279A <- as.data.frame(sapply(abierta_rep1_ASP279A, as.numeric))
abierta_rep1_ASP279B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ASP279B.xvg")
abierta_rep1_ASP279B <- as.data.frame(sapply(abierta_rep1_ASP279B, as.numeric))
abierta_rep1_ASP280A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ASP280A.xvg")
abierta_rep1_ASP280A <- as.data.frame(sapply(abierta_rep1_ASP280A, as.numeric))
abierta_rep1_ASP280B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ASP280B.xvg")
abierta_rep1_ASP280B <- as.data.frame(sapply(abierta_rep1_ASP280B, as.numeric))
abierta_rep1_CYS281A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-CYS281A.xvg")
abierta_rep1_CYS281A <- as.data.frame(sapply(abierta_rep1_CYS281A, as.numeric))
abierta_rep1_CYS281B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-CYS281B.xvg")
abierta_rep1_CYS281B <- as.data.frame(sapply(abierta_rep1_CYS281B, as.numeric))
abierta_rep1_ALA282A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ALA282A.xvg")
abierta_rep1_ALA282A <- as.data.frame(sapply(abierta_rep1_ALA282A, as.numeric))
abierta_rep1_ALA282B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ALA282B.xvg")
abierta_rep1_ALA282B <- as.data.frame(sapply(abierta_rep1_ALA282B, as.numeric))
abierta_rep1_ALA283A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ALA283A.xvg")
abierta_rep1_ALA283A <- as.data.frame(sapply(abierta_rep1_ALA283A, as.numeric))
abierta_rep1_ALA283B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ALA283B.xvg")
abierta_rep1_ALA283B <- as.data.frame(sapply(abierta_rep1_ALA283B, as.numeric))
abierta_rep1_PHE284A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-PHE284A.xvg")
abierta_rep1_PHE284A <- as.data.frame(sapply(abierta_rep1_PHE284A, as.numeric))
abierta_rep1_PHE284B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-PHE284B.xvg")
abierta_rep1_PHE284B <- as.data.frame(sapply(abierta_rep1_PHE284B, as.numeric))
abierta_rep1_SER285A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-SER285A.xvg")
abierta_rep1_SER285A <- as.data.frame(sapply(abierta_rep1_SER285A, as.numeric))
abierta_rep1_SER285B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-SER285B.xvg")
abierta_rep1_SER285B <- as.data.frame(sapply(abierta_rep1_SER285B, as.numeric))
abierta_rep1_ARG286A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ARG286A.xvg")
abierta_rep1_ARG286A <- as.data.frame(sapply(abierta_rep1_ARG286A, as.numeric))
abierta_rep1_ARG286B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-ARG286B.xvg")
abierta_rep1_ARG286B <- as.data.frame(sapply(abierta_rep1_ARG286B, as.numeric))
abierta_rep1_GLU287A <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-GLU287A.xvg")
abierta_rep1_GLU287A <- as.data.frame(sapply(abierta_rep1_GLU287A, as.numeric))
abierta_rep1_GLU287B <- readXVG("abierta-rep1/crop-200/distance/abierta-rep1-GLU287B.xvg")
abierta_rep1_GLU287B <- as.data.frame(sapply(abierta_rep1_GLU287B, as.numeric))
abierta_rep1_protein <- read.delim("abierta-rep1/crop-200/distance/abierta1_protein.xvg")


abierta_rep1 <- data.frame('time' = abierta_rep1_protein$time/1000,
                           "gln278a" = sqrt((abierta_rep1_GLN278A$`r_278_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_GLN278A$`r_278_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_GLN278A$`r_278_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "gln278b" = sqrt((abierta_rep1_GLN278B$`r_475_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_GLN278B$`r_475_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_GLN278B$`r_475_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "asp279a" = sqrt((abierta_rep1_ASP279A$`r_279_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ASP279A$`r_279_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ASP279A$`r_279_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "asp279b" = sqrt((abierta_rep1_ASP279B$`r_476_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ASP279B$`r_476_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ASP279B$`r_476_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "asp280a" = sqrt((abierta_rep1_ASP280A$`r_280_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ASP280A$`r_280_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ASP280A$`r_280_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "asp280b" = sqrt((abierta_rep1_ASP280B$`r_477_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ASP280B$`r_477_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ASP280B$`r_477_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "cys281a" = sqrt((abierta_rep1_CYS281A$`r_281_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_CYS281A$`r_281_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_CYS281A$`r_281_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "cys281b" = sqrt((abierta_rep1_CYS281B$`r_478_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_CYS281B$`r_478_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_CYS281B$`r_478_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ala282a" = sqrt((abierta_rep1_ALA282A$`r_282_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ALA282A$`r_282_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ALA282A$`r_282_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ala282b" = sqrt((abierta_rep1_ALA282B$`r_479_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ALA282B$`r_479_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ALA282B$`r_479_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ala283a" = sqrt((abierta_rep1_ALA283A$`r_283_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ALA283A$`r_283_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ALA283A$`r_283_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ala283b" = sqrt((abierta_rep1_ALA283B$`r_480_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ALA283B$`r_480_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ALA283B$`r_480_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "phe284a" = sqrt((abierta_rep1_PHE284A$`r_284_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_PHE284A$`r_284_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_PHE284A$`r_284_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "phe284b" = sqrt((abierta_rep1_PHE284B$`r_481_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_PHE284B$`r_481_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_PHE284B$`r_481_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ser285a" = sqrt((abierta_rep1_SER285A$`r_285_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_SER285A$`r_285_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_SER285A$`r_285_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "ser285b" = sqrt((abierta_rep1_SER285B$`r_482_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_SER285B$`r_482_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_SER285B$`r_482_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "arg286a" = sqrt((abierta_rep1_ARG286A$`r_286_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ARG286A$`r_286_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ARG286A$`r_286_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "arg286b" = sqrt((abierta_rep1_ARG286B$`r_483_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_ARG286B$`r_483_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_ARG286B$`r_483_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "glu287a" = sqrt((abierta_rep1_GLU287A$`r_287_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_GLU287A$`r_287_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_GLU287A$`r_287_&_C-alpha Z` - abierta_rep1_protein$z)^2),
                           "glu287b" = sqrt((abierta_rep1_GLU287B$`r_484_&_C-alpha X` - abierta_rep1_protein$x)^2 + 
                                              (abierta_rep1_GLU287B$`r_484_&_C-alpha Y` - abierta_rep1_protein$y)^2 + 
                                              (abierta_rep1_GLU287B$`r_484_&_C-alpha Z` - abierta_rep1_protein$z)^2)
                            )

norm_abierta_rep1 <- data.frame("time" = abierta_rep1$time,
                                "gln278a" = (abierta_rep1$gln278a - abierta_rep1$gln278a[1]),
                                "gln278b" = (abierta_rep1$gln278b - abierta_rep1$gln278b[1]),
                                "asp279a" = (abierta_rep1$asp279a - abierta_rep1$asp279a[1]),
                                "asp279b" = (abierta_rep1$asp279b - abierta_rep1$asp279b[1]),
                                "asp280a" = (abierta_rep1$asp280a - abierta_rep1$asp280a[1]),
                                "asp280b" = (abierta_rep1$asp280b - abierta_rep1$asp280b[1]),
                                "cys281a" = (abierta_rep1$cys281a - abierta_rep1$cys281a[1]),
                                "cys281b" = (abierta_rep1$cys281b - abierta_rep1$cys281b[1]),
                                "ala282a" = (abierta_rep1$ala282a - abierta_rep1$ala282a[1]),
                                "ala282b" = (abierta_rep1$ala282b - abierta_rep1$ala282b[1]),
                                "ala283a" = (abierta_rep1$ala283a - abierta_rep1$ala283a[1]),
                                "ala283b" = (abierta_rep1$ala283b - abierta_rep1$ala283b[1]),
                                "phe284a" = (abierta_rep1$phe284a - abierta_rep1$phe284a[1]),
                                "phe284b" = (abierta_rep1$phe284b - abierta_rep1$phe284b[1]),
                                "ser285a" = (abierta_rep1$ser285a - abierta_rep1$ser285a[1]),
                                "ser285b" = (abierta_rep1$ser285b - abierta_rep1$ser285b[1]),
                                "arg286a" = (abierta_rep1$arg286a - abierta_rep1$arg286a[1]),
                                "arg286b" = (abierta_rep1$arg286b - abierta_rep1$arg286b[1]),
                                "glu287a" = (abierta_rep1$glu287a - abierta_rep1$glu287a[1]),
                                "glu287b" = (abierta_rep1$glu287b - abierta_rep1$glu287b[1])
                                )

rm(list=ls(pattern="abierta_rep1_"))



#### ABIERTA REP 2 ####
abierta_rep2_GLN278A <- readXVG("abierta-rep2/distance/abierta-rep2-GLN278A.xvg")
abierta_rep2_GLN278A <- as.data.frame(sapply(abierta_rep2_GLN278A, as.numeric))
abierta_rep2_GLN278B <- readXVG("abierta-rep2/distance/abierta-rep2-GLN278B.xvg")
abierta_rep2_GLN278B <- as.data.frame(sapply(abierta_rep2_GLN278B, as.numeric))
abierta_rep2_ASP279A <- readXVG("abierta-rep2/distance/abierta-rep2-ASP279A.xvg")
abierta_rep2_ASP279A <- as.data.frame(sapply(abierta_rep2_ASP279A, as.numeric))
abierta_rep2_ASP279B <- readXVG("abierta-rep2/distance/abierta-rep2-ASP279B.xvg")
abierta_rep2_ASP279B <- as.data.frame(sapply(abierta_rep2_ASP279B, as.numeric))
abierta_rep2_ASP280A <- readXVG("abierta-rep2/distance/abierta-rep2-ASP280A.xvg")
abierta_rep2_ASP280A <- as.data.frame(sapply(abierta_rep2_ASP280A, as.numeric))
abierta_rep2_ASP280B <- readXVG("abierta-rep2/distance/abierta-rep2-ASP280B.xvg")
abierta_rep2_ASP280B <- as.data.frame(sapply(abierta_rep2_ASP280B, as.numeric))
abierta_rep2_CYS281A <- readXVG("abierta-rep2/distance/abierta-rep2-CYS281A.xvg")
abierta_rep2_CYS281A <- as.data.frame(sapply(abierta_rep2_CYS281A, as.numeric))
abierta_rep2_CYS281B <- readXVG("abierta-rep2/distance/abierta-rep2-CYS281B.xvg")
abierta_rep2_CYS281B <- as.data.frame(sapply(abierta_rep2_CYS281B, as.numeric))
abierta_rep2_ALA282A <- readXVG("abierta-rep2/distance/abierta-rep2-ALA282A.xvg")
abierta_rep2_ALA282A <- as.data.frame(sapply(abierta_rep2_ALA282A, as.numeric))
abierta_rep2_ALA282B <- readXVG("abierta-rep2/distance/abierta-rep2-ALA282B.xvg")
abierta_rep2_ALA282B <- as.data.frame(sapply(abierta_rep2_ALA282B, as.numeric))
abierta_rep2_ALA283A <- readXVG("abierta-rep2/distance/abierta-rep2-ALA283A.xvg")
abierta_rep2_ALA283A <- as.data.frame(sapply(abierta_rep2_ALA283A, as.numeric))
abierta_rep2_ALA283B <- readXVG("abierta-rep2/distance/abierta-rep2-ALA283B.xvg")
abierta_rep2_ALA283B <- as.data.frame(sapply(abierta_rep2_ALA283B, as.numeric))
abierta_rep2_PHE284A <- readXVG("abierta-rep2/distance/abierta-rep2-PHE284A.xvg")
abierta_rep2_PHE284A <- as.data.frame(sapply(abierta_rep2_PHE284A, as.numeric))
abierta_rep2_PHE284B <- readXVG("abierta-rep2/distance/abierta-rep2-PHE284B.xvg")
abierta_rep2_PHE284B <- as.data.frame(sapply(abierta_rep2_PHE284B, as.numeric))
abierta_rep2_SER285A <- readXVG("abierta-rep2/distance/abierta-rep2-SER285A.xvg")
abierta_rep2_SER285A <- as.data.frame(sapply(abierta_rep2_SER285A, as.numeric))
abierta_rep2_SER285B <- readXVG("abierta-rep2/distance/abierta-rep2-SER285B.xvg")
abierta_rep2_SER285B <- as.data.frame(sapply(abierta_rep2_SER285B, as.numeric))
abierta_rep2_ARG286A <- readXVG("abierta-rep2/distance/abierta-rep2-ARG286A.xvg")
abierta_rep2_ARG286A <- as.data.frame(sapply(abierta_rep2_ARG286A, as.numeric))
abierta_rep2_ARG286B <- readXVG("abierta-rep2/distance/abierta-rep2-ARG286B.xvg")
abierta_rep2_ARG286B <- as.data.frame(sapply(abierta_rep2_ARG286B, as.numeric))
abierta_rep2_GLU287A <- readXVG("abierta-rep2/distance/abierta-rep2-GLU287A.xvg")
abierta_rep2_GLU287A <- as.data.frame(sapply(abierta_rep2_GLU287A, as.numeric))
abierta_rep2_GLU287B <- readXVG("abierta-rep2/distance/abierta-rep2-GLU287B.xvg")
abierta_rep2_GLU287B <- as.data.frame(sapply(abierta_rep2_GLU287B, as.numeric))
abierta_rep2_protein <- read.delim("abierta-rep2/distance/abierta2_protein.xvg")


abierta_rep2 <- data.frame('time' = abierta_rep2_protein$time/1000,
                           "gln278a" = sqrt((abierta_rep2_GLN278A$`r_278_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_GLN278A$`r_278_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_GLN278A$`r_278_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "gln278b" = sqrt((abierta_rep2_GLN278B$`r_475_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_GLN278B$`r_475_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_GLN278B$`r_475_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "asp279a" = sqrt((abierta_rep2_ASP279A$`r_279_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ASP279A$`r_279_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ASP279A$`r_279_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "asp279b" = sqrt((abierta_rep2_ASP279B$`r_476_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ASP279B$`r_476_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ASP279B$`r_476_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "asp280a" = sqrt((abierta_rep2_ASP280A$`r_280_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ASP280A$`r_280_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ASP280A$`r_280_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "asp280b" = sqrt((abierta_rep2_ASP280B$`r_477_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ASP280B$`r_477_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ASP280B$`r_477_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "cys281a" = sqrt((abierta_rep2_CYS281A$`r_281_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_CYS281A$`r_281_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_CYS281A$`r_281_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "cys281b" = sqrt((abierta_rep2_CYS281B$`r_478_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_CYS281B$`r_478_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_CYS281B$`r_478_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ala282a" = sqrt((abierta_rep2_ALA282A$`r_282_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ALA282A$`r_282_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ALA282A$`r_282_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ala282b" = sqrt((abierta_rep2_ALA282B$`r_479_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ALA282B$`r_479_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ALA282B$`r_479_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ala283a" = sqrt((abierta_rep2_ALA283A$`r_283_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ALA283A$`r_283_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ALA283A$`r_283_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ala283b" = sqrt((abierta_rep2_ALA283B$`r_480_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ALA283B$`r_480_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ALA283B$`r_480_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "phe284a" = sqrt((abierta_rep2_PHE284A$`r_284_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_PHE284A$`r_284_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_PHE284A$`r_284_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "phe284b" = sqrt((abierta_rep2_PHE284B$`r_481_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_PHE284B$`r_481_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_PHE284B$`r_481_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ser285a" = sqrt((abierta_rep2_SER285A$`r_285_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_SER285A$`r_285_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_SER285A$`r_285_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "ser285b" = sqrt((abierta_rep2_SER285B$`r_482_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_SER285B$`r_482_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_SER285B$`r_482_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "arg286a" = sqrt((abierta_rep2_ARG286A$`r_286_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ARG286A$`r_286_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ARG286A$`r_286_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "arg286b" = sqrt((abierta_rep2_ARG286B$`r_483_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_ARG286B$`r_483_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_ARG286B$`r_483_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "glu287a" = sqrt((abierta_rep2_GLU287A$`r_287_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_GLU287A$`r_287_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_GLU287A$`r_287_&_C-alpha Z` - abierta_rep2_protein$z)^2),
                           "glu287b" = sqrt((abierta_rep2_GLU287B$`r_484_&_C-alpha X` - abierta_rep2_protein$x)^2 + 
                                              (abierta_rep2_GLU287B$`r_484_&_C-alpha Y` - abierta_rep2_protein$y)^2 + 
                                              (abierta_rep2_GLU287B$`r_484_&_C-alpha Z` - abierta_rep2_protein$z)^2)
)

norm_abierta_rep2 <- data.frame("time" = abierta_rep2$time,
                                "gln278a" = (abierta_rep2$gln278a - abierta_rep2$gln278a[1]),
                                "gln278b" = (abierta_rep2$gln278b - abierta_rep2$gln278b[1]),
                                "asp279a" = (abierta_rep2$asp279a - abierta_rep2$asp279a[1]),
                                "asp279b" = (abierta_rep2$asp279b - abierta_rep2$asp279b[1]),
                                "asp280a" = (abierta_rep2$asp280a - abierta_rep2$asp280a[1]),
                                "asp280b" = (abierta_rep2$asp280b - abierta_rep2$asp280b[1]),
                                "cys281a" = (abierta_rep2$cys281a - abierta_rep2$cys281a[1]),
                                "cys281b" = (abierta_rep2$cys281b - abierta_rep2$cys281b[1]),
                                "ala282a" = (abierta_rep2$ala282a - abierta_rep2$ala282a[1]),
                                "ala282b" = (abierta_rep2$ala282b - abierta_rep2$ala282b[1]),
                                "ala283a" = (abierta_rep2$ala283a - abierta_rep2$ala283a[1]),
                                "ala283b" = (abierta_rep2$ala283b - abierta_rep2$ala283b[1]),
                                "phe284a" = (abierta_rep2$phe284a - abierta_rep2$phe284a[1]),
                                "phe284b" = (abierta_rep2$phe284b - abierta_rep2$phe284b[1]),
                                "ser285a" = (abierta_rep2$ser285a - abierta_rep2$ser285a[1]),
                                "ser285b" = (abierta_rep2$ser285b - abierta_rep2$ser285b[1]),
                                "arg286a" = (abierta_rep2$arg286a - abierta_rep2$arg286a[1]),
                                "arg286b" = (abierta_rep2$arg286b - abierta_rep2$arg286b[1]),
                                "glu287a" = (abierta_rep2$glu287a - abierta_rep2$glu287a[1]),
                                "glu287b" = (abierta_rep2$glu287b - abierta_rep2$glu287b[1])
)

rm(list=ls(pattern="abierta_rep2_"))
#### CERRADA REP 1 ####
cerrada_rep1_GLN278A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-GLN278A.xvg")
cerrada_rep1_GLN278A <- as.data.frame(sapply(cerrada_rep1_GLN278A, as.numeric))
cerrada_rep1_GLN278B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-GLN278B.xvg")
cerrada_rep1_GLN278B <- as.data.frame(sapply(cerrada_rep1_GLN278B, as.numeric))
cerrada_rep1_ASP279A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ASP279A.xvg")
cerrada_rep1_ASP279A <- as.data.frame(sapply(cerrada_rep1_ASP279A, as.numeric))
cerrada_rep1_ASP279B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ASP279B.xvg")
cerrada_rep1_ASP279B <- as.data.frame(sapply(cerrada_rep1_ASP279B, as.numeric))
cerrada_rep1_ASP280A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ASP280A.xvg")
cerrada_rep1_ASP280A <- as.data.frame(sapply(cerrada_rep1_ASP280A, as.numeric))
cerrada_rep1_ASP280B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ASP280B.xvg")
cerrada_rep1_ASP280B <- as.data.frame(sapply(cerrada_rep1_ASP280B, as.numeric))
cerrada_rep1_CYS281A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-CYS281A.xvg")
cerrada_rep1_CYS281A <- as.data.frame(sapply(cerrada_rep1_CYS281A, as.numeric))
cerrada_rep1_CYS281B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-CYS281B.xvg")
cerrada_rep1_CYS281B <- as.data.frame(sapply(cerrada_rep1_CYS281B, as.numeric))
cerrada_rep1_ALA282A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ALA282A.xvg")
cerrada_rep1_ALA282A <- as.data.frame(sapply(cerrada_rep1_ALA282A, as.numeric))
cerrada_rep1_ALA282B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ALA282B.xvg")
cerrada_rep1_ALA282B <- as.data.frame(sapply(cerrada_rep1_ALA282B, as.numeric))
cerrada_rep1_ALA283A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ALA283A.xvg")
cerrada_rep1_ALA283A <- as.data.frame(sapply(cerrada_rep1_ALA283A, as.numeric))
cerrada_rep1_ALA283B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ALA283B.xvg")
cerrada_rep1_ALA283B <- as.data.frame(sapply(cerrada_rep1_ALA283B, as.numeric))
cerrada_rep1_PHE284A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-PHE284A.xvg")
cerrada_rep1_PHE284A <- as.data.frame(sapply(cerrada_rep1_PHE284A, as.numeric))
cerrada_rep1_PHE284B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-PHE284B.xvg")
cerrada_rep1_PHE284B <- as.data.frame(sapply(cerrada_rep1_PHE284B, as.numeric))
cerrada_rep1_SER285A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-SER285A.xvg")
cerrada_rep1_SER285A <- as.data.frame(sapply(cerrada_rep1_SER285A, as.numeric))
cerrada_rep1_SER285B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-SER285B.xvg")
cerrada_rep1_SER285B <- as.data.frame(sapply(cerrada_rep1_SER285B, as.numeric))
cerrada_rep1_ARG286A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ARG286A.xvg")
cerrada_rep1_ARG286A <- as.data.frame(sapply(cerrada_rep1_ARG286A, as.numeric))
cerrada_rep1_ARG286B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-ARG286B.xvg")
cerrada_rep1_ARG286B <- as.data.frame(sapply(cerrada_rep1_ARG286B, as.numeric))
cerrada_rep1_GLU287A <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-GLU287A.xvg")
cerrada_rep1_GLU287A <- as.data.frame(sapply(cerrada_rep1_GLU287A, as.numeric))
cerrada_rep1_GLU287B <- readXVG("cerrada-rep1/prod/distances/cerrada-rep1-GLU287B.xvg")
cerrada_rep1_GLU287B <- as.data.frame(sapply(cerrada_rep1_GLU287B, as.numeric))
cerrada_rep1_protein <- read.delim("cerrada-rep1/prod/distances/protein-rep1.xvg")


cerrada_rep1 <- data.frame('time' = cerrada_rep1_protein$time/1000,
                           "gln278a" = sqrt((cerrada_rep1_GLN278A$`r_278_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_GLN278A$`r_278_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_GLN278A$`r_278_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "gln278b" = sqrt((cerrada_rep1_GLN278B$`r_475_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_GLN278B$`r_475_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_GLN278B$`r_475_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "asp279a" = sqrt((cerrada_rep1_ASP279A$`r_279_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ASP279A$`r_279_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ASP279A$`r_279_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "asp279b" = sqrt((cerrada_rep1_ASP279B$`r_476_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ASP279B$`r_476_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ASP279B$`r_476_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "asp280a" = sqrt((cerrada_rep1_ASP280A$`r_280_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ASP280A$`r_280_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ASP280A$`r_280_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "asp280b" = sqrt((cerrada_rep1_ASP280B$`r_477_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ASP280B$`r_477_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ASP280B$`r_477_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "cys281a" = sqrt((cerrada_rep1_CYS281A$`r_281_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_CYS281A$`r_281_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_CYS281A$`r_281_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "cys281b" = sqrt((cerrada_rep1_CYS281B$`r_478_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_CYS281B$`r_478_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_CYS281B$`r_478_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ala282a" = sqrt((cerrada_rep1_ALA282A$`r_282_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ALA282A$`r_282_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ALA282A$`r_282_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ala282b" = sqrt((cerrada_rep1_ALA282B$`r_479_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ALA282B$`r_479_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ALA282B$`r_479_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ala283a" = sqrt((cerrada_rep1_ALA283A$`r_283_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ALA283A$`r_283_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ALA283A$`r_283_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ala283b" = sqrt((cerrada_rep1_ALA283B$`r_480_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ALA283B$`r_480_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ALA283B$`r_480_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "phe284a" = sqrt((cerrada_rep1_PHE284A$`r_284_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_PHE284A$`r_284_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_PHE284A$`r_284_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "phe284b" = sqrt((cerrada_rep1_PHE284B$`r_481_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_PHE284B$`r_481_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_PHE284B$`r_481_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ser285a" = sqrt((cerrada_rep1_SER285A$`r_285_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_SER285A$`r_285_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_SER285A$`r_285_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "ser285b" = sqrt((cerrada_rep1_SER285B$`r_482_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_SER285B$`r_482_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_SER285B$`r_482_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "arg286a" = sqrt((cerrada_rep1_ARG286A$`r_286_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ARG286A$`r_286_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ARG286A$`r_286_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "arg286b" = sqrt((cerrada_rep1_ARG286B$`r_483_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_ARG286B$`r_483_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_ARG286B$`r_483_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "glu287a" = sqrt((cerrada_rep1_GLU287A$`r_287_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_GLU287A$`r_287_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_GLU287A$`r_287_&_C-alpha Z` - cerrada_rep1_protein$z)^2),
                           "glu287b" = sqrt((cerrada_rep1_GLU287B$`r_484_&_C-alpha X` - cerrada_rep1_protein$x)^2 + 
                                              (cerrada_rep1_GLU287B$`r_484_&_C-alpha Y` - cerrada_rep1_protein$y)^2 + 
                                              (cerrada_rep1_GLU287B$`r_484_&_C-alpha Z` - cerrada_rep1_protein$z)^2)
)

norm_cerrada_rep1 <- data.frame("time" = cerrada_rep1$time,
                                "gln278a" = (cerrada_rep1$gln278a - cerrada_rep1$gln278a[1]),
                                "gln278b" = (cerrada_rep1$gln278b - cerrada_rep1$gln278b[1]),
                                "asp279a" = (cerrada_rep1$asp279a - cerrada_rep1$asp279a[1]),
                                "asp279b" = (cerrada_rep1$asp279b - cerrada_rep1$asp279b[1]),
                                "asp280a" = (cerrada_rep1$asp280a - cerrada_rep1$asp280a[1]),
                                "asp280b" = (cerrada_rep1$asp280b - cerrada_rep1$asp280b[1]),
                                "cys281a" = (cerrada_rep1$cys281a - cerrada_rep1$cys281a[1]),
                                "cys281b" = (cerrada_rep1$cys281b - cerrada_rep1$cys281b[1]),
                                "ala282a" = (cerrada_rep1$ala282a - cerrada_rep1$ala282a[1]),
                                "ala282b" = (cerrada_rep1$ala282b - cerrada_rep1$ala282b[1]),
                                "ala283a" = (cerrada_rep1$ala283a - cerrada_rep1$ala283a[1]),
                                "ala283b" = (cerrada_rep1$ala283b - cerrada_rep1$ala283b[1]),
                                "phe284a" = (cerrada_rep1$phe284a - cerrada_rep1$phe284a[1]),
                                "phe284b" = (cerrada_rep1$phe284b - cerrada_rep1$phe284b[1]),
                                "ser285a" = (cerrada_rep1$ser285a - cerrada_rep1$ser285a[1]),
                                "ser285b" = (cerrada_rep1$ser285b - cerrada_rep1$ser285b[1]),
                                "arg286a" = (cerrada_rep1$arg286a - cerrada_rep1$arg286a[1]),
                                "arg286b" = (cerrada_rep1$arg286b - cerrada_rep1$arg286b[1]),
                                "glu287a" = (cerrada_rep1$glu287a - cerrada_rep1$glu287a[1]),
                                "glu287b" = (cerrada_rep1$glu287b - cerrada_rep1$glu287b[1])
)

rm(list=ls(pattern="cerrada_rep1_"))
#### CERRADA REP 2 ####
cerrada_rep2_GLN278A <- readXVG("cerrada-rep2/distances/cerrada-rep2-GLN278A.xvg")
cerrada_rep2_GLN278A <- as.data.frame(sapply(cerrada_rep2_GLN278A, as.numeric))
cerrada_rep2_GLN278B <- readXVG("cerrada-rep2/distances/cerrada-rep2-GLN278B.xvg")
cerrada_rep2_GLN278B <- as.data.frame(sapply(cerrada_rep2_GLN278B, as.numeric))
cerrada_rep2_ASP279A <- readXVG("cerrada-rep2/distances/cerrada-rep2-ASP279A.xvg")
cerrada_rep2_ASP279A <- as.data.frame(sapply(cerrada_rep2_ASP279A, as.numeric))
cerrada_rep2_ASP279B <- readXVG("cerrada-rep2/distances/cerrada-rep2-ASP279B.xvg")
cerrada_rep2_ASP279B <- as.data.frame(sapply(cerrada_rep2_ASP279B, as.numeric))
cerrada_rep2_ASP280A <- readXVG("cerrada-rep2/distances/cerrada-rep2-ASP280A.xvg")
cerrada_rep2_ASP280A <- as.data.frame(sapply(cerrada_rep2_ASP280A, as.numeric))
cerrada_rep2_ASP280B <- readXVG("cerrada-rep2/distances/cerrada-rep2-ASP280B.xvg")
cerrada_rep2_ASP280B <- as.data.frame(sapply(cerrada_rep2_ASP280B, as.numeric))
cerrada_rep2_CYS281A <- readXVG("cerrada-rep2/distances/cerrada-rep2-CYS281A.xvg")
cerrada_rep2_CYS281A <- as.data.frame(sapply(cerrada_rep2_CYS281A, as.numeric))
cerrada_rep2_CYS281B <- readXVG("cerrada-rep2/distances/cerrada-rep2-CYS281B.xvg")
cerrada_rep2_CYS281B <- as.data.frame(sapply(cerrada_rep2_CYS281B, as.numeric))
cerrada_rep2_ALA282A <- readXVG("cerrada-rep2/distances/cerrada-rep2-ALA282A.xvg")
cerrada_rep2_ALA282A <- as.data.frame(sapply(cerrada_rep2_ALA282A, as.numeric))
cerrada_rep2_ALA282B <- readXVG("cerrada-rep2/distances/cerrada-rep2-ALA282B.xvg")
cerrada_rep2_ALA282B <- as.data.frame(sapply(cerrada_rep2_ALA282B, as.numeric))
cerrada_rep2_ALA283A <- readXVG("cerrada-rep2/distances/cerrada-rep2-ALA283A.xvg")
cerrada_rep2_ALA283A <- as.data.frame(sapply(cerrada_rep2_ALA283A, as.numeric))
cerrada_rep2_ALA283B <- readXVG("cerrada-rep2/distances/cerrada-rep2-ALA283B.xvg")
cerrada_rep2_ALA283B <- as.data.frame(sapply(cerrada_rep2_ALA283B, as.numeric))
cerrada_rep2_PHE284A <- readXVG("cerrada-rep2/distances/cerrada-rep2-PHE284A.xvg")
cerrada_rep2_PHE284A <- as.data.frame(sapply(cerrada_rep2_PHE284A, as.numeric))
cerrada_rep2_PHE284B <- readXVG("cerrada-rep2/distances/cerrada-rep2-PHE284B.xvg")
cerrada_rep2_PHE284B <- as.data.frame(sapply(cerrada_rep2_PHE284B, as.numeric))
cerrada_rep2_SER285A <- readXVG("cerrada-rep2/distances/cerrada-rep2-SER285A.xvg")
cerrada_rep2_SER285A <- as.data.frame(sapply(cerrada_rep2_SER285A, as.numeric))
cerrada_rep2_SER285B <- readXVG("cerrada-rep2/distances/cerrada-rep2-SER285B.xvg")
cerrada_rep2_SER285B <- as.data.frame(sapply(cerrada_rep2_SER285B, as.numeric))
cerrada_rep2_ARG286A <- readXVG("cerrada-rep2/distances/cerrada-rep2-ARG286A.xvg")
cerrada_rep2_ARG286A <- as.data.frame(sapply(cerrada_rep2_ARG286A, as.numeric))
cerrada_rep2_ARG286B <- readXVG("cerrada-rep2/distances/cerrada-rep2-ARG286B.xvg")
cerrada_rep2_ARG286B <- as.data.frame(sapply(cerrada_rep2_ARG286B, as.numeric))
cerrada_rep2_GLU287A <- readXVG("cerrada-rep2/distances/cerrada-rep2-GLU287A.xvg")
cerrada_rep2_GLU287A <- as.data.frame(sapply(cerrada_rep2_GLU287A, as.numeric))
cerrada_rep2_GLU287B <- readXVG("cerrada-rep2/distances/cerrada-rep2-GLU287B.xvg")
cerrada_rep2_GLU287B <- as.data.frame(sapply(cerrada_rep2_GLU287B, as.numeric))
cerrada_rep2_protein <- read.delim("cerrada-rep1/prod/distances/protein-rep1.xvg")


cerrada_rep2 <- data.frame('time' = cerrada_rep2_protein$time/1000,
                           "gln278a" = sqrt((cerrada_rep2_GLN278A$`r_278_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_GLN278A$`r_278_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_GLN278A$`r_278_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "gln278b" = sqrt((cerrada_rep2_GLN278B$`r_475_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_GLN278B$`r_475_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_GLN278B$`r_475_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "asp279a" = sqrt((cerrada_rep2_ASP279A$`r_279_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ASP279A$`r_279_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ASP279A$`r_279_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "asp279b" = sqrt((cerrada_rep2_ASP279B$`r_476_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ASP279B$`r_476_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ASP279B$`r_476_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "asp280a" = sqrt((cerrada_rep2_ASP280A$`r_280_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ASP280A$`r_280_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ASP280A$`r_280_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "asp280b" = sqrt((cerrada_rep2_ASP280B$`r_477_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ASP280B$`r_477_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ASP280B$`r_477_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "cys281a" = sqrt((cerrada_rep2_CYS281A$`r_281_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_CYS281A$`r_281_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_CYS281A$`r_281_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "cys281b" = sqrt((cerrada_rep2_CYS281B$`r_478_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_CYS281B$`r_478_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_CYS281B$`r_478_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ala282a" = sqrt((cerrada_rep2_ALA282A$`r_282_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ALA282A$`r_282_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ALA282A$`r_282_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ala282b" = sqrt((cerrada_rep2_ALA282B$`r_479_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ALA282B$`r_479_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ALA282B$`r_479_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ala283a" = sqrt((cerrada_rep2_ALA283A$`r_283_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ALA283A$`r_283_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ALA283A$`r_283_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ala283b" = sqrt((cerrada_rep2_ALA283B$`r_480_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ALA283B$`r_480_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ALA283B$`r_480_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "phe284a" = sqrt((cerrada_rep2_PHE284A$`r_284_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_PHE284A$`r_284_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_PHE284A$`r_284_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "phe284b" = sqrt((cerrada_rep2_PHE284B$`r_481_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_PHE284B$`r_481_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_PHE284B$`r_481_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ser285a" = sqrt((cerrada_rep2_SER285A$`r_285_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_SER285A$`r_285_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_SER285A$`r_285_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "ser285b" = sqrt((cerrada_rep2_SER285B$`r_482_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_SER285B$`r_482_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_SER285B$`r_482_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "arg286a" = sqrt((cerrada_rep2_ARG286A$`r_286_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ARG286A$`r_286_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ARG286A$`r_286_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "arg286b" = sqrt((cerrada_rep2_ARG286B$`r_483_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_ARG286B$`r_483_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_ARG286B$`r_483_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "glu287a" = sqrt((cerrada_rep2_GLU287A$`r_287_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_GLU287A$`r_287_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_GLU287A$`r_287_&_C-alpha Z` - cerrada_rep2_protein$z)^2),
                           "glu287b" = sqrt((cerrada_rep2_GLU287B$`r_484_&_C-alpha X` - cerrada_rep2_protein$x)^2 + 
                                              (cerrada_rep2_GLU287B$`r_484_&_C-alpha Y` - cerrada_rep2_protein$y)^2 + 
                                              (cerrada_rep2_GLU287B$`r_484_&_C-alpha Z` - cerrada_rep2_protein$z)^2)
)

norm_cerrada_rep2 <- data.frame("time" = cerrada_rep2$time,
                                "gln278a" = (cerrada_rep2$gln278a - cerrada_rep2$gln278a[1]),
                                "gln278b" = (cerrada_rep2$gln278b - cerrada_rep2$gln278b[1]),
                                "asp279a" = (cerrada_rep2$asp279a - cerrada_rep2$asp279a[1]),
                                "asp279b" = (cerrada_rep2$asp279b - cerrada_rep2$asp279b[1]),
                                "asp280a" = (cerrada_rep2$asp280a - cerrada_rep2$asp280a[1]),
                                "asp280b" = (cerrada_rep2$asp280b - cerrada_rep2$asp280b[1]),
                                "cys281a" = (cerrada_rep2$cys281a - cerrada_rep2$cys281a[1]),
                                "cys281b" = (cerrada_rep2$cys281b - cerrada_rep2$cys281b[1]),
                                "ala282a" = (cerrada_rep2$ala282a - cerrada_rep2$ala282a[1]),
                                "ala282b" = (cerrada_rep2$ala282b - cerrada_rep2$ala282b[1]),
                                "ala283a" = (cerrada_rep2$ala283a - cerrada_rep2$ala283a[1]),
                                "ala283b" = (cerrada_rep2$ala283b - cerrada_rep2$ala283b[1]),
                                "phe284a" = (cerrada_rep2$phe284a - cerrada_rep2$phe284a[1]),
                                "phe284b" = (cerrada_rep2$phe284b - cerrada_rep2$phe284b[1]),
                                "ser285a" = (cerrada_rep2$ser285a - cerrada_rep2$ser285a[1]),
                                "ser285b" = (cerrada_rep2$ser285b - cerrada_rep2$ser285b[1]),
                                "arg286a" = (cerrada_rep2$arg286a - cerrada_rep2$arg286a[1]),
                                "arg286b" = (cerrada_rep2$arg286b - cerrada_rep2$arg286b[1]),
                                "glu287a" = (cerrada_rep2$glu287a - cerrada_rep2$glu287a[1]),
                                "glu287b" = (cerrada_rep2$glu287b - cerrada_rep2$glu287b[1])
)

rm(list=ls(pattern="cerrada_rep2_"))
#### V160M REP 1 ####
v160m_rep1_GLN278A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-GLN278A.xvg")
v160m_rep1_GLN278A <- as.data.frame(sapply(v160m_rep1_GLN278A, as.numeric))
v160m_rep1_GLN278B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-GLN278B.xvg")
v160m_rep1_GLN278B <- as.data.frame(sapply(v160m_rep1_GLN278B, as.numeric))
v160m_rep1_ASP279A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ASP279A.xvg")
v160m_rep1_ASP279A <- as.data.frame(sapply(v160m_rep1_ASP279A, as.numeric))
v160m_rep1_ASP279B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ASP279B.xvg")
v160m_rep1_ASP279B <- as.data.frame(sapply(v160m_rep1_ASP279B, as.numeric))
v160m_rep1_ASP280A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ASP280A.xvg")
v160m_rep1_ASP280A <- as.data.frame(sapply(v160m_rep1_ASP280A, as.numeric))
v160m_rep1_ASP280B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ASP280B.xvg")
v160m_rep1_ASP280B <- as.data.frame(sapply(v160m_rep1_ASP280B, as.numeric))
v160m_rep1_CYS281A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-CYS281A.xvg")
v160m_rep1_CYS281A <- as.data.frame(sapply(v160m_rep1_CYS281A, as.numeric))
v160m_rep1_CYS281B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-CYS281B.xvg")
v160m_rep1_CYS281B <- as.data.frame(sapply(v160m_rep1_CYS281B, as.numeric))
v160m_rep1_ALA282A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ALA282A.xvg")
v160m_rep1_ALA282A <- as.data.frame(sapply(v160m_rep1_ALA282A, as.numeric))
v160m_rep1_ALA282B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ALA282B.xvg")
v160m_rep1_ALA282B <- as.data.frame(sapply(v160m_rep1_ALA282B, as.numeric))
v160m_rep1_ALA283A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ALA283A.xvg")
v160m_rep1_ALA283A <- as.data.frame(sapply(v160m_rep1_ALA283A, as.numeric))
v160m_rep1_ALA283B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ALA283B.xvg")
v160m_rep1_ALA283B <- as.data.frame(sapply(v160m_rep1_ALA283B, as.numeric))
v160m_rep1_PHE284A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-PHE284A.xvg")
v160m_rep1_PHE284A <- as.data.frame(sapply(v160m_rep1_PHE284A, as.numeric))
v160m_rep1_PHE284B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-PHE284B.xvg")
v160m_rep1_PHE284B <- as.data.frame(sapply(v160m_rep1_PHE284B, as.numeric))
v160m_rep1_SER285A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-SER285A.xvg")
v160m_rep1_SER285A <- as.data.frame(sapply(v160m_rep1_SER285A, as.numeric))
v160m_rep1_SER285B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-SER285B.xvg")
v160m_rep1_SER285B <- as.data.frame(sapply(v160m_rep1_SER285B, as.numeric))
v160m_rep1_ARG286A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ARG286A.xvg")
v160m_rep1_ARG286A <- as.data.frame(sapply(v160m_rep1_ARG286A, as.numeric))
v160m_rep1_ARG286B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-ARG286B.xvg")
v160m_rep1_ARG286B <- as.data.frame(sapply(v160m_rep1_ARG286B, as.numeric))
v160m_rep1_GLU287A <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-GLU287A.xvg")
v160m_rep1_GLU287A <- as.data.frame(sapply(v160m_rep1_GLU287A, as.numeric))
v160m_rep1_GLU287B <- readXVG("v160m-rep1/prod_200ns/distance/abierta-v160m-rep1-GLU287B.xvg")
v160m_rep1_GLU287B <- as.data.frame(sapply(v160m_rep1_GLU287B, as.numeric))
v160m_rep1_protein <- read.delim("v160m-rep1/prod_200ns/distance/v160m1_protein.xvg")


v160m_rep1 <- data.frame('time' = v160m_rep1_protein$time/1000,
                         "gln278a" = sqrt((v160m_rep1_GLN278A$`r_278_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_GLN278A$`r_278_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_GLN278A$`r_278_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "gln278b" = sqrt((v160m_rep1_GLN278B$`r_475_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_GLN278B$`r_475_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_GLN278B$`r_475_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "asp279a" = sqrt((v160m_rep1_ASP279A$`r_279_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ASP279A$`r_279_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ASP279A$`r_279_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "asp279b" = sqrt((v160m_rep1_ASP279B$`r_476_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ASP279B$`r_476_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ASP279B$`r_476_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "asp280a" = sqrt((v160m_rep1_ASP280A$`r_280_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ASP280A$`r_280_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ASP280A$`r_280_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "asp280b" = sqrt((v160m_rep1_ASP280B$`r_477_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ASP280B$`r_477_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ASP280B$`r_477_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "cys281a" = sqrt((v160m_rep1_CYS281A$`r_281_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_CYS281A$`r_281_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_CYS281A$`r_281_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "cys281b" = sqrt((v160m_rep1_CYS281B$`r_478_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_CYS281B$`r_478_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_CYS281B$`r_478_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ala282a" = sqrt((v160m_rep1_ALA282A$`r_282_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ALA282A$`r_282_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ALA282A$`r_282_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ala282b" = sqrt((v160m_rep1_ALA282B$`r_479_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ALA282B$`r_479_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ALA282B$`r_479_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ala283a" = sqrt((v160m_rep1_ALA283A$`r_283_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ALA283A$`r_283_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ALA283A$`r_283_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ala283b" = sqrt((v160m_rep1_ALA283B$`r_480_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ALA283B$`r_480_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ALA283B$`r_480_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "phe284a" = sqrt((v160m_rep1_PHE284A$`r_284_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_PHE284A$`r_284_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_PHE284A$`r_284_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "phe284b" = sqrt((v160m_rep1_PHE284B$`r_481_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_PHE284B$`r_481_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_PHE284B$`r_481_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ser285a" = sqrt((v160m_rep1_SER285A$`r_285_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_SER285A$`r_285_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_SER285A$`r_285_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "ser285b" = sqrt((v160m_rep1_SER285B$`r_482_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_SER285B$`r_482_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_SER285B$`r_482_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "arg286a" = sqrt((v160m_rep1_ARG286A$`r_286_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ARG286A$`r_286_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ARG286A$`r_286_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "arg286b" = sqrt((v160m_rep1_ARG286B$`r_483_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_ARG286B$`r_483_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_ARG286B$`r_483_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "glu287a" = sqrt((v160m_rep1_GLU287A$`r_287_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_GLU287A$`r_287_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_GLU287A$`r_287_&_C-alpha Z` - v160m_rep1_protein$z)^2),
                         "glu287b" = sqrt((v160m_rep1_GLU287B$`r_484_&_C-alpha X` - v160m_rep1_protein$x)^2 + 
                                            (v160m_rep1_GLU287B$`r_484_&_C-alpha Y` - v160m_rep1_protein$y)^2 + 
                                            (v160m_rep1_GLU287B$`r_484_&_C-alpha Z` - v160m_rep1_protein$z)^2)
)

norm_v160m_rep1 <- data.frame("time" = v160m_rep1$time,
                              "gln278a" = (v160m_rep1$gln278a - v160m_rep1$gln278a[1]),
                              "gln278b" = (v160m_rep1$gln278b - v160m_rep1$gln278b[1]),
                              "asp279a" = (v160m_rep1$asp279a - v160m_rep1$asp279a[1]),
                              "asp279b" = (v160m_rep1$asp279b - v160m_rep1$asp279b[1]),
                              "asp280a" = (v160m_rep1$asp280a - v160m_rep1$asp280a[1]),
                              "asp280b" = (v160m_rep1$asp280b - v160m_rep1$asp280b[1]),
                              "cys281a" = (v160m_rep1$cys281a - v160m_rep1$cys281a[1]),
                              "cys281b" = (v160m_rep1$cys281b - v160m_rep1$cys281b[1]),
                              "ala282a" = (v160m_rep1$ala282a - v160m_rep1$ala282a[1]),
                              "ala282b" = (v160m_rep1$ala282b - v160m_rep1$ala282b[1]),
                              "ala283a" = (v160m_rep1$ala283a - v160m_rep1$ala283a[1]),
                              "ala283b" = (v160m_rep1$ala283b - v160m_rep1$ala283b[1]),
                              "phe284a" = (v160m_rep1$phe284a - v160m_rep1$phe284a[1]),
                              "phe284b" = (v160m_rep1$phe284b - v160m_rep1$phe284b[1]),
                              "ser285a" = (v160m_rep1$ser285a - v160m_rep1$ser285a[1]),
                              "ser285b" = (v160m_rep1$ser285b - v160m_rep1$ser285b[1]),
                              "arg286a" = (v160m_rep1$arg286a - v160m_rep1$arg286a[1]),
                              "arg286b" = (v160m_rep1$arg286b - v160m_rep1$arg286b[1]),
                              "glu287a" = (v160m_rep1$glu287a - v160m_rep1$glu287a[1]),
                              "glu287b" = (v160m_rep1$glu287b - v160m_rep1$glu287b[1])
)

rm(list=ls(pattern="v160m_rep1_"))
#### V160M REP 2 ####
v160m_rep2_GLN278A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-GLN278A.xvg")
v160m_rep2_GLN278A <- as.data.frame(sapply(v160m_rep2_GLN278A, as.numeric))
v160m_rep2_GLN278B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-GLN278B.xvg")
v160m_rep2_GLN278B <- as.data.frame(sapply(v160m_rep2_GLN278B, as.numeric))
v160m_rep2_ASP279A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ASP279A.xvg")
v160m_rep2_ASP279A <- as.data.frame(sapply(v160m_rep2_ASP279A, as.numeric))
v160m_rep2_ASP279B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ASP279B.xvg")
v160m_rep2_ASP279B <- as.data.frame(sapply(v160m_rep2_ASP279B, as.numeric))
v160m_rep2_ASP280A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ASP280A.xvg")
v160m_rep2_ASP280A <- as.data.frame(sapply(v160m_rep2_ASP280A, as.numeric))
v160m_rep2_ASP280B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ASP280B.xvg")
v160m_rep2_ASP280B <- as.data.frame(sapply(v160m_rep2_ASP280B, as.numeric))
v160m_rep2_CYS281A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-CYS281A.xvg")
v160m_rep2_CYS281A <- as.data.frame(sapply(v160m_rep2_CYS281A, as.numeric))
v160m_rep2_CYS281B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-CYS281B.xvg")
v160m_rep2_CYS281B <- as.data.frame(sapply(v160m_rep2_CYS281B, as.numeric))
v160m_rep2_ALA282A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ALA282A.xvg")
v160m_rep2_ALA282A <- as.data.frame(sapply(v160m_rep2_ALA282A, as.numeric))
v160m_rep2_ALA282B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ALA282B.xvg")
v160m_rep2_ALA282B <- as.data.frame(sapply(v160m_rep2_ALA282B, as.numeric))
v160m_rep2_ALA283A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ALA283A.xvg")
v160m_rep2_ALA283A <- as.data.frame(sapply(v160m_rep2_ALA283A, as.numeric))
v160m_rep2_ALA283B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ALA283B.xvg")
v160m_rep2_ALA283B <- as.data.frame(sapply(v160m_rep2_ALA283B, as.numeric))
v160m_rep2_PHE284A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-PHE284A.xvg")
v160m_rep2_PHE284A <- as.data.frame(sapply(v160m_rep2_PHE284A, as.numeric))
v160m_rep2_PHE284B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-PHE284B.xvg")
v160m_rep2_PHE284B <- as.data.frame(sapply(v160m_rep2_PHE284B, as.numeric))
v160m_rep2_SER285A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-SER285A.xvg")
v160m_rep2_SER285A <- as.data.frame(sapply(v160m_rep2_SER285A, as.numeric))
v160m_rep2_SER285B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-SER285B.xvg")
v160m_rep2_SER285B <- as.data.frame(sapply(v160m_rep2_SER285B, as.numeric))
v160m_rep2_ARG286A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ARG286A.xvg")
v160m_rep2_ARG286A <- as.data.frame(sapply(v160m_rep2_ARG286A, as.numeric))
v160m_rep2_ARG286B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-ARG286B.xvg")
v160m_rep2_ARG286B <- as.data.frame(sapply(v160m_rep2_ARG286B, as.numeric))
v160m_rep2_GLU287A <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-GLU287A.xvg")
v160m_rep2_GLU287A <- as.data.frame(sapply(v160m_rep2_GLU287A, as.numeric))
v160m_rep2_GLU287B <- readXVG("v160m-rep2/distance/abierta-v160m-rep2-GLU287B.xvg")
v160m_rep2_GLU287B <- as.data.frame(sapply(v160m_rep2_GLU287B, as.numeric))
v160m_rep2_protein <- read.delim("v160m-rep2/distance/abierta-v160m-rep2-protein.xvg")


v160m_rep2 <- data.frame('time' = v160m_rep2_protein$time/1000,
                         "gln278a" = sqrt((v160m_rep2_GLN278A$`r_278_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_GLN278A$`r_278_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_GLN278A$`r_278_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "gln278b" = sqrt((v160m_rep2_GLN278B$`r_475_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_GLN278B$`r_475_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_GLN278B$`r_475_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "asp279a" = sqrt((v160m_rep2_ASP279A$`r_279_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ASP279A$`r_279_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ASP279A$`r_279_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "asp279b" = sqrt((v160m_rep2_ASP279B$`r_476_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ASP279B$`r_476_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ASP279B$`r_476_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "asp280a" = sqrt((v160m_rep2_ASP280A$`r_280_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ASP280A$`r_280_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ASP280A$`r_280_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "asp280b" = sqrt((v160m_rep2_ASP280B$`r_477_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ASP280B$`r_477_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ASP280B$`r_477_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "cys281a" = sqrt((v160m_rep2_CYS281A$`r_281_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_CYS281A$`r_281_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_CYS281A$`r_281_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "cys281b" = sqrt((v160m_rep2_CYS281B$`r_478_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_CYS281B$`r_478_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_CYS281B$`r_478_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ala282a" = sqrt((v160m_rep2_ALA282A$`r_282_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ALA282A$`r_282_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ALA282A$`r_282_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ala282b" = sqrt((v160m_rep2_ALA282B$`r_479_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ALA282B$`r_479_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ALA282B$`r_479_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ala283a" = sqrt((v160m_rep2_ALA283A$`r_283_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ALA283A$`r_283_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ALA283A$`r_283_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ala283b" = sqrt((v160m_rep2_ALA283B$`r_480_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ALA283B$`r_480_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ALA283B$`r_480_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "phe284a" = sqrt((v160m_rep2_PHE284A$`r_284_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_PHE284A$`r_284_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_PHE284A$`r_284_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "phe284b" = sqrt((v160m_rep2_PHE284B$`r_481_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_PHE284B$`r_481_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_PHE284B$`r_481_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ser285a" = sqrt((v160m_rep2_SER285A$`r_285_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_SER285A$`r_285_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_SER285A$`r_285_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "ser285b" = sqrt((v160m_rep2_SER285B$`r_482_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_SER285B$`r_482_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_SER285B$`r_482_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "arg286a" = sqrt((v160m_rep2_ARG286A$`r_286_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ARG286A$`r_286_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ARG286A$`r_286_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "arg286b" = sqrt((v160m_rep2_ARG286B$`r_483_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_ARG286B$`r_483_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_ARG286B$`r_483_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "glu287a" = sqrt((v160m_rep2_GLU287A$`r_287_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_GLU287A$`r_287_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_GLU287A$`r_287_&_C-alpha Z` - v160m_rep2_protein$z)^2),
                         "glu287b" = sqrt((v160m_rep2_GLU287B$`r_484_&_C-alpha X` - v160m_rep2_protein$x)^2 + 
                                            (v160m_rep2_GLU287B$`r_484_&_C-alpha Y` - v160m_rep2_protein$y)^2 + 
                                            (v160m_rep2_GLU287B$`r_484_&_C-alpha Z` - v160m_rep2_protein$z)^2)
)

norm_v160m_rep2 <- data.frame("time" = v160m_rep2$time,
                              "gln278a" = (v160m_rep2$gln278a - v160m_rep2$gln278a[1]),
                              "gln278b" = (v160m_rep2$gln278b - v160m_rep2$gln278b[1]),
                              "asp279a" = (v160m_rep2$asp279a - v160m_rep2$asp279a[1]),
                              "asp279b" = (v160m_rep2$asp279b - v160m_rep2$asp279b[1]),
                              "asp280a" = (v160m_rep2$asp280a - v160m_rep2$asp280a[1]),
                              "asp280b" = (v160m_rep2$asp280b - v160m_rep2$asp280b[1]),
                              "cys281a" = (v160m_rep2$cys281a - v160m_rep2$cys281a[1]),
                              "cys281b" = (v160m_rep2$cys281b - v160m_rep2$cys281b[1]),
                              "ala282a" = (v160m_rep2$ala282a - v160m_rep2$ala282a[1]),
                              "ala282b" = (v160m_rep2$ala282b - v160m_rep2$ala282b[1]),
                              "ala283a" = (v160m_rep2$ala283a - v160m_rep2$ala283a[1]),
                              "ala283b" = (v160m_rep2$ala283b - v160m_rep2$ala283b[1]),
                              "phe284a" = (v160m_rep2$phe284a - v160m_rep2$phe284a[1]),
                              "phe284b" = (v160m_rep2$phe284b - v160m_rep2$phe284b[1]),
                              "ser285a" = (v160m_rep2$ser285a - v160m_rep2$ser285a[1]),
                              "ser285b" = (v160m_rep2$ser285b - v160m_rep2$ser285b[1]),
                              "arg286a" = (v160m_rep2$arg286a - v160m_rep2$arg286a[1]),
                              "arg286b" = (v160m_rep2$arg286b - v160m_rep2$arg286b[1]),
                              "glu287a" = (v160m_rep2$glu287a - v160m_rep2$glu287a[1]),
                              "glu287b" = (v160m_rep2$glu287b - v160m_rep2$glu287b[1])
)

rm(list=ls(pattern="v160m_rep2_"))


#### graficas ####

#open1
ab1_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open1, aes(x = time, y = chainA_norm), color = "palegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = cys281a), color = "palegreen") +
  geom_line(data = open1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "palegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
ab1_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open1, aes(x = time, y = chainA_norm), color = "palegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = ala282a), color = "palegreen") +
  geom_line(data = open1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "palegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-abierta-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ab1_cys281, ab1_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-abierta-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open1, aes(x = time, y = chainA_norm), color = "palegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = arg286a), color = "palegreen") +
  geom_line(data = open1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep1, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "palegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286")
dev.off()



#open2
ab2_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open2, aes(x = time, y = chainA_norm), color = "limegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = cys281a), color = "limegreen") +
  geom_line(data = open2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "limegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
ab2_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open2, aes(x = time, y = chainA_norm), color = "limegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = ala282a), color = "limegreen") +
  geom_line(data = open2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "limegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-abierta-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(ab2_cys281, ab2_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-abierta-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = open2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = open2, aes(x = time, y = chainA_norm), color = "limegreen", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = arg286a), color = "limegreen") +
  geom_line(data = open2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_abierta_rep2, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "limegreen", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286") 
dev.off()

#close1
cl1_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close1, aes(x = time, y = chainA_norm), color = "blue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = cys281a), color = "blue") +
  geom_line(data = close1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "blue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
cl1_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close1, aes(x = time, y = chainA_norm), color = "blue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = ala282a), color = "blue") +
  geom_line(data = close1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "blue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-cerrada-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(cl1_cys281, cl1_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-cerrada-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close1, aes(x = time, y = chainA_norm), color = "blue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = arg286a), color = "blue") +
  geom_line(data = close1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep1, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "blue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286") 
dev.off()

#close2
cl2_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close2, aes(x = time, y = chainA_norm), color = "dodgerblue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = cys281a), color = "dodgerblue") +
  geom_line(data = close2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "dodgerblue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
cl2_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close2, aes(x = time, y = chainA_norm), color = "dodgerblue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = ala282a), color = "dodgerblue") +
  geom_line(data = close2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "dodgerblue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-cerrada-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(cl2_cys281, cl2_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-cerrada-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = close2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = close2, aes(x = time, y = chainA_norm), color = "dodgerblue", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = arg286a), color = "dodgerblue") +
  geom_line(data = close2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_cerrada_rep2, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "dodgerblue", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286") 
dev.off()
#v160m1
vm1_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m1, aes(x = time, y = chainA_norm), color = "deeppink", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = cys281a), color = "deeppink") +
  geom_line(data = v160m1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "deeppink", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
vm1_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m1, aes(x = time, y = chainA_norm), color = "deeppink", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = ala282a), color = "deeppink") +
  geom_line(data = v160m1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "deeppink", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-v160m-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(vm1_cys281, vm1_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-v160m-rep1.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m1, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m1, aes(x = time, y = chainA_norm), color = "deeppink", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = arg286a), color = "deeppink") +
  geom_line(data = v160m1, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep1, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "deeppink", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286") 
dev.off()
#v160m2
vm2_cys281 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m2, aes(x = time, y = chainA_norm), color = "hotpink1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = cys281a), color = "hotpink1") +
  geom_line(data = v160m2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = cys281b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "hotpink1", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("CYSTEINE 281") 
vm2_ala282 <- ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m2, aes(x = time, y = chainA_norm), color = "hotpink1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = ala282a), color = "hotpink1") +
  geom_line(data = v160m2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = ala282b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "hotpink1", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ALANINE 282") 

png("pics-19Sep/dist-max-bfacs-v160m-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggarrange(vm2_cys281, vm2_ala282, nrow = 2, ncol = 1)
dev.off()

png("pics-19Sep/dist-arg286-v160m-rep2.png", height = 300, width = 750, units = "mm", res = 300)
ggplot() + ylim(-0.5, 1.2) + 
  geom_line(data = v160m2, aes(x = time, y = 0), color = "black", alpha = 0.5, linetype = 1) +
  geom_line(data = v160m2, aes(x = time, y = chainA_norm), color = "hotpink1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = arg286a), color = "hotpink1") +
  geom_line(data = v160m2, aes(x = time, y = chainB_norm), color = "darkorange1", alpha = 0.5) + 
  geom_line(data = norm_v160m_rep2, aes(x = time, y = arg286b), color = "darkorange1") +
  annotate("text", x = 0, y = 1, label = "- Chain A", color = "hotpink1", size = 8) +
  annotate("text", x = 0, y = 0.8, label = "- Chain B", color = "darkorange1", size = 8) +
  a + xlab("time") + ylab("distance (nm)") + ggtitle("ARGININE 286")
dev.off()
    
