
# Dane
dane <- read.csv(".c/Desktop/dane.csv", header = T, sep = ";")
str(dane)
dane <- data.frame(dane)
VDR_C <- dane$VDR_FokI[dane$Group == 'cancer']
VDR_K <- dane$VDR_FokI[dane$Group == 'control']
BSM_C <- dane$BSM[dane$Group == 'cancer']
BSM2 <- dane[dane$BSM != '',]
BSM_K <- BSM2$BSM[BSM2$Group == 'control']

H0 <- 'Subpopulacja znajduje się w równowadze Hardyego- Weinberga'
H1 <- 'Subpopulacja nie znajduje się w równowadze Hardyego- Weinberga'  
  
# RÓWNOWAGA HARDY'EGO WEINBERGA
# 1) VDR - cancer

allel1 <- table(VDR_C)
allel1
N1 <- sum(allel1)

p1 <- as.numeric(2*allel1[1] + allel1[2]) / (2*N1)
q1 <- as.numeric(2*allel1[3] + allel1[2]) / (2*N1)
AA1 <- p1^2*N1
AG1 <- 2*p1*q1*N1
GG1 <- q1^2*N1
Chi1 <- as.numeric(((allel1[1]-AA1)^2/AA1)+((allel1[2]-AG1)^2/AG1)+((allel1[3]-GG1)^2/GG1))

if (Chi1 > qchisq(0.95,2)) {
  cat(H1)
} else {
  cat(H0)
}
p_value1 <- 1 - pchisq(Chi1, 2)
cat("Wartość p dla testu VDR - cancer:", p_value1, "\n")
# 2) VDR - control

allel2 <- table(VDR_K)
allel2
N2 <- sum(allel2)

p2 <- as.numeric(2*allel2[1] + allel2[2]) / (2*N2)
q2 <- as.numeric(2*allel2[3] + allel2[2]) / (2*N2)
AA2 <- p2^2*N2
AG2 <- 2*p2*q2*N2
GG2 <- q2^2*N2
Chi2 <- as.numeric(((allel2[1]-AA2)^2/AA2)+((allel2[2]-AG2)^2/AG2)+((allel2[3]-GG2)^2/GG2))

if (Chi2 > qchisq(0.95,2)) {
  cat(H1)
} else {
  cat(H0)
}
p_value2 <- 1 - pchisq(Chi2, 2)
cat("Wartość p dla testu VDR - control:", p_value2, "\n")
# 3) BSM - cancer

allel3 <- table(BSM_C)
allel3
N3 <- sum(allel3)

p3 <- as.numeric((2*allel3[1] + allel3[2] + allel3[3]) / (2*N3))
q3 <- as.numeric((2*allel3[4] + allel3[2] + allel3[3])/ (2*N3))
CC1 <- p3^2*N3
CT1 <- 2*p3*q3*N3
TT1 <- q3^2*N3
Chi3 <- as.numeric(((allel3[1]-CC1)^2/CC1)+((allel3[2] + allel3[3]-CT1)^2/CT1)+((allel3[4]-TT1)^2/TT1))
Chi3
if (Chi3 > qchisq(0.95,3)) {
  cat(H1)
} else {
  cat(H0)
}
p_value3 <- 1 - pchisq(Chi3, 2)
cat("Wartość p dla testu BSM - cancer:", p_value3, "\n")
# 4) BSM - control

allel4 <- table(BSM_K)
allel4
N4 <- sum(allel4)

p4 <- as.numeric((2*allel4[1] + allel4[2] + allel4[3]) / (2*N4))
q4 <- as.numeric((2*allel4[4] + allel4[2] + allel4[3])/ (2*N4))

CC2 <- p4^2*N4
CT2 <- 2*p4*q4*N4
TT2 <- q4^2*N4

Chi4 <- as.numeric(((allel4[1]-CC2)^2/CC2)+((allel4[2]+ allel4[3]-CT2)^2/CT2)+((allel4[4]-TT2)^2/TT2))
Chi4

if (Chi4 > qchisq(0.95,3)) {
  cat(H1)
} else {
  cat(H0)
}
p_value4 <- 1 - pchisq(Chi4, 2)
cat("Wartość p dla testu BSM- control:", p_value4, "\n")
# ODDS RATIO
AA_Fok <- (allel1[1]* (allel2[2]+ allel2[3]))/ ((allel1[2]+allel1[3])*allel2[1])
AG_Fok <- (allel1[2]* (allel2[1]+ allel2[3]))/ ((allel1[1]+allel1[3])*allel2[2])
GG_Fok <- (allel1[3]* (allel2[1]+ allel2[2]))/ ((allel1[2]+allel1[1])*allel2[3])

CC_Bsm <- (allel3[1]* (allel4[2]+ allel4[3]+allel4[4]))/ ((allel3[2]+allel3[3]+allel3[4])*allel4[1])
CT_Bsm <- ((allel3[2] + allel3[3])* (allel4[1]+allel4[4]))/ ((allel3[1]+allel3[4])*(allel4[3] + allel4[2]))
TT_Bsm <-(allel3[4]* (allel4[1]+ allel4[2]+allel4[3]))/ ((allel3[1]+allel3[2]+allel3[3])*allel4[4])




