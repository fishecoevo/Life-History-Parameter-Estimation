# Exploring Electrofishing Data from 1974-2018

# Loading in the data
rm(list=ls()) # Remove any lists from the environment
setwd("H:\\AWQM Biological Monitoring\\Fish Data\\PJK DB Project\\R") # Set working directory
# AWQM water chemistry yearly means
AWQM_WQ= read.table("AWQM Yearly Means 1971-2018.csv", header = TRUE, sep=",", na.strings=".") 
str(AWQM_WQ)
# AWQM fish data
AWQM_EF= read.table("MWRD AWQM Fish Data 2001-2018.csv", header = TRUE, sep=",", na.strings=".", fill= TRUE) 
str(AWQM_EF) # Checking the structure of variables
# View(AWQM_EF) # View data

# Long-term fish data
AC_EF= read.table("MWRD Fish Data 1974-2000.csv", header = TRUE, sep=",", na.strings=".", fill= TRUE) 
str(AC_EF) # Checking the structure of variables
# View(AC_EF) # View data

# Simplifying data for combining DC and AC gear
AWQM_EF2= AWQM_EF[c("Fish.Species", "TL.mm","Wt.g","Year")]

# Matching AC column names for combination
AC_EF$TL.mm= AC_EF$TL_mm
AC_EF$Wt.g= AC_EF$Wt_g
AC_EF2= AC_EF[c("Fish.Species", "TL.mm","Wt.g","Year")]

# New combined data frame
EF_Combined= rbind(AWQM_EF2, AC_EF2)


# Subsetting for species of interest
unique(EF_Combined$Fish.Species)
Combined_Bluegill= subset(EF_Combined, EF_Combined$Fish.Species== "Bluegill") 
Combined_Carp= subset(EF_Combined, EF_Combined$Fish.Species== "Carp") 
Combined_GreenSF= subset(EF_Combined, EF_Combined$Fish.Species== "Green sunfish") 
Combined_Shad= subset(EF_Combined, EF_Combined$Fish.Species== "Gizzard shad") 
Combined_Largemouth= subset(EF_Combined, EF_Combined$Fish.Species== "Largemouth bass") 
Combined_Pumpkinseed= subset(EF_Combined, EF_Combined$Fish.Species== "Pumpkinseed") 
Combined_RockBass= subset(EF_Combined, EF_Combined$Fish.Species== "Rock bass") 
Combined_Smallmouth= subset(EF_Combined, EF_Combined$Fish.Species== "Smallmouth bass") 
Combined_BM= subset(EF_Combined, EF_Combined$Fish.Species== "Bluntnose minnow") 

#####################################################################################
# Plotting weight~length relationships 
# and fitting linear/ nonlinear models 

#####
# Bluegill
str(Combined_Bluegill)
Combined_Bluegill= na.omit(Combined_Bluegill)
Combined_Bluegill= subset(Combined_Bluegill, Combined_Bluegill$TL.mm > 0)

plot(Wt.g~TL.mm, data= Combined_Bluegill)
plot(log(Wt.g)~log(TL.mm), data= Combined_Bluegill)
Bluegill_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Bluegill)
summary(Bluegill_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Bluegill, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Bluegill), col= "red", lwd= 2)
New_Bluegill= data.frame(TL.mm= seq(min(Combined_Bluegill$TL.mm), 
                                    max(Combined_Bluegill$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Bluegill)
lines(New_Bluegill$TL.mm, exp(predict(Bluegill_LW, New_Bluegill)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(70, 250, substitute(b0*x^b1, list(b0=exp(coef(Bluegill_LW)[1]), b1=coef(Bluegill_LW)[2])))
text(70, 230, substitute(plain("R-square: ") * r2, list(r2=summary(Bluegill_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Bluegill$Wt.est= 7.552729E-6*(Combined_Bluegill$TL.mm^3.199134)
Combined_Bluegill$Kn= Combined_Bluegill$Wt.g/Combined_Bluegill$Wt.est
#Combined_Bluegill= subset(Combined_Bluegill, Combined_Bluegill$Kn < 10)
plot(Combined_Bluegill$Kn~Combined_Bluegill$Year)
plot(Combined_Bluegill$Kn~Combined_Bluegill$TL.mm)

#####
# Carp
str(Combined_Carp)
Combined_Carp= na.omit(Combined_Carp)
Combined_Carp= subset(Combined_Carp, Combined_Carp$TL.mm > 0)
Combined_Carp= subset(Combined_Carp, Combined_Carp$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_Carp)
plot(log(Wt.g)~log(TL.mm), data= Combined_Carp)
Carp_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Carp)
summary(Carp_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Carp, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Carp), col= "red", lwd= 2)
New_Carp= data.frame(TL.mm= seq(min(Combined_Carp$TL.mm), 
                                max(Combined_Carp$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Carp)
lines(New_Carp$TL.mm, exp(predict(Carp_LW, New_Carp)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(250, 11200, substitute(b0*x^b1, list(b0=exp(coef(Carp_LW)[1]), b1=coef(Carp_LW)[2])))
text(250, 10000, substitute(plain("R-square: ") * r2, list(r2=summary(Carp_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Carp$Wt.est= 2.478727E-5*(Combined_Carp$TL.mm^2.930382)
Combined_Carp$Kn= Combined_Carp$Wt.g/Combined_Carp$Wt.est
#Combined_Carp= subset(Combined_Carp, Combined_Carp$Kn < 7)
plot(Combined_Carp$Kn~Combined_Carp$Year)
plot(Combined_Carp$Kn~Combined_Carp$TL.mm)

#####
# Green sunfish
str(Combined_GreenSF)
Combined_GreenSF= na.omit(Combined_GreenSF)
Combined_GreenSF= subset(Combined_GreenSF, Combined_GreenSF$TL.mm > 0)
Combined_GreenSF= subset(Combined_GreenSF, Combined_GreenSF$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_GreenSF)
plot(log(Wt.g)~log(TL.mm), data= Combined_GreenSF)
GreenSF_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_GreenSF)
summary(GreenSF_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_GreenSF, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_GreenSF), col= "red", lwd= 2)
New_GreenSF= data.frame(TL.mm= seq(min(Combined_GreenSF$TL.mm), 
                                   max(Combined_GreenSF$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_GreenSF)
lines(New_GreenSF$TL.mm, exp(predict(GreenSF_LW, New_GreenSF)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(90, 250, substitute(b0*x^b1, list(b0=exp(coef(GreenSF_LW)[1]), b1=coef(GreenSF_LW)[2])))
text(90, 220, substitute(plain("R-square: ") * r2, list(r2=summary(GreenSF_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_GreenSF$Wt.est= 1.445341E-5*(Combined_GreenSF$TL.mm^3.074909)
Combined_GreenSF$Kn= Combined_GreenSF$Wt.g/Combined_GreenSF$Wt.est
#Combined_GreenSF= subset(Combined_GreenSF, Combined_GreenSF$Kn < 10)
plot(Combined_GreenSF$Kn~Combined_GreenSF$Year)
plot(Combined_GreenSF$Kn~Combined_GreenSF$TL.mm)

#####
# Gizzard shad
str(Combined_Shad)
Combined_Shad= na.omit(Combined_Shad)
Combined_Shad= subset(Combined_Shad, Combined_Shad$TL.mm > 0)
Combined_Shad= subset(Combined_Shad, Combined_Shad$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_Shad)
plot(log(Wt.g)~log(TL.mm), data= Combined_Shad)
Shad_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Shad)
summary(Shad_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Shad, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Shad), col= "red", lwd= 2)
New_Shad= data.frame(TL.mm= seq(min(Combined_Shad$TL.mm), 
                                max(Combined_Shad$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Shad)
lines(New_Shad$TL.mm, exp(predict(Shad_LW, New_Shad)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(130, 1800, substitute(b0*x^b1, list(b0=exp(coef(Shad_LW)[1]), b1=coef(Shad_LW)[2])))
text(130, 1650, substitute(plain("R-square: ") * r2, list(r2=summary(Shad_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Shad$Wt.est= 9.101293E-6*(Combined_Shad$TL.mm^3.031049)
Combined_Shad$Kn= Combined_Shad$Wt.g/Combined_Shad$Wt.est
#Combined_Shad= subset(Combined_Shad, Combined_Shad$Kn < 10)
plot(Combined_Shad$Kn~Combined_Shad$Year)
plot(Combined_Shad$Kn~Combined_Shad$TL.mm)

#####
# Largemouth bass
str(Combined_Largemouth)
Combined_Largemouth= na.omit(Combined_Largemouth)
Combined_Largemouth= subset(Combined_Largemouth, Combined_Largemouth$TL.mm > 0)
Combined_Largemouth= subset(Combined_Largemouth, Combined_Largemouth$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_Largemouth)
plot(log(Wt.g)~log(TL.mm), data= Combined_Largemouth)
Largemouth_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Largemouth)
summary(Largemouth_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Largemouth, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Largemouth), col= "red", lwd= 2)
New_Largemouth= data.frame(TL.mm= seq(min(Combined_Largemouth$TL.mm), 
                                      max(Combined_Largemouth$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Largemouth)
lines(New_Largemouth$TL.mm, exp(predict(Largemouth_LW, New_Largemouth)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(120, 1700, substitute(b0*x^b1, list(b0=exp(coef(Largemouth_LW)[1]), b1=coef(Largemouth_LW)[2])))
text(120, 1550, substitute(plain("R-square: ") * r2, list(r2=summary(Largemouth_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Largemouth$Wt.est= 8.779749E-6*(Combined_Largemouth$TL.mm^3.085543)
Combined_Largemouth$Kn= Combined_Largemouth$Wt.g/Combined_Largemouth$Wt.est
#Combined_Largemouth= subset(Combined_Largemouth, Combined_Largemouth$Kn < 8)
plot(Combined_Largemouth$Kn~Combined_Largemouth$Year)
plot(Combined_Largemouth$Kn~Combined_Largemouth$TL.mm)

#####
# Pumpkinseed
str(Combined_Pumpkinseed)
Combined_Pumpkinseed= na.omit(Combined_Pumpkinseed)
Combined_Pumpkinseed= subset(Combined_Pumpkinseed, Combined_Pumpkinseed$TL.mm > 0)
Combined_Pumpkinseed= subset(Combined_Pumpkinseed, Combined_Pumpkinseed$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_Pumpkinseed)
plot(log(Wt.g)~log(TL.mm), data= Combined_Pumpkinseed)
Pumpkinseed_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Pumpkinseed)
summary(Pumpkinseed_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Pumpkinseed, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Pumpkinseed), col= "red", lwd= 2)
New_Pumpkinseed= data.frame(TL.mm= seq(min(Combined_Pumpkinseed$TL.mm), 
                                       max(Combined_Pumpkinseed$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Pumpkinseed)
lines(New_Pumpkinseed$TL.mm, exp(predict(Pumpkinseed_LW, New_Pumpkinseed)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(75, 500, substitute(b0*x^b1, list(b0=exp(coef(Pumpkinseed_LW)[1]), b1=coef(Pumpkinseed_LW)[2])))
text(75, 450, substitute(plain("R-square: ") * r2, list(r2=summary(Pumpkinseed_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Pumpkinseed$Wt.est= 9.839659E-6*(Combined_Pumpkinseed$TL.mm^3.166119)
Combined_Pumpkinseed$Kn= Combined_Pumpkinseed$Wt.g/Combined_Pumpkinseed$Wt.est
#Combined_Pumpkinseed= subset(Combined_Pumpkinseed, Combined_Pumpkinseed$Kn < 8)
plot(Combined_Pumpkinseed$Kn~Combined_Pumpkinseed$Year)
plot(Combined_Pumpkinseed$Kn~Combined_Pumpkinseed$TL.mm)

#####
# Rock bass
str(Combined_RockBass)
Combined_RockBass= na.omit(Combined_RockBass)
Combined_RockBass= subset(Combined_RockBass, Combined_RockBass$TL.mm > 0)
Combined_RockBass= subset(Combined_RockBass, Combined_RockBass$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_RockBass)
plot(log(Wt.g)~log(TL.mm), data= Combined_RockBass)
RockBass_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_RockBass)
summary(RockBass_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_RockBass, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_RockBass), col= "red", lwd= 2)
New_RockBass= data.frame(TL.mm= seq(min(Combined_RockBass$TL.mm), 
                                        max(Combined_RockBass$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_RockBass)
lines(New_RockBass$TL.mm, exp(predict(RockBass_LW, New_RockBass)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(90, 250, substitute(b0*x^b1, list(b0=exp(coef(RockBass_LW)[1]), b1=coef(RockBass_LW)[2])))
text(90, 220, substitute(plain("R-square: ") * r2, list(r2=summary(RockBass_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_RockBass$Wt.est= 1.812288E-5*(Combined_RockBass$TL.mm^3.014802)
Combined_RockBass$Kn= Combined_RockBass$Wt.g/Combined_RockBass$Wt.est
#Combined_RockBass= subset(Combined_RockBass, Combined_RockBass$Kn < 8)
plot(Combined_RockBass$Kn~Combined_RockBass$Year)
RB_KN_lm= lm(Combined_RockBass$Kn~Combined_RockBass$Year)
#summary(RB_KN_lm)
abline(lm(Combined_RockBass$Kn~Combined_RockBass$Year))
plot(Combined_RockBass$Kn~Combined_RockBass$TL.mm)

#####
# Smallmouth bass
str(Combined_Smallmouth)
Combined_Smallmouth= na.omit(Combined_Smallmouth)
Combined_Smallmouth= subset(Combined_Smallmouth, Combined_Smallmouth$TL.mm > 0)
Combined_Smallmouth= subset(Combined_Smallmouth, Combined_Smallmouth$Wt.g > 0)

plot(Wt.g~TL.mm, data= Combined_Smallmouth)
plot(log(Wt.g)~log(TL.mm), data= Combined_Smallmouth)
Smallmouth_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_Smallmouth)
summary(Smallmouth_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_Smallmouth, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_Smallmouth), col= "red", lwd= 2)
New_Smallmouth= data.frame(TL.mm= seq(min(Combined_Smallmouth$TL.mm), 
                                      max(Combined_Smallmouth$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_Smallmouth)
lines(New_Smallmouth$TL.mm, exp(predict(Smallmouth_LW, New_Smallmouth)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(140, 1000, substitute(b0*x^b1, list(b0=exp(coef(Smallmouth_LW)[1]), b1=coef(Smallmouth_LW)[2])))
text(140, 880, substitute(plain("R-square: ") * r2, list(r2=summary(Smallmouth_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_Smallmouth$Wt.est= 1.491204E-5*(Combined_Smallmouth$TL.mm^2.973249)
Combined_Smallmouth$Kn= Combined_Smallmouth$Wt.g/Combined_Smallmouth$Wt.est
#Combined_Smallmouth= subset(Combined_Smallmouth, Combined_Smallmouth$Kn < 5)
plot(Combined_Smallmouth$Kn~Combined_Smallmouth$Year)
plot(Combined_Smallmouth$Kn~Combined_Smallmouth$TL.mm)

#####
# Bluntnose minnow
str(Combined_BM)
Combined_BM= na.omit(Combined_BM)
Combined_BM= subset(Combined_BM, Combined_BM$TL.mm > 0)

plot(Wt.g~TL.mm, data= Combined_BM)
plot(log(Wt.g)~log(TL.mm), data= Combined_BM)
Bluegill_LW= lm(log(Wt.g)~log(TL.mm), data= Combined_BM)
summary(Bluegill_LW)
plot(log(Wt.g)~log(TL.mm), data= Combined_BM, pch= 16)
abline(lm(log(Wt.g)~log(TL.mm), data= Combined_BM), col= "red", lwd= 2)
New_Bluegill= data.frame(TL.mm= seq(min(Combined_BM$TL.mm), 
                                    max(Combined_BM$TL.mm), len=100))
plot(Wt.g~TL.mm, data= Combined_BM)
lines(New_Bluegill$TL.mm, exp(predict(Bluegill_LW, New_Bluegill)),
      col= "red", lwd= 2)
# Plotting equation and r-square value
text(50, 35, substitute(b0*x^b1, list(b0=exp(coef(Bluegill_LW)[1]), b1=coef(Bluegill_LW)[2])))
text(50, 32, substitute(plain("R-square: ") * r2, list(r2=summary(Bluegill_LW)$r.squared)))

# Calculating LeCren's relative condition factor based on this relationship
Combined_BM$Wt.est= 3.913539E-6*(Combined_BM$TL.mm^3.226072)
Combined_BM$Kn= Combined_BM$Wt.g/Combined_BM$Wt.est
#Combined_BM= subset(Combined_BM, Combined_BM$Kn < 10)
plot(Combined_BM$Kn~Combined_BM$Year)
plot(Combined_BM$Kn~Combined_BM$TL.mm)

#####################################################################################

# Subsetting for Calumet EF sites
Calumet_EF= subset(AWQM_EF, AWQM_EF$Sta.No== "WW55" | 
                      AWQM_EF$Sta.No== "WW56" |
                      AWQM_EF$Sta.No== "WW57" |    
                      AWQM_EF$Sta.No== "WW58" |
                      AWQM_EF$Sta.No== "WW59" |
                      AWQM_EF$Sta.No== "WW76")
str(Calumet_EF)
table(Calumet_EF$Fish.Species, Calumet_EF$Year)
unique(Calumet_EF$Fish.Species)

table(AWQM_WQ$Location.ID) # Missing site #57
# Subsetting for Calumet WQ sites
Calumet_WQ= subset(AWQM_WQ, AWQM_WQ$Location.ID== "55" |  
                      AWQM_WQ$Location.ID== "56" |
                      AWQM_WQ$Location.ID== "58" |
                      AWQM_WQ$Location.ID== "59" |
                      AWQM_WQ$Location.ID== "76")
str(Calumet_WQ)

# Subsetting Calumet sites for species of interest
Calumet_Bluegill= subset(Calumet_EF, Calumet_EF$Fish.Species== "Bluegill") 
Calumet_Carp= subset(Calumet_EF, Calumet_EF$Fish.Species== "Carp") 
Calumet_GreenSF= subset(Calumet_EF, Calumet_EF$Fish.Species== "Green sunfish") 
Calumet_Shad= subset(Calumet_EF, Calumet_EF$Fish.Species== "Gizzard shad") 
Calumet_Largemouth= subset(Calumet_EF, Calumet_EF$Fish.Species== "Largemouth bass") 
Calumet_Pumpkinseed= subset(Calumet_EF, Calumet_EF$Fish.Species== "Pumpkinseed") 
Calumet_RockBass= subset(Calumet_EF, Calumet_EF$Fish.Species== "Rock bass") 
Calumet_Smallmouth= subset(Calumet_EF, Calumet_EF$Fish.Species== "Smallmouth bass")


