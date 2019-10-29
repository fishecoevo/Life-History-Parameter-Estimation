# Life History Parameter Estimation Example Code
# P. Kennedy March 2019

# Parameters to be estimated are:
# 1. Asymptotic length
# 2. Gallucci-Quinn Index (omega)
# 3. Average size of largest 5% of catch
# 4. Age at maturation
# 5. Size at maturation
# 6. Instantaneous total mortality
# 7. Annual survivorship

# Code uses catch data for a fish population from a standardized netting program

# Loading in the data
rm(list=ls()) # Remove any lists from the environment
setwd("C:\\Users\\Patrick James\\Desktop") # Set working directory
# Fish population data
Fish= read.table("Life History Example Data.csv", header = TRUE, sep=",", na.strings=".") 
str(Fish) # Checking the structure of variables
# View(Fish) # View data 
Fish$Maturity= as.factor(Fish$Maturity) # Restructuring maturity status

Fish_F= subset(Fish, Fish$Sex== 2) # Subsetting for females only
############################################################################################
# Asymptotic length and Gallucci-Quinn Index:

# Asymptotic length will be estimated using the standard von Bertalanffy Growth Equation,
# as well as the largest 5% of the catch
# Early growth will be estimated using the Gallucci/ Quinn index (omega)

##### 
# Standard VBGM fit
fem.st=c(L = 800, k = 0.1) # Starting values for parameter estimate
Fem_L= nls(Fork_length~L*(1-exp(-k*(ageplus1-0))), data = Fish_F, # Intercept forced through 0
           start = fem.st, trace = T) # Equation for VBGM with fork lengths
coef(Fem_L) # Coefficients for model 
# L inf= 939.08 mm fork length
F_sum=summary(Fem_L)
F_sum$coef[3] # Standard error for L
F_sum$coef[4] # Standard error for K
F_sum$coef[7] # p-value for L
F_sum$coef[8] # p-value for K
summary(Fem_L)[3] # Sigma value for goodness of fit
F_omega= F_sum$coef[1]*F_sum$coef[2] # Gallucci-Quinn Index for early growth near age 0
# Calculating the largest 5% of catch
F_large5= mean(head(sort(Fish_F$Fork_length, decreasing=TRUE), 0.05*dim(Fish_F)[1]))
F_large5 # Mean largest 5% of females is 898.25 mm fork length

# Standard VBGM Plot
F_mean= tapply(Fish_F$Fork_length, Fish_F$ageplus1, mean, na.rm=TRUE)
F_mean= as.data.frame(F_mean)
F_num= as.numeric(rownames(F_mean)) # Changing means to num structure
For_F= data.frame(age=F_num, F_mean) 
ages= c(0:12) # Setting what age the plot x-axis should go to
plot(For_F[,1], For_F[,2], pch=16, col="white", ylab = "Fork Length (mm)", 
     xlab = "Age (years)", ylim = c(0, 1200), xlim = c(0, 12), main = "" )
points(Fish_F$Fork_length ~ Fish_F$ageplus1, col="black") # Plotting actual data points
lines(ages, predict(Fem_L, list (ageplus1 = ages)), lwd =2, col ="black") # Curve

############################################################################################
# Age and length at 50% maturity:

#####
# Subsetting data for known maturation status
Fish_Fem= subset(Fish_F, Fish_F$Maturity != "9")
str(Fish_Fem)

# Running a binomial glm for maturity status ~ age 
age_mat= glm(Fish_Fem$Maturity ~ Fish_Fem$ageplus1, family = binomial)
summary(age_mat)

p_age_mat= summary(age_mat)$coeff[8] # p-value for the age at maturity fit
p_age_mat # good fit

Age_at_mat= (((log(0.5/(1-0.5)))-age_mat$coeff[[1]])/age_mat$coeff[[2]])
Age_at_mat # Age at maturity is 2.43 years old

# Running a binomial glm for maturity status ~ fork length
Length_mat= glm(Fish_Fem$Maturity ~ Fish_Fem$Fork_length, family = binomial)
summary(Length_mat)

p_L_mat<-summary(Length_mat)$coeff[8] # p-value for the length at maturity fit
p_L_mat # good fit

L_at_mat= (((log(0.5/(1-0.5)))-Length_mat$coeff[[1]])/Length_mat$coeff[[2]])
L_at_mat # Size at maturity is 419.69 mm fork length

############################################################################################
# Robson-Chapman estimation of instantaneous total mortality: 

# Function for mode of a variable
# This will be needed for age at recruitment to gear
Mode= function(x) {
  ux= unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#####
# Age counts
age_count= table(Fish_F$ageplus1)
age_count

# Make age count a data frame
age_count= as.data.frame(age_count)
names(age_count)= c("Age","Cnt")
age_classes= nrow(age_count)
age_classes

# There may be missing ages
# Here we create a seqential vector from youngest to oldest age

# Minimum age
min_age= 0

# Max age
max_age= max(Fish_F$ageplus1)

# Age at recruitment aka modal age class
age_recruit= Mode(Fish_F$ageplus1)

all_ages= data.frame("Age"= seq(min_age, max_age))
age_classes= nrow(all_ages)

if(nrow(all_ages)!=nrow(age_count)){
  # If the number of age classes isn't the same, we need to add in
  # some zeroes in the right places
  age_count= merge(all_ages, age_count, by='Age', all.x=TRUE)
  age_count$Cnt[is.na(age_count$Cnt)] <- 0
}      

age_count$Age= as.numeric(as.character(age_count$Age))

# Trim age count so that it has only ages past the mode
age_trim= subset (age_count, age_count$Age >= age_recruit)

# Robson Chapman estimator. First Calculate coded age and coded age * count
age_trim$coded_age= age_trim$Age-age_recruit
age_trim$coded_count= age_trim$Cnt*age_trim$coded_age

# Robson Chapman parameters; S= annual survival, Z= inst. total mortality
T= sum(age_trim$coded_count)
n= sum(age_trim$Cnt)
S= T/(n+T-1) # Survivorship= 0.75
Z= round(-log(S),3) # Inst. total mortality= 0.28
A= 1-exp(-Z)

############################################################################################
# That's it for now!
# These parameters form a great basis for assessing the status of a fish population 
# This code can easily be thrown into loops to calculate values for multiple populations
# Enjoy!


