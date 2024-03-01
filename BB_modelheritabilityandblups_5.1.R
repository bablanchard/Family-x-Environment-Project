library(readxl)
library("PerformanceAnalytics")
library(car)
library(dplyr)
library(nlme)




#=============QC of the data===========
#This block will read in the master data, take out New Iberia year 2 data, calculate stalk weight, 
#replace missing stool values with survival, use individual data to estimate family data,
#remove certain families with missing data, make a csv dataset that is clean,
#and calculate phenotypic correlations in a chart

#read master
df <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/Master.csv", header = TRUE)
#take out NI Y2
df_filtered <- df[!(df$PlantYear == 2021 & df$Location == "NI"), ]
#make Bustalks numeric with the number of stalks in the bundle
df_filtered$Bustalks <- as.numeric(gsub("^\\*([0-9]+)$", "\\1", df_filtered$Bustalks))
df_filtered$Bustalks[is.na(df_filtered$Bustalks) & is.na(df_filtered$Bustalks)] <- 10
#calcualte SW
df_filtered$SW <- df_filtered$Bu.Wt / df_filtered$Bustalks
#if Stool = missing then fill in with Survival
df_filtered$Stool[is.na(df_filtered$Stool) | df_filtered$Stool == ""] <- df_filtered$Survival[is.na(df_filtered$Stool) | df_filtered$Stool == ""]
#read individual data from New Roads
individ <- read_excel("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/SeedlingIndividual.xlsx", sheet = "NRdatasheet")
#remove observations in individ that contain NAs for multiple columns
individ <- individ[complete.cases(individ[, c("Stalks", "Height", "Dia1", "Dia2", "Brix")]), ]
#add location column and PlantYear according to crop
individ$Location <- "NR"
individ$PlantYear <- ifelse(individ$Crop == 2, 2020, ifelse(individ$Crop == 1, 2021, NA))
#average Dia1 and Dia2
individ$Diam <- rowMeans(individ[, c("Dia1", "Dia2")], na.rm = TRUE)
#calculate average family diameter from individuals and apply to each observation of family to input into master
famdiam <- aggregate(Diam ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famdiam, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Diam[is.na(individ$Diam)] <- individ$Diam_avg[is.na(individ$Diam)]
names(individ)[names(individ) == "Diam_avg"] <- "famdiam"
#calculate average family Brix from individuals and apply to each observation of family to input into master
famBrix <- aggregate(Brix ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famBrix, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Brix[is.na(individ$Brix)] <- individ$famBrix[is.na(individ$Brix)]
names(individ)[names(individ) == "Brix_avg"] <- "famBrix"
#calculate average family diameter from individuals and apply to each observation of family to input into master
famHt <- aggregate(Height ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famHt, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Height[is.na(individ$Height)] <- individ$famHt[is.na(individ$Height)]
names(individ)[names(individ) == "Height_avg"] <- "famHt"
#calculate average family stalks from individuals and apply to each observation of family to input into master
famStlk <- aggregate(Stalks ~ PlantYear + Crop + Plot + Rep, data = individ, FUN = mean, na.rm = TRUE)
individ <- merge(individ, famStlk, by = c("PlantYear", "Crop", "Plot", "Rep"), suffixes = c("", "_avg"))
individ$Stalks[is.na(individ$Stalks)] <- individ$famStlk[is.na(individ$Stalks)]
names(individ)[names(individ) == "Stalks_avg"] <- "famStlk"
#remove observation for SG PlantYear2020 Plot 49
df_filtered <- df_filtered[!(df_filtered$PlantYear == 2020 & df_filtered$Location == "SG" & df_filtered$Plot == 49), ]
#remove observation for SG PlantYear2021 Plot 28
df_filtered <- df_filtered[!(df_filtered$PlantYear == 2021 & df_filtered$Location == "SG" & df_filtered$Plot == 28), ]
# Match combinations of PlantYear, Crop, Plot, Rep, and Location in famStlk_avg to df_filtered and add the famStlk column
df_filtered$famStlk <- individ$famStlk[match(paste(df_filtered$PlantYear, df_filtered$Crop, df_filtered$Plot, df_filtered$Rep, df_filtered$Location), 
                                             paste(individ$PlantYear, individ$Crop, individ$Plot, individ$Rep, individ$Location))]
#now multiple the famStalk column to the number of stools to get an estimated stalk count for those families where only individual data was taken (these will be inflated)
df_filtered$Stalk <- ifelse(is.na(df_filtered$famStlk), df_filtered$Stalk, df_filtered$famStlk * df_filtered$Stool)
# Write the dataframe to a CSV file
write.csv(df_filtered, file = "Master_1.2.csv", row.names = TRUE)

# phenotypic correlation between traits
sapply(df_filtered[, 11:22, 27], class)
df_filtered[, 11:22] <- lapply(df_filtered[, 11:22], as.numeric)
df_filtered[, 27] <- as.numeric(df_filtered[, 27])
round(cor(df_filtered[, c(11:22, 27)], use = "pairwise"), 2)
data_subset <- df_filtered[, c(11:18, 27)]
#install.packages("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(data_subset)), histogram = TRUE, pch = 1)

# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
Stalk <- df_filtered[complete.cases(df_filtered$Stalk) & df_filtered$Stalk != 0,]
Height <- df_filtered[complete.cases(df_filtered$Height) & df_filtered$Height != 0,]
Dia <- df_filtered[complete.cases(df_filtered$Dia) & df_filtered$Dia != 0,]
PlotWeight <- df_filtered[complete.cases(df_filtered$PlotWeight) & df_filtered$PlotWeight != 0,]
Brix <- df_filtered[complete.cases(df_filtered$Brix) & df_filtered$Brix != 0,]
Fiber <- df_filtered[complete.cases(df_filtered$Fiber) & df_filtered$Fiber != 0,]
SW <- df_filtered[complete.cases(df_filtered$SW) & df_filtered$SW != 0,]
Sucrose <- df_filtered[complete.cases(df_filtered$Sucrose..) & df_filtered$Sucrose.. != 0,]


hist(TRS$TRS)
hist(Stalk$Stalk)
hist(Height$Height)
hist(Dia$Dia)
hist(PlotWeight$PlotWeight)
hist(Brix$Brix)
hist(Fiber$Fiber)
hist(SW$SW)

boxplot(TRS$TRS)
boxplot(Stalk$Stalk)
boxplot(Height$Height)
boxplot(Dia$Dia)
boxplot(PlotWeight$PlotWeight)
boxplot(Brix$Brix)
boxplot(Fiber$Fiber)
boxplot(SW$SW)



#==============Modeling=================
#Create empty data frame for BLUP output
DataOutput <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutput$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutput$Row <- c(1:50)
DataOutput$Cross <- DataOutput$Entry
DataOutput <- subset(DataOutput, select = -Entry)


#this empty dataframe is for variance components
DataVarComp <- data.frame()
DataVarCompOutput <- data.frame()
HeritabilityData <- data.frame()


#===========TRS==============================
# Load the 'car' package
table(TRS$Cross)
table(TRS$Location)
table(TRS$Rep)
table(TRS$Crop)
table(TRS$Env)

table(TRS$Cross, TRS$Location)
table(TRS$Crop, TRS$Location, TRS$PlantYear)
table(TRS$Cross, TRS$Location, TRS$PlantYear)
table(TRS$Cross, TRS$Location, TRS$PlantYear, TRS$Crop)


#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
TRS <- subset(TRS, !(Cross %in% c("CP18-1611", "HB 17-3208")))

#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(TRS$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
TRSharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", TRSharmonic_mean, "\n")

#Combine columns to create the 'Env' factor
TRS <- TRS %>%
  mutate(Env = interaction(Location, PlantYear))
table(TRS$Env)

TRS$Inter <- interaction(TRS$Cross, 
                         TRS$Env)
TRS$FxCrop <- interaction(TRS$Cross, 
                         TRS$Crop)

# Create the linear mixed model using lme function and specify heterogenous variance 
TRSmod1 <- lme(fixed = TRS ~ 1,  
               random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), 
                             PlantYear = pdDiag(~1), Location = pdDiag(~1), 
                             FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
               weights = varIdent(form = ~1 | Env),
               data = TRS,
               method = 'REML')
summary(TRSmod1)
plot(TRSmod1)
residuals <- residuals(TRSmod1)

qqnorm(residuals)
qqline(residuals)

# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(TRSmod1)
# Create an empty dataframe to store the variance components
DataVarComp <- data.frame(Component = character(),
                          Variance = numeric(),
                          stringsAsFactors = FALSE)

# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "Cross",
                                               Variance = variance_components[2]))
}
if (!is.null(variance_components[4])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "Crop",
                                               Variance = variance_components[4]))
}
if (!is.null(variance_components[6])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "Env",
                                               Variance = variance_components[6]))
}
if (!is.null(variance_components[8])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "FxCrop",
                                               Variance = variance_components[8]))
}
if (!is.null(variance_components[10])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "Inter",
                                               Variance = variance_components[10]))
}
if (!is.null(variance_components[11])) {
  DataVarComp <- rbind(DataVarComp, data.frame(Component = "Residual",
                                               Variance = variance_components[11]))
}

# Print the updated DataVarComp dataframe
DataVarComp$TRS <- DataVarComp$Variance
DataVarComp$Variance <- NULL
print(DataVarComp)

# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(TRSmod1)

# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(TRSmod1)$Cross

# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)

# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")

# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_TRS <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_TRS")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_TRS"] <- "TRS"

# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(TRS = as.numeric(as.character(TRS)))


# Filter and calculate result
TRSH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(TRS, na.rm = TRUE) / sum(DataVarComp$TRS, na.rm = TRUE))




# Assuming TRSH2 contains the calculated value
value_to_assign <- TRSH2$result_value

# Create a new row with Component "Heritability" and TRS value
new_row <- data.frame(Component = "Heritability", TRS = value_to_assign)

# Add the new row to DataVarComp
DataVarComp <- rbind(DataVarComp, new_row)


#============PlotWeight====
#==============SW===========
# Load the 'car' package
table(SW$Cross)
table(SW$Location)
table(SW$Rep)
table(SW$Crop)
table(SW$Cross, SW$Location)
table(SW$Crop, SW$Location, SW$PlantYear)
table(SW$Cross, SW$Location, SW$PlantYear)
table(SW$Cross, SW$Location, SW$PlantYear, SW$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
SW <- subset(SW, !(Cross %in% c("CP18-1611", "HB 17-3208")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(SW$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
SWharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", SWharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
SW <- SW %>%
  mutate(Env = interaction(Location, PlantYear))
table(SW$Env)

SW$Inter <- interaction(SW$Cross, 
                        SW$Env)
SW$FxCrop <- interaction(SW$Cross, 
                         SW$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
SWmod1 <- lme(fixed = SW ~ 1,  
                      random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
                      weights = varIdent(form = ~1 | Env),
                      data = SW,
                      method = 'REML')
summary(SWmod1)
plot(SWmod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(SWmod1)
# Create an empty dataframe to store the variance components
DataVarComp$SW <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$SW[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$SW[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$SW[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$SW[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$SW[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$SW[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(SWmod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(SWmod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_SW <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_SW")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_SW"] <- "SW"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(SW = as.numeric(as.character(SW)))
# Filter and calculate result
SWH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(SW, na.rm = TRUE) / sum(DataVarComp$SW, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- SWH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$SW[DataVarComp$Component == "Heritability"] <- value_to_assign

#=============Brix==========
# Load the 'car' package
table(Brix$Cross)
table(Brix$Location)
table(Brix$Rep)
table(Brix$Crop)
table(Brix$Cross, Brix$Location)
table(Brix$Crop, Brix$Location, Brix$PlantYear)
table(Brix$Cross, Brix$Location, Brix$PlantYear)
table(Brix$Cross, Brix$Location, Brix$PlantYear, Brix$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
Brix <- subset(Brix, !(Cross %in% c("CP18-1611", "HB 17-3208")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(Brix$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
Brixharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", Brixharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
Brix <- Brix %>%
  mutate(Env = interaction(Location, PlantYear))
table(Brix$Env)

Brix$Inter <- interaction(Brix$Cross, 
                          Brix$Env)
Brix$FxCrop <- interaction(Brix$Cross, 
                         Brix$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
Brixmod1 <- lme(fixed = Brix ~ 1,  
              random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
              weights = varIdent(form = ~1 | Env),
              data = Brix,
              method = 'REML')
summary(Brixmod1)
plot(Brixmod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(Brixmod1)
# Create an empty dataframe to store the variance components
DataVarComp$Brix <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$Brix[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$Brix[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$Brix[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$Brix[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$Brix[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$Brix[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(Brixmod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(Brixmod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_Brix <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_Brix")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_Brix"] <- "Brix"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(Brix = as.numeric(as.character(Brix)))
# Filter and calculate result
BrixH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(Brix, na.rm = TRUE) / sum(DataVarComp$Brix, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- BrixH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$Brix[DataVarComp$Component == "Heritability"] <- value_to_assign

#===============Fiber=============
# Load the 'car' package
table(Fiber$Cross)
table(Fiber$Location)
table(Fiber$Rep)
table(Fiber$Crop)
table(Fiber$Cross, Fiber$Location)
table(Fiber$Crop, Fiber$Location, Fiber$PlantYear)
table(Fiber$Cross, Fiber$Location, Fiber$PlantYear)
table(Fiber$Cross, Fiber$Location, Fiber$PlantYear, Fiber$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
Fiber <- subset(Fiber, !(Cross %in% c("CP18-1611", "HB 17-3208")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(Fiber$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
Fiberharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", Fiberharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
Fiber <- Fiber %>%
  mutate(Env = interaction(Location, PlantYear))
table(Fiber$Env)

Fiber$Inter <- interaction(Fiber$Cross, 
                           Fiber$Env)
Fiber$FxCrop <- interaction(Fiber$Cross, 
                            Fiber$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
Fibermod1 <- lme(fixed = Fiber ~ 1,  
                random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
                weights = varIdent(form = ~1 | Env),
                data = Fiber,
                method = 'REML')
summary(Fibermod1)
plot(Fibermod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(Fibermod1)
# Create an empty dataframe to store the variance components
DataVarComp$Fiber <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$Fiber[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$Fiber[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$Fiber[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$Fiber[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$Fiber[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$Fiber[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(Fibermod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(Fibermod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_Fiber <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_Fiber")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_Fiber"] <- "Fiber"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(Fiber = as.numeric(as.character(Fiber)))
# Filter and calculate result
FiberH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(Fiber, na.rm = TRUE) / sum(DataVarComp$Fiber, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- FiberH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$Fiber[DataVarComp$Component == "Heritability"] <- value_to_assign

#==============Dia==========
#Load the 'car' package
table(Dia$Cross)
table(Dia$Location)
table(Dia$Rep)
table(Dia$Crop)
table(Dia$Cross, Dia$Location)
table(Dia$Crop, Dia$Location, Dia$PlantYear)
table(Dia$Cross, Dia$Location, Dia$PlantYear)
table(Dia$Cross, Dia$Location, Dia$PlantYear, Dia$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
Dia <- subset(Dia, !(Cross %in% c("CP18-1611")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(Dia$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
Diaharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", Diaharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
Dia <- Dia %>%
  mutate(Env = interaction(Location, PlantYear))
table(Dia$Env)

Dia$Inter <- interaction(Dia$Cross, 
                         Dia$Env)
Dia$FxCrop <- interaction(Dia$Cross, 
                          Dia$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
Diamod1 <- lme(fixed = Dia ~ 1,  
                 random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
                 weights = varIdent(form = ~1 | Env),
                 data = Dia,
                 method = 'REML')
summary(Diamod1)
plot(Diamod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(Diamod1)
# Create an empty dataframe to store the variance components
DataVarComp$Dia <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$Dia[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$Dia[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$Dia[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$Dia[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$Dia[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$Dia[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(Diamod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(Diamod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_Dia <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_Dia")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_Dia"] <- "Dia"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(Dia = as.numeric(as.character(Dia)))
# Filter and calculate result
DiaH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(Dia, na.rm = TRUE) / sum(DataVarComp$Dia, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- DiaH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$Dia[DataVarComp$Component == "Heritability"] <- value_to_assign
# Assuming Dia is the column you want to round
DataVarComp$Dia <- round(DataVarComp$Dia, digits = 7)

#==============Height=========
#Load the 'car' package
table(Height$Cross)
table(Height$Location)
table(Height$Rep)
table(Height$Crop)
table(Height$Cross, Height$Location)
table(Height$Crop, Height$Location, Height$PlantYear)
table(Height$Cross, Height$Location, Height$PlantYear)
table(Height$Cross, Height$Location, Height$PlantYear, Height$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
Height <- subset(Height, !(Cross %in% c("CP18-1611")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(Height$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
Heightharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", Heightharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
Height <- Height %>%
  mutate(Env = interaction(Location, PlantYear))
table(Height$Env)

Height$Inter <- interaction(Height$Cross, 
                            Height$Env)
Height$FxCrop <- interaction(Height$Cross, 
                             Height$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
Heightmod1 <- lme(fixed = Height ~ 1,  
               random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
               weights = varIdent(form = ~1 | Env),
               data = Height,
               method = 'REML')
summary(Heightmod1)
plot(Heightmod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(Heightmod1)
# Create an empty dataframe to store the variance components
DataVarComp$Height <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$Height[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$Height[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$Height[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$Height[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$Height[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$Height[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(Heightmod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(Heightmod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_Height <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_Height")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_Height"] <- "Height"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(Height = as.numeric(as.character(Height)))
# Filter and calculate result
HeightH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(Height, na.rm = TRUE) / sum(DataVarComp$Height, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- HeightH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$Height[DataVarComp$Component == "Heritability"] <- value_to_assign
#============Sucrose=============
# Load the 'car' package
table(Sucrose$Cross)
table(Sucrose$Location)
table(Sucrose$Rep)
table(Sucrose$Crop)
table(Sucrose$Cross, Sucrose$Location)
table(Sucrose$Crop, Sucrose$Location, Sucrose$PlantYear)
table(Sucrose$Cross, Sucrose$Location, Sucrose$PlantYear)
table(Sucrose$Cross, Sucrose$Location, Sucrose$PlantYear, Brix$Crop)
#upon evaluation of the dataframe for TRS, certain crosses need to be removed
#because they are not observed in every location
# Create a new dataframe without the unwanted rows
Sucrose <- subset(Sucrose, !(Cross %in% c("CP18-1611", "HB 17-3208")))
#the harmonic mean of the number of observations is needed to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(Sucrose$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
Sucroseharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", Sucroseharmonic_mean, "\n")
#Combine columns to create the 'Env' factor
Sucrose <- Sucrose %>%
  mutate(Env = interaction(Location, PlantYear))
table(Sucrose$Env)

Sucrose$Inter <- interaction(Sucrose$Cross, 
                             Sucrose$Env)
Sucrose$FxCrop <- interaction(Sucrose$Cross, 
                              Sucrose$Crop)
# Create the linear mixed model using lme function and specify heterogenous variance 
Sucrosemod1 <- lme(fixed = Sucrose.. ~ 1,  
                random = list(Cross = pdDiag(~1), Crop = pdDiag(~1), Env = pdDiag(~1), FxCrop = pdDiag(~1), Inter = pdDiag(~1)),
                weights = varIdent(form = ~1 | Env),
                data = Sucrose,
                method = 'REML')
summary(Sucrosemod1)
plot(Sucrosemod1)
# Obtain the variance components
# Extract the variance components from TRSlme2
variance_components <- VarCorr(Sucrosemod1)
# Create an empty dataframe to store the variance components
DataVarComp$Sucrose <- NA
# Check if variance components exist for "Cross"
if (!is.null(variance_components[2])) {
  DataVarComp$Sucrose[DataVarComp$Component == "Cross"] <- variance_components[2]
}
if (!is.null(variance_components[4])) {
  DataVarComp$Sucrose[DataVarComp$Component == "Crop"] <- variance_components[4]
}
if (!is.null(variance_components[6])) {
  DataVarComp$Sucrose[DataVarComp$Component == "Env"] <- variance_components[6]
}
if (!is.null(variance_components[8])) {
  DataVarComp$Sucrose[DataVarComp$Component == "FxCrop"] <- variance_components[8]
}
if (!is.null(variance_components[10])) {
  DataVarComp$Sucrose[DataVarComp$Component == "Inter"] <- variance_components[10]
}
if (!is.null(variance_components[11])) {
  DataVarComp$Sucrose[DataVarComp$Component == "Residual"] <- variance_components[11]
}
# Extract the fixed effects (mean) from TRSlme2
fixed_effects <- fixef(Sucrosemod1)
# Extract the BLUPs for each level of Cross from TRSlme2
BLUPs_Cross <- ranef(Sucrosemod1)$Cross
# Convert BLUPs to a dataframe and add the Cross column
BLUPs_df <- as.data.frame(BLUPs_Cross)
BLUPs_df$Cross <- rownames(BLUPs_df)
# Merge DataOutput with the BLUPs dataframe based on the "Cross" column
merged_data <- merge(DataOutput, BLUPs_df, by = "Cross")
# Calculate the predicted values by adding the fixed effects (mean) and BLUPs
merged_data$Predicted_Sucrose <- fixed_effects["(Intercept)"] + merged_data$`(Intercept)`
DataOutput <- merge(DataOutput, merged_data[, c("Cross", "Predicted_Sucrose")], by = "Cross", all.x = TRUE)
colnames(DataOutput)[colnames(DataOutput) == "Predicted_Sucrose"] <- "Sucrose"
# Convert TRS to numeric, handling non-numeric values
DataVarComp <- DataVarComp %>%
  mutate(Sucrose = as.numeric(as.character(Sucrose)))
# Filter and calculate result
SucroseH2 <- DataVarComp %>%
  filter(Component %in% c("Cross", "FxCrop")) %>%
  summarise(result_value = sum(Sucrose, na.rm = TRUE) / sum(DataVarComp$Sucrose, na.rm = TRUE))
# Assuming TRSH2 contains the calculated value
value_to_assign <- SucroseH2$result_value
# Update the specific cell in DataVarComp
DataVarComp$Sucrose[DataVarComp$Component == "Heritability"] <- value_to_assign
