library(readxl)
library("PerformanceAnalytics")
library(car)
require(bestNormalize)
library(lme4)
library(sommer)
library(dplyr)
library(statgenGxE)
library(ggplot2)
library(metan)







#=============QC of the data===========
#This block will read in the master data, take out New Iberia year 2 data, calculate stalk weight, 
#replace missing stool values with survival, use individual data to estimate family data,
#remove certain families with missing data, make a csv dataset that is clean,
#and calculate phenotypic correlations in a chart

#read master
df <- read.csv("C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/Dissertation Projects/Master.csv", header = TRUE)

#take out NI Y2
df_filtered <- df[!(df$PlantYear == 2021 & df$Location == "NI"), ]

#remove observations from master that have 0 or NA for all of the specified columns
#df_filtered <- df_filtered[rowSums(df_filtered[, c("Planted", "Survival", "PlotWeight", "TRS", "Bu.Wt")] == 0 | is.na(df_filtered[, c("Planted", "Survival", "PlotWeight", "TRS", "Bu.Wt")])) != length(c("Planted", "Survival", "PlotWeight", "TRS", "Bu.Wt")),]
#df_filtered <- df_filtered[rowSums(df_filtered[, c("Stool", "Survival", "PlotWeight", "TRS", "Bu.Wt")] == 0 | is.na(df_filtered[, c("Stool", "Survival", "PlotWeight", "TRS", "Bu.Wt")])) != length(c("Stool", "Survival", "PlotWeight", "TRS", "Bu.Wt")),]


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
#install.packages("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(df_filtered[11:22])), histogram = TRUE, pch = 1)


#===========TRS Analysis==============
# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
# Count the number of unique values of each cross in each location
unique_counts <- aggregate(TRS ~ Cross + Location, data = TRS, function(x) length(unique(x)))
print(unique_counts)

# Load the 'car' package
table(TRS$Cross)
table(TRS$Location)
table(TRS$Rep)
table(TRS$Crop)
table(TRS$Cross, TRS$Location)
table(TRS$Cross, TRS$Location, TRS$PlantYear)
table(TRS$Cross, TRS$Location, TRS$PlantYear, TRS$Crop)


# Check normality of the TRS variable and plot a histogram
shapiro.test(TRS$TRS)
hist(TRS$TRS)

# Plot a Q-Q plot of the Econ variable
qqnorm(TRS$TRS)
qqline(TRS$TRS, col = "red")

# testing for normality
# First lets check using patterns
shapiro.test(rnorm(length(TRS$TRS))) # normal distribution
shapiro.test(runif(length(TRS$TRS))) # uniform distribution
# then, 
shapiro.test(TRS$TRS)

#install.packages("bestNormalize")
TRSadj <- bestNormalize(TRS$TRS, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
TRSadj$chosen_transform
shapiro.test(TRSadj$x.t)
TRS$TRSadj <- TRSadj$x.t
head(TRS)

# Define the mixed model formula
# Use () to nest Rep within Location, and Crop within PlantYear
# Use : to specify interactions between Cross, Location, and Crop
# Use | to specify that all factors are random effects
# Use (1 | ) to include an intercept
formula1 <- TRS ~ (1) + (1 | Cross) + (1 | Rep:Location) + (1 | Location)  +  (1 | Crop:PlantYear) + 
  (1 | Cross:Location) + (1 | Cross:Crop:Location)
# Fit the mixed model using the lmer function from the lme4 package
fit1 <- lmer(formula, data = TRS, REML = TRUE)
summary(fit1)

#model 2
formula2 <- TRS ~ (1) + (1 | Cross) + (1 | Rep:Location) + (1 | Location)  +  (1 | Crop:PlantYear) + 
  (1 | Cross:Location)
# Fit the mixed model using the lmer function from the lme4 package
fit2 <- lmer(formula2, data = TRS, REML = TRUE)
summary(fit2)

# comparing the models
anova(fit1, fit2)

#more models using sommer package
blocks <- length(unique(TRS$Rep))
loc <- length(unique(TRS$Location))
crops <- length(unique(TRS$Crop))
years <- length(unique(TRS$PlantYear))

#model3
# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
TRSfitMET <- mmer(fixed = TRS ~ 1,
                  random = ~Cross + Location + Crop + Rep + PlantYear + Cross:Crop + Cross:Location,
                  rcov = ~ vsr(units),
                  data = TRS)
summary(TRSfitMET)

# model4 Fitting genotype by environment models - unstructured model (US)
TRSfitMET.US <- mmer(fixed = TRS ~ 1,
                     random = ~Cross+ Location + Crop + Rep + PlantYear + Cross:Crop + vsr(usr(Location), Cross),
                     rcov = ~ vsr(units),
                     data = TRS)
summary(TRSfitMET.US)

# model5 Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.H <- mmer(fixed = TRS ~ 1,
                       random = ~Cross + Location + Crop + Rep + PlantYear + Cross:Crop + vsr(usr(Location), Cross),
                       rcov = ~ vsr(dsr(Location), units),
                       data = TRS)
summary(TRSfitMET.US.H)

# model6 Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.H2 <- mmer(fixed = TRS ~ 1,
                        random = ~Cross + Location + Crop + Rep + PlantYear + vsr(usr(Location), Cross),
                        rcov = ~ vsr(dsr(Location), units),
                        data = TRS)
summary(TRSfitMET.US.H2)

# comparing the models
anova(TRSfitMET, TRSfitMET.US)
anova(TRSfitMET.US, TRSfitMET.US.H)
anova(TRSfitMET.US.H, TRSfitMET.US.H2)

#model6 seems to be the best so use for analysis
# Extract variance components
vc <- summary(TRSfitMET.US.H2)$varcomp

# Calculate trial-level heritability
trial_vc <- vc["Cross.TRS-TRS", "VarComp"]
NIR <- vc["NI:units.TRS-TRS", "VarComp"]
NRR <- vc["NR:units.TRS-TRS", "VarComp"]
SGR <- vc["SG:units.TRS-TRS", "VarComp"]
trial_h2 <- trial_vc / (trial_vc + NIR + NRR + SGR)
cat("Trial-level heritability:", trial_h2, "\n")

# Calculate plot-level heritability
plotNIvc <- vc["NI:Cross.TRS-TRS", "VarComp"]
plotNRvc <- vc["NR:Cross.TRS-TRS", "VarComp"]
plotSGvc <- vc["SG:Cross.TRS-TRS", "VarComp"]
plot_vc <- plotSGvc + plotNRvc + plotNIvc
plot_h2 <- plot_vc / (plot_vc +  NIR + NRR + SGR)
cat("Plot-level heritability:", plot_h2, "\n")

# predicting BLUPs 
BLUPS <- data.frame(
  MET = TRSfitMET$U$Cross$TRS,
  US = TRSfitMET.US$U$Cross$TRS,
  USH = TRSfitMET.US.H$U$Cross$TRS,
  USH2 = TRSfitMET.US.H2$U$Cross$TRS)

head(BLUPS)
cor(BLUPS)



#modify the BLUPS dataframe
BLUPS$Cross <- rownames(BLUPS)
BLUPS$Cross <- sub("^Cross", "", BLUPS$Cross)
BLUPS <- dplyr::rename(BLUPS, BLUP_USH2 = USH2)
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(BLUPS$BLUP_USH2)
blup_sd <- sd(BLUPS$BLUP_USH2)
# Create a new column for highlighting
BLUPS$SignificantDiff <- ifelse(BLUPS$BLUP_USH2 > (blup_mean + blup_sd), "Greater", ifelse(BLUPS$BLUP_USH2 < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
BLUPS <- BLUPS[order(BLUPS$BLUP_USH2),]
# Create a bar plot with highlighted crosses and their BLUPs
ggplot(BLUPS, aes(x=reorder(Cross, BLUP_USH2), y=BLUP_USH2, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses Across All Locations", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))




# predicting the overal performance
predict.mmer(TRSfitMET.US.H2, D="Cross")$pvals

BLUPs <- predict.mmer(object = TRSfitMET.US.H2, classify = "Cross")$pvals



# to predict BLUPs per environment, we need to sum the main effect and the interaction for that specific case 
# For instance, G in ideal N
BLUPPERENV <- data.frame(BLUPNI = TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`NI:Cross`$TRS)
BLUPPERENV$BLUPNR <- TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`NR:Cross`$TRS
BLUPPERENV$BLUPSG <- TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`SG:Cross`$TRS
cor(BLUPPERENV)

BLUPPERENV$Cross <- gsub("^Cross", "", rownames(BLUPPERENV))
# Merge BLUPPERENV and BLUPS by matching values in the Cross column
results <- merge(BLUPPERENV, BLUPS, by = "Cross")

# Rename BLUP_USH2 as BLUPALL
names(results)[names(results) == "BLUP_USH2"] <- "BLUPALL"

# Keep only desired columns
results <- results[, c("Cross", "BLUPALL", "BLUPNI", "BLUPNR", "BLUPSG")]

#correlations of BLUPs
chart.Correlation(as.matrix(na.omit(results[2:5])), histogram = TRUE, pch = 1)

#bar plot for NI
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(results$BLUPNI)
blup_sd <- sd(results$BLUPNI)
# Create a new column for highlighting
results$SignificantDiff <- ifelse(results$BLUPNI > (blup_mean + blup_sd), "Greater", ifelse(results$BLUPNI < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
results <- results[order(results$BLUPNI),]
# Create a bar plot with highlighted crosses and their BLUPs
ggplot(results, aes(x=reorder(Cross, BLUPNI), y=BLUPNI, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses in New Iberia", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

#bar plot for NR
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(results$BLUPNR)
blup_sd <- sd(results$BLUPNR)
# Create a new column for highlighting
results$SignificantDiff <- ifelse(results$BLUPNR > (blup_mean + blup_sd), "Greater", ifelse(results$BLUPNR < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
results <- results[order(results$BLUPNR),]
# Create a bar plot with highlighted crosses and their BLUPs
ggplot(results, aes(x=reorder(Cross, BLUPNR), y=BLUPNR, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses in New Roads", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

#bar plot for SG
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(results$BLUPSG)
blup_sd <- sd(results$BLUPSG)
# Create a new column for highlighting
results$SignificantDiff <- ifelse(results$BLUPSG > (blup_mean + blup_sd), "Greater", ifelse(results$BLUPSG < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
results <- results[order(results$BLUPSG),]
# Create a bar plot with highlighted crosses and their BLUPs
ggplot(results, aes(x=reorder(Cross, BLUPSG), y=BLUPSG, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses in St. Gabriel", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))

#################################### MET analysis #########################
####################### TRS Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)


## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = TRS, genotype = "Cross", trial = "Location")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "TRS")

TRSminusone <- TRS[TRS$Cross != "CP18-1611", ]#remove one cross because only present in one location

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = TRSminusone, genotype = "Cross", trial = "Location")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "TRS")
summary(dropsFW)

## Create a box plot of dropsTD.
## Color the boxes based on the variable scenarioFull.
## Plot in  descending order.
plot(dropsTD, plotType = "box", traits = "TRS", colorTrialBy = "trial",
     orderBy = "descending")

# let's take a look at the output
names(dropsFW)

## Create a scatter plot of dropsTD.
## Color the genotypes based on the variable geneticGroup.
## Color the histograms for trials based on the variable scenarioFull.
plot(dropsTD, plotType = "scatter", traits = "TRS", 
     colorTrialBy = "trial")

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

# Create an empty column 'GenMean' in results
results$GenMean <- NA
results$SE_GenMean <- NA
results$Rank <- NA
results$Sens <- NA
results$SE_Sens <- NA
results$MSdeviation <- NA

# Loop through each row of results
for (i in 1:nrow(results)) {
  
  # Find the index of the matching genotype in dropsFW$estimates
  idx <- match(results$Cross[i], dropsFW$estimates$Genotype)
  
  # If a matching genotype is found, set the corresponding 'GenMean' value in results
  if (!is.na(idx)) {
    results$GenMean[i] <- dropsFW$estimates$GenMean[idx]
    results$SE_GenMean[i] <- dropsFW$estimates$SE_GenMean[idx]
    results$Rank[i] <- dropsFW$estimates$Rank[idx]
    results$Sens[i] <- dropsFW$estimates$Sens[idx]
    results$SE_Sens[i] <- dropsFW$estimates$SE_Sens[idx]
    results$MSdeviation[i] <- dropsFW$estimates$MSdeviation[idx]
    
  }
}

############################# TRS GGE-Biplot Analysis ######################################

model.gge <- gge(TRSminusone, Location, Cross, TRS, svp = "symmetrical")

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where
(d <- plot(model.gge, type = 4)) #Discriminativeness vs. representativeness
(e <- plot(model.gge, type = 5, sel_env = "NI")) #examine an environment
(f <- plot(model.gge, type = 6)) #ranking environments
(g <- plot(model.gge, type = 7, sel_gen = "CP18-1041")) #examine a geno
(h <- plot(model.gge, type = 8)) #ranking genos
(i <- plot(model.gge, type = 9, sel_gen1 = "XL17-226", sel_gen2 = "CP18-1041")) #compare two genos
(j <- plot(model.gge, type = 10)) #relationship among environments


#===========Multi-Trait Selection using BLUPs=================use harmonic mean for multiple traits
#Multi-"Trait" Selection for BLUPs in NI
rank_BLUPALL <- rank(results$BLUPALL)
rank_BLUPNI <- rank(results$BLUPNI)
rank_cutoff_BLUPALL <- quantile(rank_BLUPALL, 0.85)
rank_cutoff_BLUPNI <- quantile(rank_BLUPNI, 0.85)
top_15_percent <- results[(rank_BLUPALL >= rank_cutoff_BLUPALL) & (rank_BLUPNI >= rank_cutoff_BLUPNI), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- results$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_BLUPALL >= rank_cutoff_BLUPALL & rank_BLUPNI >= rank_cutoff_BLUPNI, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(results$BLUPALL)
blupall_sd <- sd(results$BLUPALL)
blupNI_mean <- mean(results$BLUPNI)
blupNI_sd <- sd(results$BLUPNI)

#find high outliers
blupall_outliers <- results$BLUPALL > (blupall_mean + blupall_sd)
blupNI_outliers <- results$BLUPNI > (blupNI_mean + blupNI_sd)
blupall_outlier_cross <- results$Cross[blupall_outliers]
blupNI_outlier_cross <- results$Cross[blupNI_outliers]
print(blupall_outlier_cross)
print(blupNI_outlier_cross)


point_colors <- ifelse(results$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(results$Cross %in% results$Cross[blupall_outliers | blupNI_outliers], "#E67E22", "#2E86C1"))




# Scatterplot with colors modified for outliers and a gridded background
ggplot(results, aes(x = BLUPALL, y = BLUPNI, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Correlation between BLUPALL and BLUPNI",
       x = "BLUPALL",
       y = "BLUPNI") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)

#Multi-"Trait" Selection for BLUPs in NR
rank_BLUPALL <- rank(results$BLUPALL)
rank_BLUPNR <- rank(results$BLUPNR)
rank_cutoff_BLUPALL <- quantile(rank_BLUPALL, 0.85)
rank_cutoff_BLUPNR <- quantile(rank_BLUPNR, 0.85)
top_15_percent <- results[(rank_BLUPALL >= rank_cutoff_BLUPALL) & (rank_BLUPNR >= rank_cutoff_BLUPNR), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- results$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_BLUPALL >= rank_cutoff_BLUPALL & rank_BLUPNR >= rank_cutoff_BLUPNR, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(results$BLUPALL)
blupall_sd <- sd(results$BLUPALL)
blupNR_mean <- mean(results$BLUPNR)
blupNR_sd <- sd(results$BLUPNR)

#find high outliers
blupall_outliers <- results$BLUPALL > (blupall_mean + blupall_sd)
blupNR_outliers <- results$BLUPNR > (blupNR_mean + blupNR_sd)
blupall_outlier_cross <- results$Cross[blupall_outliers]
blupNR_outlier_cross <- results$Cross[blupNR_outliers]
print(blupall_outlier_cross)
print(blupNR_outlier_cross)


point_colors <- ifelse(results$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(results$Cross %in% results$Cross[blupall_outliers | blupNR_outliers], "#E67E22", "#2E86C1"))




# Scatterplot with colors modified for outliers and a gridded background
ggplot(results, aes(x = BLUPALL, y = BLUPNR, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Correlation between BLUPALL and BLUPNR",
       x = "BLUPALL",
       y = "BLUPNR") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)


#Multi-"Trait" Selection for BLUPs in SG
rank_BLUPALL <- rank(results$BLUPALL)
rank_BLUPSG <- rank(results$BLUPSG)
rank_cutoff_BLUPALL <- quantile(rank_BLUPALL, 0.85)
rank_cutoff_BLUPSG <- quantile(rank_BLUPSG, 0.85)
top_15_percent <- results[(rank_BLUPALL >= rank_cutoff_BLUPALL) & (rank_BLUPSG >= rank_cutoff_BLUPSG), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- results$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_BLUPALL >= rank_cutoff_BLUPALL & rank_BLUPSG >= rank_cutoff_BLUPSG, "seagreen3", "royalblue3")


#calculate mean and sd BLUPs
blupall_mean <- mean(results$BLUPALL)
blupall_sd <- sd(results$BLUPALL)
blupSG_mean <- mean(results$BLUPSG)
blupSG_sd <- sd(results$BLUPSG)

#find high outliers
blupall_outliers <- results$BLUPALL > (blupall_mean + blupall_sd)
blupSG_outliers <- results$BLUPSG > (blupSG_mean + blupSG_sd)
blupall_outlier_cross <- results$Cross[blupall_outliers]
blupSG_outlier_cross <- results$Cross[blupSG_outliers]
print(blupall_outlier_cross)
print(blupSG_outlier_cross)


point_colors <- ifelse(results$Cross %in% top_15_percent$Cross, "#1B4F72",
                       ifelse(results$Cross %in% results$Cross[blupall_outliers | blupSG_outliers], "#E67E22", "#2E86C1"))




# Scatterplot with colors modified for outliers and a gridded background
ggplot(results, aes(x = BLUPALL, y = BLUPSG, label = Cross, color = point_colors)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(values = c("seagreen3", "#2E86C1", "#E67E22"),
                     labels = c("Top 15%", "Bottom 85%", "High Outliers"),
                     name = "Quantile") +
  labs(title = "Correlation between BLUPALL and BLUPSG",
       x = "BLUPALL",
       y = "BLUPSG") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 18),
        legend.position = "bottom",
        panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "gray70", size = 0.2),
        panel.border = element_blank()) +
  geom_text(aes(color = point_colors, col = "black"), check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)
