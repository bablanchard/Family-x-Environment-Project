#=============QC of the data===========
library(readxl)

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










#=============For TRS======================
# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
# Count the number of unique values of each cross in each location
unique_counts <- aggregate(TRS ~ Cross + Location, data = TRS, function(x) length(unique(x)))
print(unique_counts)

# Check normality of the TRS variable and plot a histogram
shapiro.test(TRS$TRS)
hist(TRS$TRS)

# Plot a Q-Q plot of the Econ variable
qqnorm(TRS$TRS)
qqline(TRS$TRS, col = "red")

# Load the 'car' package
library(car)
table(TRS$Cross, TRS$Location)


#take out the two varieties that have not been tested in each location
#df_FB15M <- df_FB15M[!(df_FB15M$Variety %in% c("FLB15-0444", "FLB15-0652")), ]

#table(df_FB15M$Variety, df_FB15M$Loc)



library(lme4)
library(car)

# Fit the linear mixed model
fit1 <- lmer(TRS ~ 1 + Cross + (1|Location), data = TRS)
summary(fit1)
anova(fit1)

# Extract the residuals
resid <- residuals(fit1)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit1, resid = TRUE, sqrt(abs(resid(.))) ~ fitted(.))




library(emmeans)
library(ggplot2)
# Compute pairwise comparisons using emmeans
emmeans(fit1, pairwise ~ Cross, adjust = "tukey")

# Calculate EMMs and CIs for each level of the "Variety" factor
emm <- emmeans(fit1, "Cross", type = "response", infer = c(TRUE, TRUE))

# Get the EMMs and CIs as a data frame
emm_df <- as.data.frame(summary(emm))

# Define the significance level (alpha)
alpha <- 0.05

# Create a new column that designates significant differences
avg_emmean <- mean(emm_df$emmean, na.rm = TRUE)
emm_df$significant <- ifelse(emm_df$p.value < alpha/2 & emm_df$emmean > avg_emmean, "higher",
                             ifelse(emm_df$p.value < alpha/2 & emm_df$emmean < avg_emmean, "lower", "not significant"))



# plot the confidence intervals
ggplot(emm_df, aes(x = Cross, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  labs(x = "TRS", y = "Cross", title = "Variety TRS EMM Confidence Intervals") +
  coord_flip()

# Create a color palette for the groups
palette <- c("higher" = "blue", "lower" = "red", "not significant" = "black")


# Plot the LSMs with colored error bars
ggplot(emm_df, aes(x = emmean, y = Cross, xmin = lower.CL, xmax = upper.CL, color = significant)) +
  geom_point() +
  geom_errorbarh(height = 0.5) +
  scale_color_manual(values = palette) +
  labs(x = "TRS", y = "Cross", title = "Cross TRS EMM Confidence Intervals") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Plot the EMMs with error bars
plot(emm, xlab = "TRS", ylab = "Cross", ylim = c(0, NA), main = "Estimated Marginal Means",
     type = "response", main.args = list(cex = 0.8))



library(lsmeans)
library(emmeans)
library(ggplot2)


# Compute the LS means and CI for the Variety factor
lsmeans <- lsmeans(fit1, "Cross")
lsmeans_df <- as.data.frame(lsmeans)

# plot the confidence intervals
ggplot(lsmeans_df, aes(x = Cross, y = lsmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  labs(x = "TRS", y = "Cross", title = "Cross TRS LSM Confidence Intervals") +
  coord_flip()

#average lsmean
avg_lsmean <- mean(lsmeans_df$lsmean, na.rm = TRUE)

# Create a column to indicate if the Cross is statistically different from mean
lsmeans_df$diff_from_mean <- ifelse(lsmeans_df$lsmean - avg_lsmean > lsmeans_df$SE * qt(0.975, df = lsmeans_df$df), "Higher",
                                    ifelse(avg_lsmean - lsmeans_df$lsmean > lsmeans_df$SE * qt(0.975, df = lsmeans_df$df), "Lower", "Not different"))


# Create a color palette for the groups
palette <- c("Higher" = "blue", "Lower" = "red", "Not different" = "black")

# Plot the LSMs with colored error bars
ggplot(lsmeans_df, aes(x = lsmean, y = Cross, xmin = lower.CL, xmax = upper.CL, color = diff_from_mean)) +
  geom_point() +
  geom_errorbarh(height = 0.5) +
  scale_color_manual(values = palette) +
  labs(x = "TRS", y = "Cross", title = "Cross TRS LSM Confidence Intervals") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# Pairwise comparisons of the Variety factor using Tukey's method
pairs(lsmeans(fit1, "Cross"), adjust = "tukey")




install.packages("contrast")
library(tibble)
library(multcomp)
library(tidyverse)
library(broom)




boxplot(TRS$TRS, col = "red")


# outlier detection and elimination; similar to first class
#outlier <- names(outlierTest(TRS_model)$p)
#TRS[outlier, "TRS"] <- NA


#For TRS
TRS_model <- lmer(TRS ~ Cross + (1|Location) + Crop + PlantYear + (1|Rep:Location) + (1|Cross:Location), data = TRS)
emm_table <- emmeans(TRS_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger")
# Get the summary statistics of the EMMs
TRS_avgs <- summary(emm_table, infer = c(TRUE, TRUE, TRUE, FALSE))

# Convert the TRS_avgs object to a data frame
TRS_avgs <- as.data.frame(TRS_avgs)



TRS_emm <- emmeans(TRS_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger", data = TRS)
pairs(TRS_emm, adjust = "sidak")
TRS_mc <- glht(TRS_model, linfct = mcp(Cross = "Tukey"))
summary(TRS_mc, test = adjusted(type = "bonferroni"))
TRS_summary_table <- summary(TRS_mc, test = adjusted(type = "bonferroni"))
TRS_mc_df <- tidy(TRS_summary_table)
TRS_output <- capture.output(summary(TRS_mc, test = adjusted(type = "bonferroni")))
#write.csv(TRS_output, file = "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_comps.csv", row.names = FALSE)
#write.csv(TRS_avgs, "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_avgs.csv", row.names=FALSE)
#TRS_data <- read.csv("C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_comps.csv")
#TRS_subset_data <- subset(TRS_data, grepl("CP96-1252", TRS_data$x))
#colnames(TRS_subset_data)[1] <- 'TRS_comp'
#write.csv(TRS_subset_data, file = "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TSA_comps.csv", row.names = FALSE, append = TRUE)


#################################### MET analysis #########################
library(sommer)

blocks <- length(unique(TRS$Rep))
loc <- length(unique(TRS$Location))
crops <- length(unique(TRS$Crop))
years <- length(unique(TRS$PlantYear))

# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
TRSfitMET <- mmer(fixed = TRS ~ 1 + Location + Crop + PlantYear + Cross:Crop,
               random = ~Cross + Rep  + Cross:Location,
               rcov = ~ vsr(units),
               data = TRS)

summary(TRSfitMET)
# Broad-sense heritability
vpredict(TRSfitMET, h2 ~ V1 / ( V1 + V3/loc + V4/(loc*blocks*crops*years) ) ) # trials level; replicates is not part of the heritability equation
vpredict(TRSfitMET, h2 ~ V1 / ( V1 + V3 + V4 ) ) # plot level; Vs are the order in the model


# Fitting genotype by environment models - unstructured model (US)
TRSfitMET.US <- mmer(fixed = TRS ~ 1 + Location + Crop + PlantYear + Cross:Crop,
                  random = ~Cross + Rep + vsr(usr(Location), Cross),
                  rcov = ~ vsr(units),
                  data = TRS)
summary(TRSfitMET.US)

# Broad-sense heritability
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + V9/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + V9 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.H <- mmer(fixed = TRS ~ 1 + Location + Crop + PlantYear + Cross:Crop,
                    random = ~Cross + Rep + vsr(usr(Location), Cross),
                    rcov = ~ vsr(dsr(Location), units),
                    data = TRS)

summary(TRSfitMET.US.H)
# Broad-sense heritability
vpredict(TRSfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + (V9 + V10 + V11)/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + (V9 + V10 + V11) ) ) # plot level

# comparing the models
anova(TRSfitMET, TRSfitMET.US)
anova(TRSfitMET.US, TRSfitMET.US.H)
#TRSfitMET.US.H seems to perform best

# predicting BLUPs per environment
BLUPS <- data.frame(
  MET = TRSfitMET$U$Cross$TRS,
  US = TRSfitMET.US$U$Cross$TRS,
  USH = TRSfitMET.US.H$U$Cross)

head(BLUPS)
cor(BLUPS)

####################### Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)

library(statgenGxE)

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

# let's take a look at the output
names(dropsFW)

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")


############################# GGE-Biplot Analysis ######################################
library(metan)

model.gge <- gge(TRSminusone, Location, Cross, TRS, svp = "symmetrical")

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where

######## the end ##########

