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



# phenotypic correlation between traits
sapply(df_filtered[, 11:22, 27], class)
df_filtered[, 11:22] <- lapply(df_filtered[, 11:22], as.numeric)
df_filtered[, 27] <- as.numeric(df_filtered[, 27])
round(cor(df_filtered[, c(11:22, 27)], use = "pairwise"), 2)
#install.packages("PerformanceAnalytics")
library("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(df_filtered[11:22])), histogram = TRUE, pch = 1)







#=============For TRS======================

# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
# Count the number of unique values of each cross in each location
unique_counts <- aggregate(TRS ~ Cross + Location, data = TRS, function(x) length(unique(x)))
print(unique_counts)

# Load the 'car' package
library(car)
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
require(bestNormalize)
TRSadj <- bestNormalize(TRS$TRS, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
TRSadj$chosen_transform
shapiro.test(TRSadj$x.t)
TRS$TRSadj <- TRSadj$x.t
head(TRS)

#take out the two varieties that have not been tested in each location
#df_FB15M <- df_FB15M[!(df_FB15M$Variety %in% c("FLB15-0444", "FLB15-0652")), ]

#table(df_FB15M$Variety, df_FB15M$Loc)


library(lme4)
library(car)

# Define the mixed model formula
# Use () to nest Rep within Location, and Crop within PlantYear
# Use : to specify interactions between Cross, Location, and Crop
# Use | to specify that all factors are random effects
# Use (1 | ) to include an intercept
formula <- TRS ~ (1) + (1 | Cross) + (1 | Rep:Location) + (1 | Location)  +  (1 | Crop:PlantYear) + 
  (1 | Cross:Location) + (1 | Cross:Crop:Location)

# Fit the mixed model using the lmer function from the lme4 package
fit1 <- lmer(formula, data = TRS, REML = TRUE)

# Print the model summary
summary(fit1)
anova(fit1)

# Extract the residuals
resid <- residuals(fit1)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit1, resid = TRUE, sqrt(abs(resid(.))) ~ fitted(.))

#model 2
formula2 <- TRS ~ (1) + (1 | Cross) + (1 | Rep:Location) + (1 | Location)  +  (1 | Crop:PlantYear) + 
  (1 | Cross:Location)

# Fit the mixed model using the lmer function from the lme4 package
fit2 <- lmer(formula2, data = TRS, REML = TRUE)

# Print the model summary
summary(fit2)
anova(fit2)

# Extract the residuals
resid <- residuals(fit2)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit2, resid = TRUE, sqrt(abs(resid(.))) ~ fitted(.))


# comparing the models
anova(fit1, fit2)
#no significant difference so use fit1

# Trial level heritability
Var_G <- as.numeric(VarCorr(fit1)[3])  # genetic variance
Var_GE <- as.numeric(VarCorr(fit1)[2]) + as.numeric(VarCorr(fit1)[1]) + as.numeric(VarCorr(fit1)[4]) + as.numeric(VarCorr(fit1)[5]) + as.numeric(VarCorr(fit1)[6])# genotype-by-environment variance
Var_residual <- sigma(fit1)^2  # residual variance
e <- 3  # number of environments
r <- 2  # number of replicates per environment
H2_trial <- Var_G / (Var_G + Var_GE/e + Var_residual/e*r)

# Plot level heritability
n <- 2  # number of plots per genotype per environment
H2_plot <- Var_G / (Var_G + Var_residual/n)

# Print the results
cat("Trial level heritability: ", H2_trial, "\n")
cat("Plot level heritability: ", H2_plot, "\n")

#################################### MET analysis #########################
library(sommer)

blocks <- length(unique(TRS$Rep))
loc <- length(unique(TRS$Location))
crops <- length(unique(TRS$Crop))
years <- length(unique(TRS$PlantYear))

# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
TRSfitMET <- mmer(fixed = TRS ~ 1,
                  random = ~Cross + Location + Crop + Rep + PlantYear + Cross:Crop + Cross:Location,
                  rcov = ~ vsr(units),
                  data = TRS)

summary(TRSfitMET)
# Broad-sense heritability
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + V9/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V2 + V3 + V4) + V5 + V6 + V7 + V8+ V9 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US)
TRSfitMET.US <- mmer(fixed = TRS ~ 1,
                     random = ~Cross+ Location + Crop + Rep + PlantYear + Cross:Crop + vsr(usr(Location), Cross),
                     rcov = ~ vsr(units),
                     data = TRS)
summary(TRSfitMET.US)

# Broad-sense heritability
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + V9/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + V9 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.H <- mmer(fixed = TRS ~ 1,
                       random = ~Cross + Location + Crop + Rep + PlantYear + Cross:Crop + vsr(usr(Location), Cross),
                       rcov = ~ vsr(dsr(Location), units),
                       data = TRS)

summary(TRSfitMET.US.H)
# Broad-sense heritability
vpredict(TRSfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + (V9 + V10 + V11)/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + (V9 + V10 + V11) ) ) # plot level

# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.H2 <- mmer(fixed = TRS ~ 1,
                       random = ~Cross + Location + Crop + Rep + PlantYear + vsr(usr(Location), Cross),
                       rcov = ~ vsr(dsr(Location), units),
                       data = TRS)

summary(TRSfitMET.US.H2)
# Broad-sense heritability
vpredict(TRSfitMET.US.H2, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + (V9 + V10 + V11)/(loc*blocks*crops*years) ) ) # trials level
vpredict(TRSfitMET.US.H2, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + (V9 + V10 + V11) ) ) # plot level




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


# comparing the models
anova(TRSfitMET, TRSfitMET.US)
anova(TRSfitMET.US, TRSfitMET.US.H)
anova(TRSfitMET.US.H, TRSfitMET.US.H2)

#TRSfitMET.US.H2 seems to perform best





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
# Create a bar plot with highlighted crosses
library(ggplot2)
ggplot(BLUPS, aes(x=reorder(Cross, BLUP_USH2), y=BLUP_USH2, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
theme(axis.text.x = element_text(angle=90, hjust=1))




# predicting the overal performance
predict.mmer(TRSfitMET.US.H2, D="Cross")$pvals


# to predict BLUPs per environment, we need to sum the main effect and the interaction for that specific case 
# For instance, G in ideal N
BLUPPERENV <- data.frame(BLUPNI = TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`NI:Cross`$TRS)
BLUPPERENV$BLUPNR <- TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`NR:Cross`$TRS
BLUPPERENV$BLUPSG <- TRSfitMET.US.H2$U$Cross$TRS + TRSfitMET.US.H2$U$`SG:Cross`$TRS
cor(BLUPPERENV)




library(dplyr)

# assuming your data frames are called "results" and "BLUPS"
merged_df <- merge(results, BLUPS, by = "Cross") %>%
  select(Cross, BLUP_USH2) %>%
  rename(BLUPALL = BLUP_USH2)

# assign the updated data frame back to the "results" variable
results <- merged_df





# Fitting genotype by environment models - trying more models
TRSfitMET2 <- mmer(fixed = TRS ~ 1 + Location + PlantYear + Crop,
                        random = ~Male + Female + Cross:Location,
                        rcov = ~ vsr(units),
                        data = TRS)

summary(TRSfitMET2)






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

####################### TRS Stability and adaptability - Finlay-Wilkinson #########################
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
library(metan)

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



#=============For TRS (only genetic effects)======================
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
table(TRS$Cross)


#take out the two varieties that have not been tested in each location
#df_FB15M <- df_FB15M[!(df_FB15M$Variety %in% c("FLB15-0444", "FLB15-0652")), ]

#table(df_FB15M$Variety, df_FB15M$Loc)



library(lme4)
library(car)

# Fit the linear mixed model
fit1.2 <- lm(TRS ~ 1 + Cross, data = TRS)
summary(fit1.2)
anova(fit1.2)

# Extract the residuals
resid <- residuals(fit1.2)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit1.2, which = 1, resid = sqrt(abs(resid(fit1.2))) ~ fitted(fit1.2))




library(emmeans)
library(ggplot2)
# Compute pairwise comparisons using emmeans
emmeans(fit1.2, pairwise ~ Cross, adjust = "tukey")

# Calculate EMMs and CIs for each level of the "Variety" factor
emm <- emmeans(fit1.2, "Cross", type = "response", infer = c(TRUE, TRUE))

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
lsmeans <- lsmeans(fit1.2, "Cross")
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
pairs(lsmeans(fit1.2, "Cross"), adjust = "tukey")




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
TRS_model2 <- lm(TRS ~ 1 + Cross, data = TRS)
emm_table <- emmeans(TRS_model2, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger")
# Get the summary statistics of the EMMs
TRS_avgs <- summary(emm_table, infer = c(TRUE, TRUE, TRUE, FALSE))

# Convert the TRS_avgs object to a data frame
TRS_avgs <- as.data.frame(TRS_avgs)



TRS_emm <- emmeans(TRS_model2, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger", data = TRS)
pairs(TRS_emm, adjust = "sidak")
TRS_mc <- glht(TRS_model2, linfct = mcp(Cross = "Tukey"))
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

####################### TRS Stability and adaptability - Finlay-Wilkinson #########################
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


############################# TRS GGE-Biplot Analysis ######################################
library(metan)

model.gge <- gge(TRSminusone, Location, Cross, TRS, svp = "symmetrical")

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where
#=============For Stubbling======================
# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
Stool <- df_filtered[complete.cases(df_filtered$Stool) & df_filtered$Stool != 0,]
# Count the number of unique values of each cross in each location
unique_counts <- aggregate(Stool ~ Cross + Location, data = Stool, function(x) length(unique(x)))
print(unique_counts)

# Check normality of the Stool variable and plot a histogram
shapiro.test(Stool$Stool)
hist(Stool$Stool)

# Plot a Q-Q plot of the Stool variable
qqnorm(Stool$Stool)
qqline(Stool$Stool, col = "red")

# Load the 'car' package
library(car)
table(Stool$Cross, Stool$Location)


#take out the two varieties that have not been tested in each location
#df_FB15M <- df_FB15M[!(df_FB15M$Variety %in% c("FLB15-0444", "FLB15-0652")), ]

#table(df_FB15M$Variety, df_FB15M$Loc)



library(lme4)
library(car)

# Fit the linear mixed model
fit2 <- lmer(Stool ~ 1 + Cross + (1|Location), data = Stool)
summary(fit2)
anova(fit2)

# Extract the residuals
resid <- residuals(fit2)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit2, resid = TRUE, sqrt(abs(resid(.))) ~ fitted(.))




library(emmeans)
library(ggplot2)
# Compute pairwise comparisons using emmeans
emmeans(fit2, pairwise ~ Cross, adjust = "tukey")

# Calculate EMMs and CIs for each level of the "Variety" factor
emm <- emmeans(fit2, "Cross", type = "response", infer = c(TRUE, TRUE))

# Get the EMMs and CIs as a data frame
emm_df <- as.data.frame(summary(emm))

# Define the significance level (alpha)
alpha <- 0.05

# Create a new column that designates significant differences
avg_emmean <- mean(emm_df$emmean, na.rm = TRUE)
emm_df$significant <- ifelse(emm_df$p.value < alpha/2 & emm_df$emmean > avg_emmean, "higher",
                             ifelse(emm_df$p.value < alpha/2 & emm_df$emmean < avg_emmean, "lower", "not significant"))



# plot the confidence intervals
ggplot(emm_df, aes(x = Stool, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  labs(x = "Stool", y = "Cross", title = "Variety Stool EMM Confidence Intervals") +
  coord_flip()

# Create a color palette for the groups
palette <- c("higher" = "blue", "lower" = "red", "not significant" = "black")


# Plot the LSMs with colored error bars
ggplot(emm_df, aes(x = emmean, y = Cross, xmin = lower.CL, xmax = upper.CL, color = significant)) +
  geom_point() +
  geom_errorbarh(height = 0.5) +
  scale_color_manual(values = palette) +
  labs(x = "Stool", y = "Cross", title = "Cross Stool EMM Confidence Intervals") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Plot the EMMs with error bars
plot(emm, xlab = "Stool", ylab = "Cross", ylim = c(0, NA), main = "Estimated Marginal Means",
     type = "response", main.args = list(cex = 0.8))



library(lsmeans)
library(emmeans)
library(ggplot2)


# Compute the LS means and CI for the Variety factor
lsmeans <- lsmeans(fit2, "Cross")
lsmeans_df <- as.data.frame(lsmeans)

# plot the confidence intervals
ggplot(lsmeans_df, aes(x = Cross, y = lsmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_pointrange() +
  labs(x = "Stool", y = "Cross", title = "Cross Stool LSM Confidence Intervals") +
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
  labs(x = "Stool", y = "Cross", title = "Cross Stool LSM Confidence Intervals") +
  coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



# Pairwise comparisons of the Variety factor using Tukey's method
pairs(lsmeans(fit2, "Cross"), adjust = "tukey")




library(tibble)
library(multcomp)
library(tidyverse)
library(broom)




boxplot(Stool$Stool, col = "red")


# outlier detection and elimination; similar to first class
#outlier <- names(outlierTest(Stool_model)$p)
#Stool[outlier, "Stool"] <- NA


#For TRS
Stool_model <- lmer(Stool ~ Cross + (1|Location) + Crop + PlantYear + (1|Rep:Location) + (1|Cross:Location), data = Stool)
emm_table <- emmeans(Stool_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger")
# Get the summary statistics of the EMMs
Stool_avgs <- summary(emm_table, infer = c(TRUE, TRUE, TRUE, FALSE))

# Convert the Stool_avgs object to a data frame
Stool_avgs <- as.data.frame(Stool_avgs)



Stool_emm <- emmeans(Stool_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger", data = Stool)
pairs(Stool_emm, adjust = "sidak")
Stool_mc <- glht(Stool_model, linfct = mcp(Cross = "Tukey"))
summary(Stool_mc, test = adjusted(type = "bonferroni"))
Stool_summary_table <- summary(Stool_mc, test = adjusted(type = "bonferroni"))
Stool_mc_df <- tidy(Stool_summary_table)
Stool_output <- capture.output(summary(Stool_mc, test = adjusted(type = "bonferroni")))
#write.csv(TRS_output, file = "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_comps.csv", row.names = FALSE)
#write.csv(TRS_avgs, "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_avgs.csv", row.names=FALSE)
#TRS_data <- read.csv("C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TRS_comps.csv")
#TRS_subset_data <- subset(TRS_data, grepl("CP96-1252", TRS_data$x))
#colnames(TRS_subset_data)[1] <- 'TRS_comp'
#write.csv(TRS_subset_data, file = "C:/Users/F851129/OneDrive - Florida Crystals Corporation/Desktop/TSA_comps.csv", row.names = FALSE, append = TRUE)


#################################### MET analysis #########################
library(sommer)

blocks <- length(unique(Stool$Rep))
loc <- length(unique(Stool$Location))
crops <- length(unique(Stool$Crop))
years <- length(unique(Stool$PlantYear))

# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
StoolfitMET <- mmer(fixed = Stool ~ 1 + Location + Crop + PlantYear + Cross:Crop,
                  random = ~Cross + Rep  + Cross:Location,
                  rcov = ~ vsr(units),
                  data = Stool)

summary(StoolfitMET)
# Broad-sense heritability
vpredict(StoolfitMET, h2 ~ V1 / ( V1 + V3/loc + V4/(loc*blocks*crops*years) ) ) # trials level; replicates is not part of the heritability equation
vpredict(StoolfitMET, h2 ~ V1 / ( V1 + V3 + V4 ) ) # plot level; Vs are the order in the model


# Fitting genotype by environment models - unstructured model (US)
StoolfitMET.US <- mmer(fixed = Stool ~ 1 + Location + Crop + PlantYear + Cross:Crop,
                     random = ~Cross + Rep + vsr(usr(Location), Cross),
                     rcov = ~ vsr(units),
                     data = Stool)
summary(StoolfitMET.US)

# Broad-sense heritability
vpredict(StoolfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + V9/(loc*blocks*crops*years) ) ) # trials level
vpredict(StoolfitMET.US, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + V9 ) ) # plot level


# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
StoolfitMET.US.H <- mmer(fixed = Stool ~ 1 + Location + Crop + PlantYear + Cross:Crop,
                       random = ~Cross + Rep + vsr(usr(Location), Cross),
                       rcov = ~ vsr(dsr(Location), units),
                       data = Stool)

summary(StoolfitMET.US.H)
# Broad-sense heritability
vpredict(StoolfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7)/loc + (V9 + V10 + V11)/(loc*blocks*crops*years) ) ) # trials level
vpredict(StoolfitMET.US.H, h2 ~ V1 / ( V1 + (V4 + V6 + V7) + (V9 + V10 + V11) ) ) # plot level

# comparing the models
anova(StoolfitMET, StoolfitMET.US)
anova(StoolfitMET.US, StoolfitMET.US.H)
#

# predicting BLUPs per environment
BLUPS <- data.frame(
  MET = StoolfitMET$U$Cross$TRS,
  US = StoolfitMET.US$U$Cross$TRS,
  USH = StoolfitMET.US.H$U$Cross)

head(BLUPS)
cor(BLUPS)

####################### Stubbling Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)

library(statgenGxE)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = Stool, genotype = "Cross", trial = "Location")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "Stool")

Stoolminusone <- Stool[Stool$Cross != "CP18-1611", ]#remove one cross because only present in one location

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = Stoolminusone, genotype = "Cross", trial = "Location")

## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "Stool")
summary(dropsFW)

# let's take a look at the output
names(dropsFW)

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")


############################# Stubbling GGE-Biplot Analysis ######################################
library(metan)

model.gge <- gge(Stool, Location, Cross, Stool, svp = "symmetrical")

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where

#=============For TRS for New Iberia======================

# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS <- df_filtered[complete.cases(df_filtered$TRS) & df_filtered$TRS != 0,]
NI <- TRS[TRS$Location == "NI",]
# Count the number of unique values of each cross in each location
unique_counts <- aggregate(TRS ~ Cross + Location, data = NI, function(x) length(unique(x)))
print(unique_counts)

# Load the 'car' package
library(car)
table(NI$Cross)
table(NI$Rep)
table(NI$Crop)
table(NI$Cross, NI$Crop)


# Check normality of the TRS variable and plot a histogram
shapiro.test(NI$TRS)
hist(NI$TRS)

# Plot a Q-Q plot of the Econ variable
qqnorm(NI$TRS)
qqline(NI$TRS, col = "red")

# testing for normality
# First lets check using patterns
shapiro.test(rnorm(length(NI$TRS))) # normal distribution
shapiro.test(runif(length(NI$TRS))) # uniform distribution
# then, 
shapiro.test(NI$TRS)

#install.packages("bestNormalize")
require(bestNormalize)
TRSadj <- bestNormalize(NI$TRS, standardize = FALSE, allow_orderNorm = TRUE, out_of_sample = FALSE)
TRSadj$chosen_transform
shapiro.test(TRSadj$x.t)
NI$TRSadj <- TRSadj$x.t
head(NI)

#take out the two varieties that have not been tested in each location
#df_FB15M <- df_FB15M[!(df_FB15M$Variety %in% c("FLB15-0444", "FLB15-0652")), ]

#table(df_FB15M$Variety, df_FB15M$Loc)


library(lme4)
library(car)

# Define the mixed model formula
# Use () to nest Rep within Location, and Crop within PlantYear
# Use : to specify interactions between Cross, Location, and Crop
# Use | to specify that all factors are random effects
# Use (1 | ) to include an intercept
formula <- TRS ~ (1) + (1 | Cross) + (1 | Rep)  +  (1 | Crop)

# Fit the mixed model using the lmer function from the lme4 package
fit1 <- lmer(formula, data = NI, REML = TRUE)

# Print the model summary
summary(fit1)
anova(fit1)

# Extract the residuals
resid <- residuals(fit1)

# Test for normality of residuals
shapiro.test(resid)

# Test for homogeneity of variances by plotting
plot(fit1, resid = TRUE, sqrt(abs(resid(.))) ~ fitted(.))



# Trial level heritability
Var_G <- as.numeric(VarCorr(fit1)[1])  # genetic variance
Var_GE <- as.numeric(VarCorr(fit1)[2]) + as.numeric(VarCorr(fit1)[3]) + as.numeric(VarCorr(fit1)[4])# genotype-by-environment variance
Var_residual <- sigma(fit1)^2  # residual variance
e <- 3  # number of environments
r <- 2  # number of replicates per environment
H2_trial <- Var_G / (Var_G + Var_GE/e + Var_residual/e*r)

# Plot level heritability
n <- 2  # number of plots per genotype per environment
H2_plot <- Var_G / (Var_G + Var_residual/n)

# Print the results
cat("Trial level heritability: ", H2_trial, "\n")
cat("Plot level heritability: ", H2_plot, "\n")

#################################### MET analysis #########################
library(sommer)

blocks <- length(unique(NI$Rep))
crops <- length(unique(NI$Crop))

# Fitting genotype by environment models - with a common variance (diagnonal model); everything is non-related bc no GRM
TRSfitMETNI <- mmer(fixed = TRS ~ 1,
                  random = ~Cross + Crop + Rep,
                  rcov = ~ vsr(units),
                  data = NI)

summary(TRSfitMETNI)
# Broad-sense heritability
vpredict(TRSfitMETNI, h2 ~ V1 / ( V1 + V2/crops + V4/(crops*blocks) ) ) # trials level; replicates is not part of the heritability equation
vpredict(TRSfitMETNI, h2 ~ V1 / ( V1 + V2 + V4 ) ) # plot level; Vs are the order in the model


# Fitting genotype by environment models - unstructured model (US)
TRSfitMET.USNI <- mmer(fixed = TRS ~ 1,
                     random = ~Cross + Crop + Rep + vsr(usr(Cross), Crop),
                     rcov = ~ vsr(units),
                     data = NI)
summary(TRSfitMET.USNI)

# Broad-sense heritability
vpredict(TRSfitMET.USNI, h2 ~ V1 / ( V1 + V3/loc + V4/(loc*blocks) ) ) # trials level; replicates is not part of the heritability equation
vpredict(TRSfitMET.USNI, h2 ~ V1 / ( V1 + V3 + V4 ) ) # plot level; Vs are the order in the model


# Fitting genotype by environment models - unstructured model (US) + heterogeneous variance
TRSfitMET.US.HNI <- mmer(fixed = TRS ~ 1,
                       random = ~Cross +  Crop + Rep + vsr(usr(Cross), Crop),
                       rcov = ~ vsr(dsr(Cross), units),
                       data = NI)

summary(TRSfitMET.US.HNI)
# Broad-sense heritability
vpredict(TRSfitMET.US.HNI, h2 ~ V1 / ( V1 + V3/loc + V4/(loc*blocks) ) ) # trials level; replicates is not part of the heritability equation
vpredict(TRSfitMET.US.HNI, h2 ~ V1 / ( V1 + V3 + V4 ) ) # plot level; Vs are the order in the model

# comparing the models
anova(TRSfitMETNI, TRSfitMET.USNI)
anova(TRSfitMET.USNI, TRSfitMET.US.HNI)
#TRSfitMET.US.H seems to perform best



# predicting BLUPs 
BLUPS <- data.frame(
  MET = TRSfitMETNI$U$Cross$TRS)

head(BLUPS)
cor(BLUPS)

#modify the BLUPS dataframe
BLUPS$Cross <- rownames(BLUPS)
BLUPS$Cross <- sub("^Cross", "", BLUPS$Cross)
BLUPS <- dplyr::rename(BLUPS, BLUPNI = MET)
# Calculate mean and standard deviation of BLUPs
blup_mean <- mean(BLUPS$BLUPNI)
blup_sd <- sd(BLUPS$BLUPNI)
# Create a new column for highlighting
BLUPS$SignificantDiff <- ifelse(BLUPS$BLUPNI > (blup_mean + blup_sd), "Greater", ifelse(BLUPS$BLUPNI < (blup_mean - blup_sd), "Less", "Same"))
# Order the dataframe by BLUPs
BLUPS <- BLUPS[order(BLUPS$BLUPNI),]
# Create a bar plot with highlighted crosses
library(ggplot2)
ggplot(BLUPS, aes(x=reorder(Cross, BLUPNI), y=BLUPNI, fill=SignificantDiff)) +
  geom_bar(stat="identity") +
  geom_hline(yintercept = blup_mean, color="royalblue3", linetype="dashed", size=1) +
  labs(title="BLUPs of Crosses", x="Cross", y="BLUP") +
  scale_fill_manual(values=c("Greater"="seagreen3", "Less"="tomato3", "Same"="grey")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1))



# predicting the overal performance
predict.mmer(TRSfitMETNI, D="Cross")$pvals



library(dplyr)

# assuming your data frames are called "results" and "BLUPS"
merged_df <- merge(results, BLUPS, by = "Cross") %>%
  select(Cross, BLUPNI, BLUPALL, GenMean, SE_GenMean, Rank, Sens, SE_Sens, MSdeviation)

# assign the updated data frame back to the "results" variable
results <- merged_df


library(emmeans)
library(ggplot2)
# Compute pairwise comparisons using emmeans
emmeans(TRS_model, pairwise ~ Cross, adjust = "tukey")

# Calculate EMMs and CIs for each level of the "Variety" factor
emm <- emmeans(TRS_model, "Cross", type = "response", infer = c(TRUE, TRUE))

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
lsmeans <- lsmeans(TRS_model, "Cross")
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
pairs(lsmeans(TRS_model, "Cross"), adjust = "tukey")




install.packages("contrast")
library(tibble)
library(multcomp)
library(tidyverse)
library(broom)




boxplot(NI$TRS, col = "red")


# outlier detection and elimination; similar to first class
#outlier <- names(outlierTest(TRS_model)$p)
#TRS[outlier, "TRS"] <- NA


#For TRS
TRS_model <- lmer(TRS ~ Cross + Crop + (1|Rep), data = NI)
emm_table <- emmeans(TRS_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger")
# Get the summary statistics of the EMMs
TRS_avgs <- summary(emm_table, infer = c(TRUE, TRUE, TRUE, FALSE))

# Convert the TRS_avgs object to a data frame
TRS_avgs <- as.data.frame(TRS_avgs)



TRS_emm <- emmeans(TRS_model, ~ Cross, adjust = "sidak", type = "response", DDF = "Kenward-Roger", data = NI)
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

#Multi-Trait Selection using BLUPs 
BLUPS <- colnames(results[, c(2, 3)])

#correlation between BLUPS across all locstions and of NI
correlation <- cor(results$BLUPALL, results$BLUPNI, use = "pairwise.complete.obs", method = "pearson")

rank_BLUPALL <- rank(results$BLUPALL)
rank_BLUPNI <- rank(results$BLUPNI)
rank_cutoff_BLUPALL <- quantile(rank_BLUPALL, 0.85)
rank_cutoff_BLUPNI <- quantile(rank_BLUPNI, 0.85)
top_15_percent <- results[(rank_BLUPALL >= rank_cutoff_BLUPALL) & (rank_BLUPNI >= rank_cutoff_BLUPNI), ]
cat("Number of individuals in the top 15%:", nrow(top_15_percent))
top_varieties <- results$Cross %in% top_15_percent$Cross
point_colors <- ifelse(rank_BLUPALL >= rank_cutoff_BLUPALL & rank_BLUPNI >= rank_cutoff_BLUPNI, "seagreen3", "royalblue3")


ggplot(results, aes(x = BLUPALL, y = BLUPNI, label = Cross)) +
  geom_point(size = 3, alpha = 0.8, aes(color = point_colors)) +
  geom_smooth(method = "lm", se = FALSE, color = "#34495E") +
  scale_color_manual(name = "Quantile", values = c("royalblue3", "seagreen3"), labels = c("Bottom 85%", "Top 15%")) +
  labs(title = "Correlation between BLUPALL and BLUPNI",
  x = "BLUPALL",
  y = "BLUPNI") +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
  axis.title = element_text(size = 20),
  axis.text = element_text(size = 18),
  legend.position = c(0.85, 0.2),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank()) +
  geom_text(color = "black", check_overlap = TRUE, size = 3, nudge_x = 0.3, nudge_y = 0.3)

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
                     labels = c("Top 15%", "Bottom 85%", "Outliers"),
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

