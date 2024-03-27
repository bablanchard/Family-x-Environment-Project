library(readxl)
library("PerformanceAnalytics")
library(car)
library(dplyr)
library(nlme)
library(Matrix)
library(ggplot2)
library(lme4)
library(tidyr)


#read master
df_filtered <- read.csv("///Master_1.2.csv", header = TRUE)


#=============== phenotypic correlation between traits===============
sapply(df_filtered[, 11:22, 27], class)
df_filtered[, 11:22] <- lapply(df_filtered[, 11:22], as.numeric)
df_filtered[, 27] <- as.numeric(df_filtered[, 27])
round(cor(df_filtered[, c(11:22, 27)], use = "pairwise"), 2)
data_subset <- df_filtered[, c(11:18, 27)]
#install.packages("PerformanceAnalytics")
chart.Correlation(as.matrix(na.omit(data_subset)), histogram = TRUE, pch = 1)

#============make CY from PlotWeight and convert TRS to kg per Mg and make SY=========
# Convert plot area from square feet to hectares
plot_area_hectares <- 252 * 0.000009290304

# Convert weight from pounds to metric tons (megagrams)
df_filtered$PlotWeight_megagrams <- df_filtered$PlotWeight * 0.00045359237

# Convert weight to megagrams per hectare
df_filtered$CY <- df_filtered$PlotWeight_megagrams / plot_area_hectares

# Convert lbs to kg
df_filtered$TRS <- df_filtered$TRS * 453.592

# Convert tons to Megagrams
df_filtered$TRS <- df_filtered$TRS * (1 / 907.185)

# Rename column header if needed
colnames(df_filtered)[colnames(df_filtered) == "TRS"] <- "TRS_kg_per_Mg"

# Calculate SY (kilograms per hectare)
df_filtered$SY <- df_filtered$TRS_kg_per_Mg * df_filtered$CY

# Convert pounds to kilograms (1 pound = 0.453592 kilograms)
df_filtered$SW_kg <- df_filtered$SW * 0.453592
df_filtered$SW <- df_filtered$SW_kg
#=============Create Datasets and quickly view boxplots, histograms, and groups=======
SY <- df_filtered[complete.cases(df_filtered$SY) & df_filtered$SY != 0,]
# Subset the df_filtered dataframe to remove all observations with NA or 0 values in the TRS column
TRS_kg_per_Mg <- df_filtered[complete.cases(df_filtered$TRS_kg_per_Mg) & df_filtered$TRS_kg_per_Mg != 0,]
# Create a new dataframe without the unwanted rows with TRS
TRS_kg_per_Mg <- subset(TRS_kg_per_Mg, !(Cross %in% c("CP18-1611", "HB 17-3208")))
CY <- df_filtered[complete.cases(df_filtered$CY) & df_filtered$CY != 0,]
Stalk <- df_filtered[complete.cases(df_filtered$Stalk) & df_filtered$Stalk != 0,]
Height <- df_filtered[complete.cases(df_filtered$Height) & df_filtered$Height != 0,]
Dia <- df_filtered[complete.cases(df_filtered$Dia) & df_filtered$Dia != 0,]
PlotWeight <- df_filtered[complete.cases(df_filtered$PlotWeight) & df_filtered$PlotWeight != 0,]
Brix <- df_filtered[complete.cases(df_filtered$Brix) & df_filtered$Brix != 0,]
Fiber <- df_filtered[complete.cases(df_filtered$Fiber) & df_filtered$Fiber != 0,]
SW <- df_filtered[complete.cases(df_filtered$SW) & df_filtered$SW != 0,]
Sucrose <- df_filtered[complete.cases(df_filtered$Sucrose..) & df_filtered$Sucrose.. != 0,]

hist(TRS_kg_per_Mg$TRS_kg_per_Mg)
hist(Stalk$Stalk)
hist(Height$Height)
hist(Dia$Dia)
hist(PlotWeight$PlotWeight)
hist(Brix$Brix)
hist(Fiber$Fiber)
hist(SW$SW)

boxplot(TRS_kg_per_Mg$TRS_kg_per_Mg)
boxplot(Stalk$Stalk)
boxplot(Height$Height)
boxplot(Dia$Dia)
boxplot(PlotWeight$PlotWeight)
boxplot(Brix$Brix)
boxplot(Fiber$Fiber)
boxplot(SW$SW)


# Load the 'car' package and examine groups
table(TRS_kg_per_Mg$Cross)
table(TRS_kg_per_Mg$Location)
table(TRS_kg_per_Mg$Rep)
table(TRS_kg_per_Mg$Crop)
table(TRS_kg_per_Mg$Env)

table(TRS_kg_per_Mg$Cross, TRS_kg_per_Mg$Location)
table(TRS_kg_per_Mg$Crop, TRS_kg_per_Mg$Location, TRS_kg_per_Mg$PlantYear)
table(TRS_kg_per_Mg$Cross, TRS_kg_per_Mg$Location, TRS_kg_per_Mg$PlantYear)
table(TRS_kg_per_Mg$Cross, TRS_kg_per_Mg$Location, TRS_kg_per_Mg$PlantYear, TRS_kg_per_Mg$Crop)

#the harmonic mean of the number of observations csn be used to estimate heritability
# Calculate the number of observations per Cross
cross_counts <- table(TRS_kg_per_Mg$Cross)
# Filter out zero counts (if any)
non_zero_counts <- cross_counts[cross_counts > 0]
# Calculate the harmonic mean of the number of observations per Cross
TRSharmonic_mean <- 1 / mean(1 / non_zero_counts)
# Print the results
cat("Number of observations per Cross:\n")
print(cross_counts)
cat("\nHarmonic mean of the number of observations per Cross:", TRSharmonic_mean, "\n")


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


#===========TRS Sample Model and Outputs==============================
# 'TRS' is the response variable, and 'TRS' is the predictor variable
# Define the model with fixed and random effects
model <- lmer(TRS_kg_per_Mg ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                (1|Cross:Crop:Location) + (1|Cross:Crop:Location:PlantYear),
              data = TRS_kg_per_Mg)

# Fit the model
model_fit <- summary(model)

# Obtain variance components
var_components <- VarCorr(model)

# Obtain BLUPs (Best Linear Unbiased Predictors)
blups <- ranef(model)

# Print results
print(model_fit)
print(var_components)
print(blups)


library(ggplot2)

# Assuming var_components is a list with 16 components

# Extract variance components using VarCorr function
variances <- sapply(VarCorr(model), function(x) attr(x, "stddev")^2)

# Create a data frame for plotting
plot_data <- data.frame(Component = names(variances), Variance = variances)

# Create a stacked bar plot
ggplot(plot_data, aes(x = "", y = Variance, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(x = "", y = "Variance", fill = "Component") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))  # Reverse legend order for better readability


# Extract variance components using VarCorr function
var_components <- VarCorr(model)

# Find the variance component for the intercept of 'Cross' and 'Cross:PlantYear'
vc_cross <- as.numeric(var_components$Cross)
vc_cross_plantyear <- as.numeric(var_components$`Cross:PlantYear`)
vc_cross_Location <- as.numeric(var_components$`Cross:Location`)
vc_cross_Crop <- as.numeric(var_components$`Cross:Crop`)
vc_cross_plantyear_location <- as.numeric(var_components$`Cross:PlantYear:Location`)
vc_cross_crop_plantyear <- as.numeric(var_components$`Cross:Crop:PlantYear`)
vc_cross_crop_location <- as.numeric(var_components$`Cross:Crop:Location`)
vc_cross_crop_location_plantyear <- as.numeric(var_components$`Cross:Crop:Location:PlantYear`)
vc_residual <- as.numeric(sigma(model)^2)


# Calculate broad sense heritability
broad_heritability <- vc_cross / (vc_cross + vc_cross_plantyear/2 + vc_cross_Location/3 + vc_cross_Crop/4 +
                                    vc_cross_plantyear_location/5 + vc_cross_crop_plantyear/7 + vc_cross_crop_location/10 +
                                    vc_cross_crop_location_plantyear/15 + vc_residual/30)

# Print broad sense heritability
print(broad_heritability)






#=============Remake Empty Dataframes to collect BLUPs and relevant VCs=========
# Define datasets and their corresponding predicted variables
# Create empty data frames to store results
DataOutput <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))
DataVarComp <- data.frame(VC = c("Cross", "Cross:PlantYear", "Cross:Location", 
                                 "Cross:Crop", "Cross:PlantYear:Location", 
                                 "Cross:Crop:PlantYear", "Cross:Crop:Location", 
                                 "Cross:Crop:Location:PlantYear", "Residual", "Heritability"))
#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutput$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutput$Row <- c(1:50)
DataOutput$Cross <- DataOutput$Entry
DataOutput <- subset(DataOutput, select = -Entry)

#===========Define linear mixed models for each trait===========
#Define model with fixed and random effects
SYmodel <- lmer(SY ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                  (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                  (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                  (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                  (1|Cross:Crop:Location:PlantYear),
                data = SY)

SYvar_components <- VarCorr(SYmodel)
SYmodel_fit <- summary(SYmodel)

TRSmodel <- lmer(TRS_kg_per_Mg ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                   (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                   (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                   (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                   (1|Cross:Crop:Location:PlantYear),
                 data = TRS_kg_per_Mg)

TRSvar_components <- VarCorr(TRSmodel)
TRSmodel_fit <- summary(TRSmodel)


CYmodel <- lmer(CY ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                  (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                  (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                  (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                  (1|Cross:Crop:Location:PlantYear),
                data = CY)

CYvar_components <- VarCorr(CYmodel)
CYmodel_fit <- summary(CYmodel)


Fibermodel <- lmer(Fiber ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                     (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                     (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                     (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                     (1|Cross:Crop:Location:PlantYear),
                   data = Fiber)

Fibervar_components <- VarCorr(Fibermodel)
summary(Fibermodel)

Diamodel <- lmer(Dia ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                   (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                   (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                   (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                   (1|Cross:Crop:Location:PlantYear),
                 data = Dia)

Diavar_components <- VarCorr(Diamodel)
Diamodel_fit <- summary(Diamodel)


Stalkmodel <- lmer(Stalk ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                     (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                     (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                     (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                     (1|Cross:Crop:Location:PlantYear),
                   data = Stalk)

Stalkvar_components <- VarCorr(Stalkmodel)
summary(Stalkmodel)

Heightmodel <- lmer(Height ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                      (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                      (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                      (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                      (1|Cross:Crop:Location:PlantYear),
                    data = Height)

Heightvar_components <- VarCorr(Heightmodel)
summary(Heightmodel)


Brixmodel <- lmer(Brix ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                    (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                    (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                    (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                    (1|Cross:Crop:Location:PlantYear),
                  data = Brix)

Brixvar_components <- VarCorr(Brixmodel)
summary(Brixmodel)


SWmodel <- lmer(SW ~ 1 + (1|Cross) + (1|PlantYear) + (1|Location) + (1|Crop) + (1|PlantYear:Location) +
                  (1|Rep:(PlantYear:Location)) + (1|Crop:PlantYear) + (1|Crop:Location) + 
                  (1|Crop:PlantYear:Location) + (1|Crop:(Rep:(PlantYear:Location))) + (1|Cross:PlantYear) +
                  (1|Cross:Location) + (1|Cross:Crop) + (1|Cross:PlantYear:Location) + (1|Cross:Crop:PlantYear) + 
                  (1|Cross:Crop:Location:PlantYear),
                data = SW)

SWvar_components <- VarCorr(SWmodel)
summary(SWmodel)


#============Use LMMs to extract BLUPs for each Cross overall and for each trait==========
# List of model objects
models <- list(SYmodel, TRSmodel, CYmodel, Fibermodel, Diamodel, Stalkmodel, Heightmodel, Brixmodel, SWmodel)
model_names <- c("SY", "TRS", "CY", "Fiber", "Dia", "Stalk", "Height", "Brix", "SW")  # Names of the models



# Loop through each model and extract overall BLUPs for each trait
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract the predicted variable name from the model
  predicted_variable <- all.vars(formula(model))[1]
  
  # Extract BLUPs from the model
  blups <- ranef(model)$Cross
  
  # Loop through each unique value of Cross
  for (cross_value in unique(DataOutput$Cross)) {
    # Check if the cross_value exists in the BLUPs
    if (as.character(cross_value) %in% rownames(blups)) {
      # Extract BLUP for the current cross_value from the model
      blup <- blups[as.character(cross_value), "(Intercept)"]
      
      # Extract the predicted values from the model
      predicted_values <- predict(model, newdata = DataOutput, re.form = NA)
      
      # Calculate the overall mean of the predicted variable
      overall_mean <- mean(predicted_values, na.rm = TRUE)
      
      # Add the BLUP to the overall mean value of the predicted variable
      adjusted_blup <- blup + overall_mean
    } else {
      # If cross_value does not exist in the BLUPs, assign NA to BLUP
      adjusted_blup <- NA
    }
    
    # Add adjusted BLUP to new column in DataOutput
    DataOutput[[paste0(model_name, "_BLUP")]][DataOutput$Cross == cross_value] <- adjusted_blup
  }
}
write.csv(DataOutput, file = "DataOutputOverallBlups.csv", row.names = FALSE)


#==============Extract BLUPs per location===========
# Define datasets and their corresponding predicted variables
# Create empty data frames to store results
DataOutputPerLoc <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutputPerLoc$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutputPerLoc$Row <- c(1:50)
DataOutputPerLoc$Cross <- DataOutputPerLoc$Entry
DataOutputPerLoc <- subset(DataOutputPerLoc, select = -Entry)


# Loop through each model and extract overall BLUPs for each trait
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract the predicted variable name from the model
  predicted_variable <- all.vars(formula(model))[1]
  
  # Extract BLUPs from the model
  blups <- ranef(model)$'Cross:Location'
  # Convert rownames to a new column and split into two columns
  blups$Cross_Location <- rownames(blups)
  blups <- separate(blups, Cross_Location, into = c("Cross", "Location"), sep = ":", remove = TRUE)
  # Pivot wider to create separate columns for each unique value in 'Location'
  blups <- pivot_wider(blups, names_from = Location, values_from = "(Intercept)")
  
  
  
  # Loop through each unique value of Cross
  for (cross_value in unique(DataOutputPerLoc$Cross)) {
    # Check if the cross_value exists in the BLUPs
    if (as.character(cross_value) %in% blups$Cross) {
      # Extract BLUPs for the current cross_value from the model
      blup_ni <- blups[blups$Cross == as.character(cross_value), "NI"]
      blup_nr <- blups[blups$Cross == as.character(cross_value), "NR"]
      blup_sg <- blups[blups$Cross == as.character(cross_value), "SG"]
      
      # Extract the predicted values from the model
      predicted_values <- predict(model, newdata = DataOutputPerLoc, re.form = NA)
      # Calculate the overall mean of the predicted variable
      overall_mean <- mean(predicted_values, na.rm = TRUE)
      
      # Calculate the adjusted BLUPs by adding the overall mean to each BLUP
      adjusted_blup_ni <- blup_ni + overall_mean
      adjusted_blup_nr <- blup_nr + overall_mean
      adjusted_blup_sg <- blup_sg + overall_mean
    } else {
      # If cross_value does not exist in the BLUPs, assign NA to BLUP
      adjusted_blup_ni <- NA
      adjusted_blup_nr <- NA
      adjusted_blup_sg <- NA
    }
    
    # Add adjusted BLUP to new column in DataOutput
    DataOutputPerLoc[[paste0(model_name, "_BLUPNI")]][DataOutputPerLoc$Cross == cross_value] <- adjusted_blup_ni
    DataOutputPerLoc[[paste0(model_name, "_BLUPNR")]][DataOutputPerLoc$Cross == cross_value] <- adjusted_blup_nr
    DataOutputPerLoc[[paste0(model_name, "_BLUPSG")]][DataOutputPerLoc$Cross == cross_value] <- adjusted_blup_sg
    
  }
}


# Get column names starting from the third column
cols_to_extract <- colnames(DataOutputPerLoc)[3:ncol(DataOutputPerLoc)]

# Initialize an empty list to store column lists
column_lists <- list()

# Loop through each column to extract lists
for (col_name in cols_to_extract) {
  # Check if the column is a list
  if (is.list(DataOutputPerLoc[[col_name]])) {
    # Extract the list column and store it in the list
    column_lists[[col_name]] <- DataOutputPerLoc[[col_name]]
  }
}

# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each list in column_lists
for (i in seq_along(column_lists)) {
  # Create a data frame from the current list
  df <- data.frame(column = unlist(column_lists[[i]]))
  # Add the data frame to the list
  data_frames[[i]] <- df
}

# Combine the data frames into a single data frame
combined_df <- do.call(cbind, data_frames)

# Optionally, you can rename the columns to match the names of the original lists
colnames(combined_df) <- names(column_lists)

# Extract DataOutputPerLoc$Cross
Cross <- DataOutputPerLoc$Cross

# Combine DataOutputPerLoc$Cross with combined_data
combined_df <- cbind(Cross, combined_df)
DataOutputPerLoc <- combined_df



write.csv(DataOutputPerLoc, file = "DataOutputLocBlups.csv", row.names = F)


#============examine family ranking per location=========
# Select relevant columns (performance metrics)
selected_columns <- c("Cross", "SY_BLUPNI", "SY_BLUPNR", "SY_BLUPSG")

# Prepare data in long format for plotting
long_data <- DataOutputPerLoc %>%
  select(all_of(selected_columns)) %>%
  pivot_longer(cols = -Cross, names_to = "PerformanceMetric", values_to = "Value")
long_data$PerformanceMetric <- gsub("SY_", "", long_data$PerformanceMetric)

parallel_plot <- ggplot(long_data, aes(x = PerformanceMetric, y = Value, color = Cross, group = Cross)) +
  geom_line(show.legend = FALSE) +  # Remove legend from lines
  geom_point(show.legend = FALSE) +  # Remove legend from points
  labs(title = ,
       x = "Locations", y = "Sugar Yield (kg/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust size of x-axis labels
        axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels
        axis.title = element_text(size = 14, face = "bold"),  # Adjust size and style of axis titles
        plot.title = element_text(size = 16, face = "bold"))  # Adjust size and style of plot title

# Print the modified plot
print(parallel_plot)




data_subset <- DataOutputPerLoc[, c(2:4)]

library(psych)

pairs.panels(data_subset,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals

#for TRS
# Select relevant columns (performance metrics)
selected_columns <- c("Cross", "TRS_BLUPNI", "TRS_BLUPNR", "TRS_BLUPSG")

# Prepare data in long format for plotting
long_data <- DataOutputPerLoc %>%
  select(all_of(selected_columns)) %>%
  pivot_longer(cols = -Cross, names_to = "PerformanceMetric", values_to = "Value")
long_data$PerformanceMetric <- gsub("TRS_", "", long_data$PerformanceMetric)

parallel_plot <- ggplot(long_data, aes(x = PerformanceMetric, y = Value, color = Cross, group = Cross)) +
  geom_line(show.legend = FALSE) +  # Remove legend from lines
  geom_point(show.legend = FALSE) +  # Remove legend from points
  labs(title = "Comparison of Performance Metrics by Family (Cross)",
       x = "Locations", y = "TRS (kg/Mg)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust size of x-axis labels
        axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels
        axis.title = element_text(size = 14, face = "bold"),  # Adjust size and style of axis titles
        plot.title = element_text(size = 16, face = "bold"))  # Adjust size and style of plot title

# Print the modified plot
print(parallel_plot)




data_subset <- DataOutputPerLoc[, c(5:7)]
pairs.panels(data_subset,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals

#for CY
# Select relevant columns (performance metrics)
selected_columns <- c("Cross", "CY_BLUPNI", "CY_BLUPNR", "CY_BLUPSG")

# Prepare data in long format for plotting
long_data <- DataOutputPerLoc %>%
  select(all_of(selected_columns)) %>%
  pivot_longer(cols = -Cross, names_to = "PerformanceMetric", values_to = "Value")
long_data$PerformanceMetric <- gsub("CY_", "", long_data$PerformanceMetric)

parallel_plot <- ggplot(long_data, aes(x = PerformanceMetric, y = Value, color = Cross, group = Cross)) +
  geom_line(show.legend = FALSE) +  # Remove legend from lines
  geom_point(show.legend = FALSE) +  # Remove legend from points
  labs(title = "",
       x = "Locations", y = "CY (Mg/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust size of x-axis labels
        axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels
        axis.title = element_text(size = 14, face = "bold"),  # Adjust size and style of axis titles
        plot.title = element_text(size = 16, face = "bold"))  # Adjust size and style of plot title

# Print the modified plot
print(parallel_plot)




data_subset <- DataOutputPerLoc[, c(8:10)]
pairs.panels(data_subset,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals


#==============extract blups per crop=========
# Define datasets and their corresponding predicted variables
# Create empty data frames to store results
DataOutputPerCrop <- data.frame(matrix(vector(),50,1, dimnames=list(c(), c("Entry"))))

#fill empty dataframe with 1-300 so that the cbind will work later on
DataOutputPerCrop$Entry <- unique(df_filtered[,2]) #fill in Entry numbers
DataOutputPerCrop$Row <- c(1:50)
DataOutputPerCrop$Cross <- DataOutputPerCrop$Entry
DataOutputPerCrop <- subset(DataOutputPerCrop, select = -Entry)


# Loop through each model and extract overall BLUPs for each trait
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract the predicted variable name from the model
  predicted_variable <- all.vars(formula(model))[1]
  
  # Extract BLUPs from the model
  blups <- ranef(model)$'Cross:Crop'
  # Convert rownames to a new column and split into two columns
  blups$Cross_Crop <- rownames(blups)
  blups <- separate(blups, Cross_Crop, into = c("Cross", "Crop"), sep = ":", remove = TRUE)
  # Pivot wider to create separate columns for each unique value in 'Location'
  blups <- pivot_wider(blups, names_from = Crop, values_from = "(Intercept)")
  
  
  
  # Loop through each unique value of Cross
  for (cross_value in unique(DataOutputPerCrop$Cross)) {
    # Check if the cross_value exists in the BLUPs
    if (as.character(cross_value) %in% blups$Cross) {
      # Extract BLUPs for the current cross_value from the model
      blup_pc <- blups[blups$Cross == as.character(cross_value), "0"]
      blup_1r <- blups[blups$Cross == as.character(cross_value), "1"]
      blup_2r <- blups[blups$Cross == as.character(cross_value), "2"]
      blup_3r <- blups[blups$Cross == as.character(cross_value), "3"]
      
      
      # Extract the predicted values from the model
      predicted_values <- predict(model, newdata = DataOutputPerCrop, re.form = NA)
      # Calculate the overall mean of the predicted variable
      overall_mean <- mean(predicted_values, na.rm = TRUE)
      
      # Calculate the adjusted BLUPs by adding the overall mean to each BLUP
      adjusted_blup_pc <- blup_pc + overall_mean
      adjusted_blup_1r <- blup_1r + overall_mean
      adjusted_blup_2r <- blup_2r + overall_mean
      adjusted_blup_3r <- blup_3r + overall_mean
      
    } else {
      # If cross_value does not exist in the BLUPs, assign NA to BLUP
      adjusted_blup_pc <- NA
      adjusted_blup_1r <- NA
      adjusted_blup_2r <- NA
      adjusted_blup_3r <- NA
    }
    
    # Add adjusted BLUP to new column in DataOutput
    DataOutputPerCrop[[paste0(model_name, "_BLUPPC")]][DataOutputPerCrop$Cross == cross_value] <- adjusted_blup_pc
    DataOutputPerCrop[[paste0(model_name, "_BLUP1R")]][DataOutputPerCrop$Cross == cross_value] <- adjusted_blup_1r
    DataOutputPerCrop[[paste0(model_name, "_BLUP2R")]][DataOutputPerCrop$Cross == cross_value] <- adjusted_blup_2r
    DataOutputPerCrop[[paste0(model_name, "_BLUP3R")]][DataOutputPerCrop$Cross == cross_value] <- adjusted_blup_3r
    
  }
}


# Get column names starting from the third column
cols_to_extract <- colnames(DataOutputPerCrop)[3:ncol(DataOutputPerCrop)]

# Initialize an empty list to store column lists
column_lists <- list()

# Loop through each column to extract lists
for (col_name in cols_to_extract) {
  # Check if the column is a list
  if (is.list(DataOutputPerCrop[[col_name]])) {
    # Extract the list column and store it in the list
    column_lists[[col_name]] <- DataOutputPerCrop[[col_name]]
  }
}


# Initialize an empty list to store data frames
data_frames <- list()

# Loop through each list in column_lists
for (i in seq_along(column_lists)) {
  # Create a data frame from the current list
  df <- data.frame(column = unlist(column_lists[[i]]))
  # Add the data frame to the list
  data_frames[[i]] <- df
}

# Combine the data frames into a single data frame
combined_df <- do.call(cbind, data_frames)

# Optionally, you can rename the columns to match the names of the original lists
colnames(combined_df) <- names(column_lists)

# Extract DataOutputPerLoc$Cross
Cross <- DataOutputPerCrop$Cross

# Combine DataOutputPerLoc$Cross with combined_data
combined_df <- cbind(Cross, combined_df)
DataOutputPerCrop <- combined_df

write.csv(DataOutputPerCrop, file = "DataOutputCropBlups.csv", row.names = F)
#============examine family ranking per crop=========
library(ggplot2)

# Sample data (replace with your actual data)
# DataOutputPerLoc <- your_data_frame

# Select relevant columns (performance metrics)
selected_columns <- c("Cross", "SY_BLUPPC", "SY_BLUP1R", "SY_BLUP2R", "SY_BLUP3R")

# Prepare data in long format for plotting
long_data <- DataOutputPerCrop %>%
  select(all_of(selected_columns)) %>%
  pivot_longer(cols = -Cross, names_to = "PerformanceMetric", values_to = "Value")
long_data$PerformanceMetric <- gsub("SY_", "", long_data$PerformanceMetric)
# Define the desired order of performance metrics
desired_order <- c("BLUPPC", "BLUP1R", "BLUP2R", "BLUP3R")

# Convert PerformanceMetric to a factor with the desired order
long_data$PerformanceMetric <- factor(long_data$PerformanceMetric, levels = desired_order)

# Sort long_data by PerformanceMetric
long_data <- long_data[order(long_data$PerformanceMetric), ]

parallel_plot <- ggplot(long_data, aes(x = PerformanceMetric, y = Value, color = Cross, group = Cross)) +
  geom_line(show.legend = FALSE) +  # Remove legend from lines
  geom_point(show.legend = FALSE) +  # Remove legend from points
  labs(title = "Comparison of Performance Metrics by Family (Cross)",
       x = "Ratoon Crops", y = "Sugar Yield (kg/ha)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Adjust size of x-axis labels
        axis.text.y = element_text(size = 12),  # Adjust size of y-axis labels
        axis.title = element_text(size = 14, face = "bold"),  # Adjust size and style of axis titles
        plot.title = element_text(size = 16, face = "bold"))  # Adjust size and style of plot title

# Print the modified plot
print(parallel_plot)




data_subset <- DataOutputPerCrop[, c(2:5)]

library(psych)

pairs.panels(data_subset,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = FALSE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = TRUE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,           # pch symbol
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             hist.col = 4,       # Histograms color
             stars = TRUE,       # If TRUE, adds significance level with stars
             ci = TRUE)          # If TRUE, adds confidence intervals
#==============Use LMMs to extract VCs from each model and input into DataVarComp dataframe=======
# Loop through each model and extract relevant VCs for heritability
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract variance components using VarCorr function
  var_components <- VarCorr(model)
  
  # Find the variance component for each model
  for (j in seq_len(nrow(DataVarComp))) {
    vc_name <- as.character(DataVarComp$VC[j])
    
    # Check if the variance component exists in var_components
    if (vc_name %in% names(var_components)) {
      # Extract the variance component value
      vc_value <- var_components[[vc_name]]
    } else {
      # If variance component does not exist, assign NA
      vc_value <- NA
    }
    
    # Create a new column in DataVarComp with model_name as the column name
    DataVarComp[[model_name]][j] <- vc_value
  }
}


# Loop through each model and calculate residuals
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract variance components using VarCorr function
  var_components <- VarCorr(model)
  
  # Find the residual error for each model
  residual_error <- sigma(model)^2
  
  # Fill in the residual error for the corresponding model in DataVarComp
  DataVarComp[[model_name]][DataVarComp$VC == "Residual"] <- residual_error
}

# Loop through each row of the dataframe to change NAs to 0s
for (row_index in 1:nrow(DataVarComp)) {
  # Check if the current row is the row where VC is Heritability
  if (DataVarComp$VC[row_index] == "Heritability") {
    next  # Skip to the next iteration if the row is for Heritability
  } else {
    # Replace NA values with 0 in the current row
    DataVarComp[row_index, is.na(DataVarComp[row_index, ])] <- 0
  }
}

#extract values from the dataframe to calculate the broad sense heritability for each trait
for (col_index in 2:ncol(DataVarComp)) {
  if (is.na(DataVarComp[nrow(DataVarComp), col_index])) {
    # Extract values from the first and ninth rows of the current column
    value_row1 <- as.numeric(DataVarComp[1, col_index])
    value_row2 <- as.numeric(DataVarComp[2, col_index])
    value_row3 <- as.numeric(DataVarComp[3, col_index])
    value_row4 <- as.numeric(DataVarComp[4, col_index])
    value_row5 <- as.numeric(DataVarComp[5, col_index])
    value_row6 <- as.numeric(DataVarComp[6, col_index])
    value_row7 <- as.numeric(DataVarComp[7, col_index])
    value_row8 <- as.numeric(DataVarComp[8, col_index])
    value_row9 <- as.numeric(DataVarComp[9, col_index])
    
    # Print the extracted values for debugging
    cat("Values extracted for calculation in column", col_index, ":\n")
    cat("Row 1:", value_row1, "\n")
    cat("Row 2:", value_row2, "\n")
    cat("Row 3:", value_row3, "\n")
    cat("Row 4:", value_row4, "\n")
    cat("Row 5:", value_row5, "\n")
    cat("Row 6:", value_row6, "\n")
    cat("Row 7:", value_row7, "\n")
    cat("Row 8:", value_row8, "\n")
    cat("Row 9:", value_row9, "\n")
    
    # Calculate heritability by dividing the value in Row1 by the value in Row9
    heritability <- value_row1 / (value_row1 + value_row2/2 + value_row3/3 + value_row4/4 +
                                    value_row5/5 + value_row6/7 + value_row7/10 + value_row8/15 +
                                    value_row9/30)
    
    # Print the calculated heritability for debugging
    cat("Calculated heritability in column", col_index, ":", heritability, "\n")
    
    # Fill in the calculated heritability value in the corresponding cell
    DataVarComp[nrow(DataVarComp), col_index] <- heritability
    
    # Print the updated value in the corresponding cell for debugging
    cat("Updated value in column", col_index, ":", DataVarComp[nrow(DataVarComp), col_index], "\n\n")
  }
}

# Define a function to convert numbers from scientific notation to decimal notation with 9 decimal points
format_decimal <- function(x) {
  format(as.numeric(x), scientific = FALSE, digits = 9)
}

# Apply the format_decimal function to all elements in the columns after the first one
DataVarComp[, -1] <- lapply(DataVarComp[, -1], format_decimal)

#=============create complete dataset that contains all VCs with traits melted in a trait column==========
# Initialize the combined_var_components dataframe
combined_var_components <- data.frame(Model = character(),
                                      VC = character(),
                                      Value = numeric(),
                                      SD = numeric(),
                                      stringsAsFactors = FALSE)

# Loop through each model
for (i in seq_along(models)) {
  model <- models[[i]]
  model_name <- model_names[i]
  
  # Extract variance components using VarCorr function
  var_components <- VarCorr(model)
  
  # Extract variance components manually and create a dataframe
  var_df <- data.frame(
    Model = model_name,
    VC = names(var_components),
    Value = unlist(var_components),
    SD = sapply(var_components, function(x) sqrt(x[1, 1]))
  )
  
  # Combine variance components for all models
  combined_var_components <- rbind(combined_var_components, var_df)
}

# Calculate LowerBoundCI
combined_var_components$LowerBoundCI <- combined_var_components$Value - (1.96 * combined_var_components$SD)
combined_var_components$UpperBoundCI <- combined_var_components$Value + (1.96 * combined_var_components$SD)
# Create a new column 'Significant' in combined_var_components dataframe
combined_var_components$Significant <- ifelse(combined_var_components$LowerBoundCI > 0, "Yes", "No")


# Convert VC column to factor to ensure correct order in the plot
combined_var_components$VC <- factor(combined_var_components$VC, levels = unique(combined_var_components$VC))

#==========create various stacked bar plots to display the portions of VCs in each trait model=======
# Create stacked bar plot to display all model VCs
ggplot(combined_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Variance Components of All Models", x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sort the dataframe by Model and VC columns
combined_var_components <- combined_var_components[order(combined_var_components$Model, combined_var_components$VC),]

# Perform scaling transformation within each Model group
combined_var_components <- within(combined_var_components, {
  scaled_Value <- ave(Value, Model, FUN = function(x) scale(x, center = TRUE, scale = TRUE))
})

# Plot the scaled values with modified aesthetics
ggplot(combined_var_components, aes(x = Model, y = scaled_Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = bquote(bold("Variance Components of All Models")),
       x = bquote(bold("Model")), 
       y = bquote(bold("Scaled Value")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))



# Filter combined_var_components to include only TRS, CY, and Height models
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("CY", "Height", "Dia"), ]

# Create stacked bar plot
ggplot(filtered_var_components, aes(x = Model, y = scaled_Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = bquote(bold("Variance Components of CY, Height, and Dia")),
       x = bquote(bold("Model")), 
       y = bquote(bold("Scaled Value")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

# Create stacked bar plot with labels within the stacks
ggplot(filtered_var_components, aes(x = Model, y = scaled_Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = "Variance Components of TRS and CY", x = "Model", y = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#For just TRS
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("TRS"), ]

ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_manual(values = ifelse(filtered_var_components$VC %in% c(""), "black", scales::hue_pal()(length(unique(filtered_var_components$VC))))) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("Variance")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just Brix
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("Brix"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just Fiber
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("Fiber"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just CY
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("CY"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("Variance")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just Stalks
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("Stalk"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just Height
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("Height"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  scale_fill_manual(values = ifelse(filtered_var_components$VC %in% c(""), "black", scales::hue_pal()(length(unique(filtered_var_components$VC))))) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just Dia
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("Dia"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "none",  # Remove the legend
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

#For just SW
filtered_var_components <- combined_var_components[combined_var_components$Model %in% c("SW"), ]
ggplot(filtered_var_components, aes(x = Model, y = Value, fill = VC)) +
  geom_bar(stat = "identity", position = "stack", color = "white") +
  geom_text(aes(label = VC), position = position_stack(vjust = 0.5), color = "white", size = 3) +
  labs(title = bquote(bold("")),
       x = bquote(bold("")), 
       y = bquote(bold("")),
       fill = bquote(bold("Variance Components"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))



#==========using the DataVarComp dataframe to calculate GCV, PCV, CV=========
# Create a new row with the desired values
new_row <- c("OverallMean", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
# Bind the new row to the dataframe
DataVarComp <- rbind(DataVarComp, new_row)
# Iterate through each model and its corresponding name to get overall mean
for (i in seq_along(models)) {
  # Extract the model name and object
  model_name <- model_names[i]
  model <- models[[i]]
  
  # Extract overall mean for the trait from the model
  overall_mean <- fixef(model)["(Intercept)"]
  
  # Update DataVarComp with the overall mean for the corresponding trait
  DataVarComp[DataVarComp$VC == "OverallMean", model_name] <- overall_mean
}
new_row <- c("GCV", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
DataVarComp <- rbind(DataVarComp, new_row)
#extract values from the dataframe to calculate the genetic coefficient of Variation
for (col_index in 2:ncol(DataVarComp)) {
  if (is.na(DataVarComp[nrow(DataVarComp), col_index])) {
    # Extract values from the first and ninth rows of the current column
    value_row1 <- as.numeric(DataVarComp[1, col_index])
    value_row11 <- as.numeric(DataVarComp[11, col_index])
    
    
    # Print the extracted values for debugging
    cat("Values extracted for calculation in column", col_index, ":\n")
    cat("Row 1:", value_row1, "\n")
    cat("Row 11:", value_row2, "\n")
    
    
    # Calculate GCV by dividing the value in Row1 by the value in Row9
    GCV <- ((sqrt(value_row1)) / (value_row11)) * 100
    
    # Print the calculated GCV for debugging
    cat("Calculated GCV in column", col_index, ":", GCV, "\n")
    
    # Fill in the calculated GCV value in the corresponding cell
    DataVarComp[nrow(DataVarComp), col_index] <- GCV
    
    # Print the updated value in the corresponding cell for debugging
    cat("Updated value in column", col_index, ":", DataVarComp[nrow(DataVarComp), col_index], "\n\n")
  }
}
new_row <- c("PCV", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
DataVarComp <- rbind(DataVarComp, new_row)
#extract values from the dataframe to calculate the phenotypic coefficient of variation
for (col_index in 2:ncol(DataVarComp)) {
  if (is.na(DataVarComp[nrow(DataVarComp), col_index])) {
    # Extract values from the first and ninth rows of the current column
    value_row1 <- as.numeric(DataVarComp[1, col_index])
    value_row2 <- as.numeric(DataVarComp[2, col_index])
    value_row3 <- as.numeric(DataVarComp[3, col_index])
    value_row4 <- as.numeric(DataVarComp[4, col_index])
    value_row5 <- as.numeric(DataVarComp[5, col_index])
    value_row6 <- as.numeric(DataVarComp[6, col_index])
    value_row7 <- as.numeric(DataVarComp[7, col_index])
    value_row8 <- as.numeric(DataVarComp[8, col_index])
    value_row9 <- as.numeric(DataVarComp[9, col_index])
    value_row11 <- as.numeric(DataVarComp[11, col_index])
    
    # Print the extracted values for debugging
    cat("Values extracted for calculation in column", col_index, ":\n")
    cat("Row 1:", value_row1, "\n")
    cat("Row 2:", value_row2, "\n")
    cat("Row 3:", value_row3, "\n")
    cat("Row 4:", value_row4, "\n")
    cat("Row 5:", value_row5, "\n")
    cat("Row 6:", value_row6, "\n")
    cat("Row 7:", value_row7, "\n")
    cat("Row 8:", value_row8, "\n")
    cat("Row 9:", value_row9, "\n")
    cat("Row 11:", value_row11, "\n")
    
    
    # Calculate PCV by dividing the value in Row1 by the value in Row9
    PCV <- ((sqrt(value_row1 + value_row2/2 + value_row3/3 + value_row4/4 +
                    value_row5/5 + value_row6/7 + value_row7/10 + value_row8/15 +
                    value_row9/30)) / (value_row11)) * 100
    
    # Print the calculated PCV for debugging
    cat("Calculated PCV in column", col_index, ":", PCV, "\n")
    
    # Fill in the calculated PCV value in the corresponding cell
    DataVarComp[nrow(DataVarComp), col_index] <- PCV
    
    # Print the updated value in the corresponding cell for debugging
    cat("Updated value in column", col_index, ":", DataVarComp[nrow(DataVarComp), col_index], "\n\n")
  }
}

new_row <- c("SD", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
DataVarComp <- rbind(DataVarComp, new_row)
# Loop through each model and collect the standard deviation in a row in DataVarComp
for (i in seq_along(models)) {
  # Extract the residuals from the model
  response_variable <- residuals(models[[i]])
  
  # Calculate the standard deviation of the response variable
  model_sd <- sd(response_variable, na.rm = TRUE)
  
  # Get the model name
  model_name <- model_names[i]
  
  # Find the column index in DataVarComp dataframe corresponding to the model name
  col_index <- which(names(DataVarComp) == model_name)
  
  # If the column index is found, update the cell in the "SD" row of DataVarComp
  if (length(col_index) > 0) {
    DataVarComp[DataVarComp$VC == "SD", col_index] <- model_sd
  } else {
    # Print a warning if the column index is not found
    warning(sprintf("Column index not found for model %s", model_name))
  }
}

new_row <- c("CV", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
DataVarComp <- rbind(DataVarComp, new_row)
# Find the row index where DataVarComp$VC == "CV"
cv_row_index <- which(DataVarComp$VC == "CV")
# Get the column names where the values are NA in the "CV" row
cv_columns <- colnames(DataVarComp)[is.na(DataVarComp[cv_row_index,])]

# Convert columns in cv_columns to numeric
DataVarComp[, cv_columns] <- apply(DataVarComp[, cv_columns], 2, as.numeric)

# Loop through each column where CV needs to be calculated
for (cv_column in cv_columns) {
  # Get the column index for SD and OverallMean
  col_sd <- which(DataVarComp$VC == "SD")
  col_mean <- which(DataVarComp$VC == "OverallMean")
  
  # Calculate the coefficient of variation (CV)
  cv <- (DataVarComp[col_sd, cv_column] / DataVarComp[col_mean, cv_column]) * 100
  
  # Update the corresponding cell in the "CV" row of DataVarComp
  DataVarComp[cv_row_index, cv_column] <- cv
}

# Define a function to convert numbers from scientific notation to decimal notation with 9 decimal points
format_decimal <- function(x) {
  format(as.numeric(x), scientific = FALSE, digits = 9)
}

# Apply the format_decimal function to all elements in the columns after the first one
DataVarComp[, -1] <- lapply(DataVarComp[, -1], format_decimal)
#============make a visual of the many diagnostics created with DataVarComp=============
# Filter DataVarComp for rows where VC is Heritability, GCV, PCV, or CV
filtered_data <- DataVarComp[DataVarComp$VC %in% c("Heritability", "GCV", "PCV", "CV"), ]


# Transpose the filtered_data
flipped_data <- t(filtered_data[, -1])  # Exclude the first column (VC) before transposing
colnames(flipped_data)[colnames(flipped_data) == "10"] <- "Heritability"
colnames(flipped_data)[colnames(flipped_data) == "12"] <- "GCV"
colnames(flipped_data)[colnames(flipped_data) == "13"] <- "PCV"
colnames(flipped_data)[colnames(flipped_data) == "15"] <- "CV"


str(flipped_data)

# Convert to dataframe
flipped_data <- as.data.frame(flipped_data, stringsAsFactors = FALSE)

# Change all columns to numeric
flipped_data[] <- lapply(flipped_data, as.numeric)
# Create the bar plot using ggplot2
# Create the bar plot using ggplot2
ggplot(data = NULL, aes(x = rownames(flipped_data), y = flipped_data$PCV, fill = rownames(flipped_data))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(flipped_data$PCV, 2)), vjust = -0.5) +  # Add text labels with values rounded to 2 decimal places
  labs(title = "Phenotypic Coefficient of Variation",
       x = "Traits",
       y = "PCV") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = FALSE)  # Remove the legend


#=========use DataVarComp to simulate genetic gains=========
new_row <- c("GG", NA, NA, NA, NA, NA, NA, NA, NA, NA)  # Adjust for the number of columns in your dataframe
DataVarComp <- rbind(DataVarComp, new_row)
#extract values from the dataframe to calculate the phenotypic coefficient of variation
for (col_index in 2:ncol(DataVarComp)) {
  if (is.na(DataVarComp[nrow(DataVarComp), col_index])) {
    # Extract values from the first and ninth rows of the current column
    H <- as.numeric(DataVarComp[10, col_index])
    GCV <- as.numeric(DataVarComp[12, col_index])
    
    # Print the extracted values for debugging
    cat("Values extracted for calculation in column", col_index, ":\n")
    cat("Row 10:", H, "\n")
    cat("Row 12:", GCV, "\n")
    
    
    # Calculate PCV by dividing the value in Row1 by the value in Row9
    GG <- (2.3695 * sqrt(H) * GCV)
    
    # Print the calculated PCV for debugging
    cat("Calculated GG in column", col_index, ":", GG, "\n")
    
    # Fill in the calculated PCV value in the corresponding cell
    DataVarComp[nrow(DataVarComp), col_index] <- GG
    
    # Print the updated value in the corresponding cell for debugging
    cat("Updated value in column", col_index, ":", DataVarComp[nrow(DataVarComp), col_index], "\n\n")
  }
}

# Filter DataVarComp for rows where VC is Heritability, GCV, PCV, or CV
filtered_data <- DataVarComp[DataVarComp$VC %in% c("GG"), ]


# Transpose the filtered_data
flipped_data <- t(filtered_data[, -1])  # Exclude the first column (VC) before transposing
colnames(flipped_data)[colnames(flipped_data) == "16"] <- "GG"


str(flipped_data)

# Convert to dataframe
flipped_data <- as.data.frame(flipped_data, stringsAsFactors = FALSE)

# Change all columns to numeric
flipped_data[] <- lapply(flipped_data, as.numeric)
# Create the bar plot using ggplot2
# Create the bar plot using ggplot2
ggplot(data = NULL, aes(x = rownames(flipped_data), y = flipped_data$GG, fill = rownames(flipped_data))) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(flipped_data$GG, 5)), vjust = -0.5) +  # Add text labels with values rounded to 2 decimal places
  labs(title = "Genetic Gain per cycle in standard deviations from the population mean",
       x = "Traits",
       y = "Standard Deviations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = FALSE)  # Remove the legend

#================GGE Biplots==============
####################### Stability and adaptability - Finlay-Wilkinson #########################
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)

library(statgenGxE)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = SY, genotype = "Cross", trial = "Location", year = "Crop")
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "SY")
summary(dropsFW)

# let's take a look at the output
names(dropsFW)
dropsFW$estimates
dropsFW$anova
dropsFW$envEffs
dropsFW$TD
dropsFW$fittedGeno
dropsFW$trait
dropsFW$nGeno
dropsFW$nEnv
dropsFW$tol
dropsFW$iter

## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")
## Plot in  descending order.
plot(dropsTD, plotType = "box", traits = "SY", colorTrialBy = "trial",
     orderBy = "descending")
## Color the histograms for trials based on the variable scenarioFull.
plot(dropsTD, plotType = "scatter", traits = "SY", 
     colorTrialBy = "trial", 
     trialOrder = c("NR", "NI", "SG"))
## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "SY",
                           nestingFactor = "year")
summary(dropsVarComp)
vc(dropsVarComp)
herit(dropsVarComp)
plot(dropsVarComp)

#================GGE-Biplot Analysis ========================
library(metan)
# Assuming SY is your data frame and Cross is the column you want to modify
SYGGE <- SY
SYGGE$Cross <- sub(".*-", "", SYGGE$Cross)
TRSGGE <- TRS_kg_per_Mg
TRSGGE$Cross <- sub(".*-", "", TRSGGE$Cross)
CYGGE <- CY
CYGGE$Cross <- sub(".*-", "", CYGGE$Cross)

model.ggeCY <- gge(CYGGE, Location, Cross, CY, svp = "symmetrical")
model.ggeTRS <- gge(TRSGGE, Location, Cross, TRS_kg_per_Mg, svp = "symmetrical")

model.gge$SY

(a <- plot(model.gge, type = 1)) # basic plot
(b <- plot(model.gge, type = 2)) # Mean performance vs. stability
(c <- plot(model.gge, type = 3)) # Which-won-where

d <- plot(model.ggeTRS, col.gen = "royalblue", col.env = "seagreen4", size.text.env = 4, plot_theme = theme_metan(grid = "both"), type = 4) # descriminativeness vs. representativeness
print(d)
e <- plot(model.ggeCY, col.gen = "royalblue", col.env = "seagreen4", size.text.env = 4, plot_theme = theme_metan(grid = "both"), type = 4) # descriminativeness vs. representativeness
print(e)
arrange_ggplot(d, e, tag_levels = "a")

(e <- plot(model.gge, type = 5)) # 
(f <- plot(model.gge, type = 6)) # 
(g <- plot(model.gge, type = 7)) # 
(h <- plot(model.gge, type = 8)) # 
(i <- plot(model.gge, type = 9)) #
(i <- plot(model.ggeCY, type = 10)) # 
(i <- plot(model.ggeTRS, type = 10)) # 