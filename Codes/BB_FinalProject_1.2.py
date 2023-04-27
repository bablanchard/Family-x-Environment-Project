# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 13:39:31 2023

@author: BABlanchard
"""



import pyreadr
import csv
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from keras import models, layers
from sklearn import metrics, linear_model
from keras.models import Sequential
from keras.layers import Dense

from sklearn.metrics import classification_report, confusion_matrix, accuracy_score  
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
import time
import seaborn as sns
import os
import itertools
from sklearn.model_selection import GridSearchCV
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import r2_score
from tensorflow.keras.callbacks import EarlyStopping


# Specify path to RDS file
path = r'C:\Users\BABlanchard\OneDrive - LSU AgCenter\Documents\grad semester 7\EXST 7087\Final_Project\BB_pheno_1.1.rds'

# Load data from RDS file into a Pandas dataframe
phenotypic_data = pyreadr.read_r(path)[None]




# Train a random forest model to predict TRS from SNP markers
...

# Print the first few rows of the dataframe
print(phenotypic_data.head())




# Specify path to CSV file
path = 'C:/Users/BABlanchard/OneDrive - LSU AgCenter/Documents/grad semester 7/EXST 7087/Final_Project/ImputedSNPs_original.csv'

# Read CSV file into a pandas dataframe
df = pd.read_csv(path)
# Rename the first column as 'Clone ID'
df = df.rename(columns={'Clone': 'Clone ID'})

# Read CSV file into a list of rows
with open(path, 'r') as file:
    reader = csv.reader(file)
    rows = list(reader)
    
# Get the unique values of the 'Clone ID' column from both dataframes
unique_ids_df = df['Clone ID'].unique()
unique_ids_pheno = phenotypic_data['Clone ID'].unique()

# Print the unique values of 'Clone ID' from both dataframes
print("Unique values in df: ", unique_ids_df)
print("Unique values in phenotypic_data: ", unique_ids_pheno)

# Check if both dataframes have the same unique values of 'Clone ID'
if set(unique_ids_df) == set(unique_ids_pheno):
    print("Both dataframes have the same unique values of Clone ID")
else:
    print("The unique values of Clone ID are different between the two dataframes")

# Merge the dataframes by the 'Clone ID' column and keep only the common observations
merged_df = pd.merge(df, phenotypic_data, on='Clone ID', how='inner')





#%%==== move the variables to the first position in the dataframe; get summary statistics; write csv
cols = list(merged_df)
cols = ['TRS', 'Brix', 'Fiber', 'Diam', 'Bu Wt', 'Count', 'Plot', 'Clone ID', 'Trat', 'Col', 'Row', 'Location', 'Crop', 'PlantYear'] + [col for col in merged_df.columns if col not in ['TRS', 'Brix', 'Fiber', 'Diam', 'Bu Wt', 'Count', 'Plot', 'Clone ID', 'Trat', 'Col', 'Row', 'Location', 'Crop', 'PlantYear']]
merged_df = merged_df[cols]
merged_df = merged_df.loc[:, cols]

# Select the columns to get summary statistics
summary_cols = ['TRS', 'Brix', 'Fiber', 'Diam', 'Bu Wt', 'Count']
summary_df = merged_df[summary_cols].describe()

# Print the summary statistics
print(summary_df)

# Write the summary statistics to a CSV file
summary_df.to_csv('summary_statistics.csv')


##make all needed df to predict
# Remove observations with missing TRS values
TRS = merged_df.dropna(subset=['TRS'])
TRS = TRS.drop(columns = ["Clone ID", "Plot", "Trat", "Brix", "Fiber", "Moisture", "Pol. Reading", "Purity", "Sucrose %", "MD Index", "Spec Residual Outlier", "Bu Wt", "Count", "Height(inches)", "Dia1",
                                "Dia2", "Dia3", "Diam"])

Brix = merged_df.dropna(subset=['Brix'])
Brix = Brix.drop(columns = ["Clone ID", "Plot", "Trat", "TRS", "Fiber", "Moisture", "Pol. Reading", "Purity", "Sucrose %", "MD Index", "Spec Residual Outlier", "Bu Wt", "Count", "Height(inches)", "Dia1",
                                "Dia2", "Dia3", "Diam",])

Fiber = merged_df.dropna(subset=['Fiber'])
Fiber = Fiber.drop(columns = ["Clone ID", "Plot", "Trat", "Brix", "TRS", "Moisture", "Pol. Reading", "Purity", "Sucrose %", "MD Index", "Spec Residual Outlier", "Bu Wt", "Count", "Height(inches)", "Dia1",
                                "Dia2", "Dia3", "Diam",])

Diam = merged_df.dropna(subset=['Diam'])
Diam = Diam.drop(columns = ["Clone ID", "Plot", "Trat", "Brix", "Fiber", "Moisture", "Pol. Reading", "Purity", "Sucrose %", "MD Index", "Spec Residual Outlier", "Bu Wt", "Count", "Height(inches)", "TRS", "Dia1",
                                "Dia2", "Dia3"])


BuWt = merged_df.dropna(subset=['Bu Wt'])
BuWt = BuWt.drop(columns = ["Clone ID", "Plot", "Trat", "Brix", "Fiber", "Moisture", "Pol. Reading", "Purity", "Sucrose %", "MD Index", "Spec Residual Outlier", "Height(inches)", "Count", "TRS", "Dia1",
                                "Dia2", "Dia3", "Diam"])

#%%====results dataframe
# define the columns and rows for the empty dataframe
columns = ['TRS', 'Brix', 'Fiber', 'Diameter', 'BundleWeight']
rows = ['NNMSE', 'NNRMSE', 'NNR2', 'LRMSE', 'LRRMSE', 'LRR2', 'RFMSE', 'RFRMSE', 'RFR2']

# create an empty dataframe with the specified columns and rows
results = pd.DataFrame(index=rows, columns=columns)

# print the resulting dataframe
print(results)

##for each individual trait

#%%====TRS
# Assuming your input data is stored in a pandas DataFrame called 'df'
le = LabelEncoder()
TRS['Location'] = le.fit_transform(TRS['Location'])  # replace 'column_name' with the name of the column containing the non-numeric value
#Split the dataframe into training and testing
X = TRS.drop('TRS', axis=1)
y = TRS['TRS']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# Build neural network model
model = Sequential()
model.add(Dense(1024, activation='relu', input_dim=X_train.shape[1]))
model.add(Dense(512, activation='relu'))
model.add(Dense(256, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(1, activation='linear'))


# Compile model
model.compile(loss='mse', optimizer='adam', metrics=['mse'])

# Define early stopping criteria
early_stopping = EarlyStopping(monitor='loss', patience=50)

# Train the model without a validation set
#model.fit(X_train, y_train, epochs=1000, batch_size=60, callbacks=[early_stopping])
# Train model and save MSE values after each epoch
history = model.fit(X_train, y_train, epochs=1000, batch_size=60)
trsnnmse_values = history.history['mse']

# Evaluate model on test data
test_loss, test_mse = model.evaluate(X_test, y_test)

# Print test loss and mse
print("Test Loss:", test_loss)
print("Test MSE:", test_mse)

# Make predictions on test data
y_pred = model.predict(X_test)

# Calculate R^2 score on test data
r2 = r2_score(y_test, y_pred)
print("R^2 Score:", r2)

# Get feature importances
feature_importances = model.get_weights()[0]


# Print results
import math

print("MSE:", test_mse)
print("Predictions R2:", r2)
test_rmse = math.sqrt(test_mse)
print("Test RMSE:", test_rmse)
# Map feature importances to feature names
importances_by_feature = dict(zip(X_train.columns, feature_importances))

# Print feature importances for each input column
for feature_name, importance in importances_by_feature.items():
    print(f"{feature_name}: {importance}")

# Get indices that sort feature importances in descending order
sorted_idx = np.argsort(feature_importances)[::-1]

# Get names of features sorted by importance
feature_names = X_train.columns[sorted_idx]

# Print top 10 important input columns
print("Top 10 important input columns:")
print(feature_names[:10])



# Linear Regression Model
trslin_reg_mses = []
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    lin_reg = LinearRegression()
    lin_reg.fit(X_train, y_train)
    lin_reg_preds = lin_reg.predict(X_test)
    lin_reg_mse = mean_squared_error(y_test, lin_reg_preds)
    trslin_reg_mses.append(lin_reg_mse)

lrtrsmeanmse = sum(trslin_reg_mses) / 50
lin_reg_rmse = math.sqrt(lrtrsmeanmse)
lin_reg_r2 = r2_score(y_test, lin_reg_preds)
print(f"Linear Regression Mean Squared Error: {lin_reg_mse}")
print(f"Linear Regression Root Mean Squared Error: {lin_reg_rmse}")
print(f"Linear Regression R2 Score: {lin_reg_r2}")





# Random Forest Model
from sklearn.ensemble import RandomForestRegressor

# Define the number of iterations
num_iterations = 5

# Initialize an empty list to store the MSEs for each iteration
rand_forest_mses = []
r2_scores = []

for i in range(num_iterations):
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
    
    # Initialize the random forest model
    rand_forest = RandomForestRegressor(n_estimators=250, random_state=42)
    
    # Train the model on the training data
    rand_forest.fit(X_train, y_train)
    
    # Make predictions on the test data and calculate the MSE
    rand_forest_preds = rand_forest.predict(X_test)
    rand_forest_mse = mean_squared_error(y_test, rand_forest_preds)
    rand_forest_r2 = r2_score(y_test, rand_forest_preds)

    
    # Append the MSE to the list of MSEs
    rand_forest_mses.append(rand_forest_mse)
    r2_scores.append(rand_forest_r2)

# Print the list of MSEs
print(rand_forest_mses)

# Calculate mean MSE of all iterations
mean_rf_mse = sum(rand_forest_mses) / num_iterations
rf_rmse = math.sqrt(mean_rf_mse)
print(f"Random Forest R2 Score: {np.mean(r2_scores)}")

# Print mean MSE
print(f"Mean Random Forest MSE: {mean_rf_mse}")




#fill in TRS results for NN
results.loc['NNMSE', 'TRS'] = test_mse
results.loc['NNRMSE', 'TRS'] = test_rmse
results.loc['NNR2', 'TRS'] = r2

#fill in TRS results for LR
results.loc['LRMSE', 'TRS'] = lrtrsmeanmse
results.loc['LRRMSE', 'TRS'] = lin_reg_rmse
results.loc['LRR2', 'TRS'] = lin_reg_r2

#fill in TRS results for RF
results.loc['RFMSE', 'TRS'] = mean_rf_mse
results.loc['RFRMSE', 'TRS'] = rf_rmse
results.loc['RFR2', 'TRS'] = np.mean(r2_scores)

#perform t-tests to determine significant differences
from scipy.stats import ttest_rel, wilcoxon
from scipy.stats import ttest_ind

# Assume `mse`, `rmse`, and `r2` are arrays of the performance metrics for each model
# Compute Welch's t-test between nn and linreg
t_statistic_mse12, p_value_mse12 = ttest_ind(trsnnmse_values, trslin_reg_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse12}")
print(f"p-value: {p_value_mse12}")

# Compute Welch's t-test between nn and rf
t_statistic_mse13, p_value_mse13 = ttest_ind(trsnnmse_values, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse13}")
print(f"p-value: {p_value_mse13}")
# Compute Welch's t-test between linreg and rf
t_statistic_mse23, p_value_mse23 = ttest_ind(trslin_reg_mses, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse23}")
print(f"p-value: {p_value_mse23}")

#plotting the results
import seaborn as sns
import matplotlib.pyplot as plt

# combine the three arrays into one DataFrame
data = {'Model': ['Neural Network'] * 1000 + ['Linear Regression'] * 50 + ['Random Forest'] * 5,
        'MSE': list(trsnnmse_values) + list(trslin_reg_mses) + list(rand_forest_mses)}

# create box plot
sns.boxplot(x='Model', y='MSE', data=data)
plt.title('Comparison of MSEs for Three Models')
plt.xlabel('Model')
plt.ylabel('Mean Squared Error')

# highlight significant differences
if p_value_mse12 < 0.05:
    plt.plot([0.15, 0.85], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse13 < 0.05:
    plt.plot([-0.15, 0.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(1.0, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse23 < 0.05:
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')

plt.show()

# create a dictionary to hold the results
compare = {
    'Model 1': ['Neural Network', 'Neural Network', 'Linear Regression'],
    'Model 2': ['Linear Regression', 'Random Forest', 'Random Forest'],
    't-statistic': [t_statistic_mse12, t_statistic_mse13, t_statistic_mse23],
    'p-value': [p_value_mse12, p_value_mse13, p_value_mse23]
}

# create a DataFrame from the dictionary
TRSresults_df = pd.DataFrame.from_dict(compare)

# format the p-values as scientific notation with two decimal places
TRSresults_df['p-value'] = TRSresults_df['p-value'].apply(lambda x: f'{x:.2e}')

# display the results
print(TRSresults_df)

TRSresults_df.to_csv('TRScomparison.csv', index=False)

#%%====Brix
# Assuming your input data is stored in a pandas DataFrame called 'df'
le = LabelEncoder()
Brix['Location'] = le.fit_transform(Brix['Location'])  # replace 'column_name' with the name of the column containing the non-numeric value
#Split the dataframe into training and testing
X = Brix.drop('Brix', axis=1)
y = Brix['Brix']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)



#X_train = np.array(X_train).astype(np.float32)
#y_train = np.array(y_train).astype(np.float32)

# Build neural network model
model = Sequential()
model.add(Dense(1024, activation='relu', input_dim=X_train.shape[1]))
model.add(Dense(512, activation='relu'))
model.add(Dense(256, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(1, activation='linear'))


# Compile model
model.compile(loss='mse', optimizer='adam', metrics=['mse'])

# Define early stopping criteria
#early_stopping = EarlyStopping(monitor='loss', patience=50)

# Train the model without a validation set
#model.fit(X_train, y_train, epochs=1000, batch_size=60, callbacks=[early_stopping])
# Train model and save MSE values after each epoch
history = model.fit(X_train, y_train, epochs=1000, batch_size=60)
brixnnmse_values = history.history['mse']

# Evaluate model on test data
test_loss, test_mse = model.evaluate(X_test, y_test)

# Print test loss and mse
print("Test Loss:", test_loss)
print("Test MSE:", test_mse)

# Make predictions on test data
y_pred = model.predict(X_test)

# Calculate R^2 score on test data
r2 = r2_score(y_test, y_pred)
print("R^2 Score:", r2)

# Get feature importances
feature_importances = model.get_weights()[0]


# Print results
import math

print("MSE:", test_mse)
print("Predictions R2:", r2)
test_rmse = math.sqrt(test_mse)
print("Test RMSE:", test_rmse)
# Map feature importances to feature names
importances_by_feature = dict(zip(X_train.columns, feature_importances))

# Print feature importances for each input column
for feature_name, importance in importances_by_feature.items():
    print(f"{feature_name}: {importance}")

# Get indices that sort feature importances in descending order
sorted_idx = np.argsort(feature_importances)[::-1]

# Get names of features sorted by importance
feature_names = X_train.columns[sorted_idx]

# Print top 10 important input columns
print("Top 10 important input columns:")
print(feature_names[:10])



# Linear Regression Model
brixlin_reg_mses = []
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    lin_reg = LinearRegression()
    lin_reg.fit(X_train, y_train)
    lin_reg_preds = lin_reg.predict(X_test)
    lin_reg_mse = mean_squared_error(y_test, lin_reg_preds)
    brixlin_reg_mses.append(lin_reg_mse)

lrbrixmeanmse = sum(brixlin_reg_mses) / 50
lin_reg_rmse = math.sqrt(lrbrixmeanmse)
lin_reg_r2 = r2_score(y_test, lin_reg_preds)
print(f"Linear Regression Mean Squared Error: {lin_reg_mse}")
print(f"Linear Regression Root Mean Squared Error: {lin_reg_rmse}")
print(f"Linear Regression R2 Score: {lin_reg_r2}")




# Random Forest Model
from sklearn.ensemble import RandomForestRegressor

# Define the number of iterations
num_iterations = 5

# Initialize an empty list to store the MSEs for each iteration
rand_forest_mses = []
r2_scores = []

for i in range(num_iterations):
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
    
    # Initialize the random forest model
    rand_forest = RandomForestRegressor(n_estimators=250, random_state=42)
    
    # Train the model on the training data
    rand_forest.fit(X_train, y_train)
    
    # Make predictions on the test data and calculate the MSE
    rand_forest_preds = rand_forest.predict(X_test)
    rand_forest_mse = mean_squared_error(y_test, rand_forest_preds)
    rand_forest_r2 = r2_score(y_test, rand_forest_preds)

    
    # Append the MSE to the list of MSEs
    rand_forest_mses.append(rand_forest_mse)
    r2_scores.append(rand_forest_r2)

# Print the list of MSEs
print(rand_forest_mses)

# Calculate mean MSE of all iterations
mean_rf_mse = sum(rand_forest_mses) / num_iterations
rf_rmse = math.sqrt(mean_rf_mse)

print(f"Random Forest R2 Score: {np.mean(r2_scores)}")

# Print mean MSE
print(f"Mean Random Forest MSE: {mean_rf_mse}")

#fill in TRS results for NN
results.loc['NNMSE', 'Brix'] = test_mse
results.loc['NNRMSE', 'Brix'] = test_rmse
results.loc['NNR2', 'Brix'] = r2

#fill in TRS results for LR
results.loc['LRMSE', 'Brix'] = lrbrixmeanmse
results.loc['LRRMSE', 'Brix'] = lin_reg_rmse
results.loc['LRR2', 'Brix'] = lin_reg_r2

#fill in TRS results for RF
results.loc['RFMSE', 'Brix'] = mean_rf_mse
results.loc['RFRMSE', 'Brix'] = rf_rmse
results.loc['RFR2', 'Brix'] = np.mean(r2_scores)

#perform t-tests to determine significant differences
from scipy.stats import ttest_rel, wilcoxon
from scipy.stats import ttest_ind

# Assume `mse`, `rmse`, and `r2` are arrays of the performance metrics for each model
# Compute Welch's t-test between nn and linreg
t_statistic_mse12, p_value_mse12 = ttest_ind(brixnnmse_values, brixlin_reg_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse12}")
print(f"p-value: {p_value_mse12}")

# Compute Welch's t-test between nn and rf
t_statistic_mse13, p_value_mse13 = ttest_ind(brixnnmse_values, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse13}")
print(f"p-value: {p_value_mse13}")
# Compute Welch's t-test between linreg and rf
t_statistic_mse23, p_value_mse23 = ttest_ind(brixlin_reg_mses, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse23}")
print(f"p-value: {p_value_mse23}")

#plotting the results
import seaborn as sns
import matplotlib.pyplot as plt

# combine the three arrays into one DataFrame
data = {'Model': ['Neural Network'] * 1000 + ['Linear Regression'] * 50 + ['Random Forest'] * 5,
        'MSE': list(brixnnmse_values) + list(brixlin_reg_mses) + list(rand_forest_mses)}

# create box plot
sns.boxplot(x='Model', y='MSE', data=data)
plt.title('Comparison of MSEs for Three Models')
plt.xlabel('Model')
plt.ylabel('Mean Squared Error')

# highlight significant differences
if p_value_mse12 < 0.05:
    plt.plot([0.15, 0.85], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse13 < 0.05:
    plt.plot([-0.15, 0.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(1.0, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse23 < 0.05:
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')

plt.show()

# create a dictionary to hold the results
compare = {
    'Model 1': ['Neural Network', 'Neural Network', 'Linear Regression'],
    'Model 2': ['Linear Regression', 'Random Forest', 'Random Forest'],
    't-statistic': [t_statistic_mse12, t_statistic_mse13, t_statistic_mse23],
    'p-value': [p_value_mse12, p_value_mse13, p_value_mse23]
}

# create a DataFrame from the dictionary
Brixresults_df = pd.DataFrame.from_dict(compare)

# format the p-values as scientific notation with two decimal places
Brixresults_df['p-value'] = Brixresults_df['p-value'].apply(lambda x: f'{x:.2e}')

# display the results
print(Brixresults_df)

Brixresults_df.to_csv('Brixcomparison.csv', index=False)





#%%====Fiber
# Assuming your input data is stored in a pandas DataFrame called 'df'
le = LabelEncoder()
Fiber['Location'] = le.fit_transform(Fiber['Location'])  # replace 'column_name' with the name of the column containing the non-numeric value
#Split the dataframe into training and testing
X = Fiber.drop('Fiber', axis=1)
y = Fiber['Fiber']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)



#X_train = np.array(X_train).astype(np.float32)
#y_train = np.array(y_train).astype(np.float32)

# Build neural network model
model = Sequential()
model.add(Dense(1024, activation='relu', input_dim=X_train.shape[1]))
model.add(Dense(512, activation='relu'))
model.add(Dense(256, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(1, activation='linear'))


# Compile model
model.compile(loss='mse', optimizer='adam', metrics=['mse'])

# Define early stopping criteria
#early_stopping = EarlyStopping(monitor='loss', patience=50)

# Train the model without a validation set
#model.fit(X_train, y_train, epochs=1000, batch_size=60, callbacks=[early_stopping])
# Train model and save MSE values after each epoch
history = model.fit(X_train, y_train, epochs=1000, batch_size=60)
fibnnmse_values = history.history['mse']

# Evaluate model on test data
test_loss, test_mse = model.evaluate(X_test, y_test)

# Print test loss and mse
print("Test Loss:", test_loss)
print("Test MSE:", test_mse)

# Make predictions on test data
y_pred = model.predict(X_test)

# Calculate R^2 score on test data
r2 = r2_score(y_test, y_pred)
print("R^2 Score:", r2)

# Get feature importances
feature_importances = model.get_weights()[0]


# Print results
import math

print("MSE:", test_mse)
print("Predictions R2:", r2)
test_rmse = math.sqrt(test_mse)
print("Test RMSE:", test_rmse)




# Linear Regression Model
fiblin_reg_mses = []
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    lin_reg = LinearRegression()
    lin_reg.fit(X_train, y_train)
    lin_reg_preds = lin_reg.predict(X_test)
    lin_reg_mse = mean_squared_error(y_test, lin_reg_preds)
    fiblin_reg_mses.append(lin_reg_mse)

lrfibmeanmse = sum(fiblin_reg_mses) / 50
lin_reg_rmse = math.sqrt(lrfibmeanmse)
lin_reg_r2 = r2_score(y_test, lin_reg_preds)
print(f"Linear Regression Mean Squared Error: {lin_reg_mse}")
print(f"Linear Regression Root Mean Squared Error: {lin_reg_rmse}")
print(f"Linear Regression R2 Score: {lin_reg_r2}")





# Random Forest Model
from sklearn.ensemble import RandomForestRegressor

# Define the number of iterations
num_iterations = 5

# Initialize an empty list to store the MSEs for each iteration
rand_forest_mses = []
r2_scores = []

for i in range(num_iterations):
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
    
    # Initialize the random forest model
    rand_forest = RandomForestRegressor(n_estimators=250, random_state=42)
    
    # Train the model on the training data
    rand_forest.fit(X_train, y_train)
    
    # Make predictions on the test data and calculate the MSE
    rand_forest_preds = rand_forest.predict(X_test)
    rand_forest_mse = mean_squared_error(y_test, rand_forest_preds)
    rand_forest_r2 = r2_score(y_test, rand_forest_preds)

    
    # Append the MSE to the list of MSEs
    rand_forest_mses.append(rand_forest_mse)
    r2_scores.append(rand_forest_r2)

# Print the list of MSEs
print(rand_forest_mses)

# Calculate mean MSE of all iterations
mean_rf_mse = sum(rand_forest_mses) / num_iterations
rf_rmse = math.sqrt(mean_rf_mse)

print(f"Random Forest R2 Score: {np.mean(r2_scores)}")

# Print mean MSE
print(f"Mean Random Forest MSE: {mean_rf_mse}")

#fill in Fiber results for NN
results.loc['NNMSE', 'Fiber'] = test_mse
results.loc['NNRMSE', 'Fiber'] = test_rmse
results.loc['NNR2', 'Fiber'] = r2

#fill in TRS results for LR
results.loc['LRMSE', 'Fiber'] = lrfibmeanmse
results.loc['LRRMSE', 'Fiber'] = lin_reg_rmse
results.loc['LRR2', 'Fiber'] = lin_reg_r2

#fill in TRS results for RF
results.loc['RFMSE', 'Fiber'] = mean_rf_mse
results.loc['RFRMSE', 'Fiber'] = rf_rmse
results.loc['RFR2', 'Fiber'] = np.mean(r2_scores)

#perform t-tests to determine significant differences
from scipy.stats import ttest_rel, wilcoxon
from scipy.stats import ttest_ind

# Assume `mse`, `rmse`, and `r2` are arrays of the performance metrics for each model
# Compute Welch's t-test between nn and linreg
t_statistic_mse12, p_value_mse12 = ttest_ind(fibnnmse_values, fiblin_reg_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse12}")
print(f"p-value: {p_value_mse12}")

# Compute Welch's t-test between nn and rf
t_statistic_mse13, p_value_mse13 = ttest_ind(fibnnmse_values, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse13}")
print(f"p-value: {p_value_mse13}")
# Compute Welch's t-test between linreg and rf
t_statistic_mse23, p_value_mse23 = ttest_ind(fiblin_reg_mses, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse23}")
print(f"p-value: {p_value_mse23}")

#plotting the results
import seaborn as sns
import matplotlib.pyplot as plt

# combine the three arrays into one DataFrame
data = {'Model': ['Neural Network'] * 1000 + ['Linear Regression'] * 50 + ['Random Forest'] * 5,
        'MSE': list(fibnnmse_values) + list(fiblin_reg_mses) + list(rand_forest_mses)}

# create box plot
sns.boxplot(x='Model', y='MSE', data=data)
plt.title('Comparison of MSEs for Three Models')
plt.xlabel('Model')
plt.ylabel('Mean Squared Error')

# highlight significant differences
if p_value_mse12 < 0.05:
    plt.plot([0.15, 0.85], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse13 < 0.05:
    plt.plot([-0.15, 0.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(1.0, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse23 < 0.05:
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')

plt.show()

# create a dictionary to hold the results
compare = {
    'Model 1': ['Neural Network', 'Neural Network', 'Linear Regression'],
    'Model 2': ['Linear Regression', 'Random Forest', 'Random Forest'],
    't-statistic': [t_statistic_mse12, t_statistic_mse13, t_statistic_mse23],
    'p-value': [p_value_mse12, p_value_mse13, p_value_mse23]
}

# create a DataFrame from the dictionary
Fibresults_df = pd.DataFrame.from_dict(compare)

# format the p-values as scientific notation with two decimal places
Fibresults_df['p-value'] = Fibresults_df['p-value'].apply(lambda x: f'{x:.2e}')

# display the results
print(Fibresults_df)

Fibresults_df.to_csv('Fibercomparison.csv', index=False)
#%%====Diameter
# Assuming your input data is stored in a pandas DataFrame called 'df'
le = LabelEncoder()
Diam['Location'] = le.fit_transform(Diam['Location'])  # replace 'column_name' with the name of the column containing the non-numeric value
#Split the dataframe into training and testing
X = Diam.drop('Diam', axis=1)
y = Diam['Diam']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)



#X_train = np.array(X_train).astype(np.float32)
#y_train = np.array(y_train).astype(np.float32)

# Build neural network model
model = Sequential()
model.add(Dense(64, activation='relu', input_dim=X_train.shape[1]))
model.add(Dense(32, activation='relu'))
model.add(Dense(1, activation='linear'))



# Compile model
model.compile(loss='mse', optimizer='adam', metrics=['mse'])

# Define early stopping criteria
#early_stopping = EarlyStopping(monitor='loss', patience=50)

# Train the model without a validation set
#model.fit(X_train, y_train, epochs=1000, batch_size=130, callbacks=[early_stopping])
# Train model
history = model.fit(X_train, y_train, epochs=300, batch_size=130)
diannmse_values = history.history['mse']

# Evaluate model on test data
test_loss, test_mse = model.evaluate(X_test, y_test)

# Print test loss and mse
print("Test Loss:", test_loss)
print("Test MSE:", test_mse)

# Make predictions on test data
y_pred = model.predict(X_test)

# Calculate R^2 score on test data
r2 = r2_score(y_test, y_pred)
print("R^2 Score:", r2)

# Get feature importances
feature_importances = model.get_weights()[0]


# Print results
import math

print("MSE:", test_mse)
print("Predictions R2:", r2)
test_rmse = math.sqrt(test_mse)
print("Test RMSE:", test_rmse)




# Linear Regression Model
dialin_reg_mses = []
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    lin_reg = LinearRegression()
    lin_reg.fit(X_train, y_train)
    lin_reg_preds = lin_reg.predict(X_test)
    lin_reg_mse = mean_squared_error(y_test, lin_reg_preds)
    dialin_reg_mses.append(lin_reg_mse)

lrdiameanmse = sum(dialin_reg_mses) / 50
lin_reg_rmse = math.sqrt(lrdiameanmse)
lin_reg_r2 = r2_score(y_test, lin_reg_preds)
print(f"Linear Regression Mean Squared Error: {lin_reg_mse}")
print(f"Linear Regression Root Mean Squared Error: {lin_reg_rmse}")
print(f"Linear Regression R2 Score: {lin_reg_r2}")

import matplotlib.pyplot as plt

plt.hist(y, bins=50)
plt.xlabel('Target')
plt.ylabel('Frequency')
plt.show()





# Random Forest Model
from sklearn.ensemble import RandomForestRegressor

# Define the number of iterations
num_iterations = 5

# Initialize an empty list to store the MSEs for each iteration
rand_forest_mses = []
r2_scores = []

for i in range(num_iterations):
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
    
    # Initialize the random forest model
    rand_forest = RandomForestRegressor(n_estimators=250, random_state=42)
    
    # Train the model on the training data
    rand_forest.fit(X_train, y_train)
    
    # Make predictions on the test data and calculate the MSE
    rand_forest_preds = rand_forest.predict(X_test)
    rand_forest_mse = mean_squared_error(y_test, rand_forest_preds)
    rand_forest_r2 = r2_score(y_test, rand_forest_preds)

    
    # Append the MSE to the list of MSEs
    rand_forest_mses.append(rand_forest_mse)
    r2_scores.append(rand_forest_r2)

# Print the list of MSEs
print(rand_forest_mses)

# Calculate mean MSE of all iterations
mean_rf_mse = sum(rand_forest_mses) / num_iterations
rf_rmse = math.sqrt(mean_rf_mse)

print(f"Random Forest R2 Score: {np.mean(r2_scores)}")

# Print mean MSE
print(f"Mean Random Forest MSE: {mean_rf_mse}")

#fill in Diameter results for NN
results.loc['NNMSE', 'Diameter'] = test_mse
results.loc['NNRMSE', 'Diameter'] = test_rmse
results.loc['NNR2', 'Diameter'] = r2

#fill in Diameter results for LR
results.loc['LRMSE', 'Diameter'] = lrdiameanmse
results.loc['LRRMSE', 'Diameter'] = lin_reg_rmse
results.loc['LRR2', 'Diameter'] = lin_reg_r2

#fill in Diameter results for RF
results.loc['RFMSE', 'Diameter'] = mean_rf_mse
results.loc['RFRMSE', 'Diameter'] = rf_rmse
results.loc['RFR2', 'Diameter'] = np.mean(r2_scores)

#perform t-tests to determine significant differences
from scipy.stats import ttest_rel, wilcoxon
from scipy.stats import ttest_ind

# Assume `mse`, `rmse`, and `r2` are arrays of the performance metrics for each model
# Compute Welch's t-test between nn and linreg
t_statistic_mse12, p_value_mse12 = ttest_ind(diannmse_values, dialin_reg_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse12}")
print(f"p-value: {p_value_mse12}")

# Compute Welch's t-test between nn and rf
t_statistic_mse13, p_value_mse13 = ttest_ind(diannmse_values, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse13}")
print(f"p-value: {p_value_mse13}")
# Compute Welch's t-test between linreg and rf
t_statistic_mse23, p_value_mse23 = ttest_ind(dialin_reg_mses, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse23}")
print(f"p-value: {p_value_mse23}")

#plotting the results
import seaborn as sns
import matplotlib.pyplot as plt

# combine the three arrays into one DataFrame
data = {'Model': ['Neural Network'] * 300 + ['Linear Regression'] * 50 + ['Random Forest'] * 5,
        'MSE': list(diannmse_values) + list(dialin_reg_mses) + list(rand_forest_mses)}

# create box plot
sns.boxplot(x='Model', y='MSE', data=data)
plt.title('Comparison of MSEs for Three Models')
plt.xlabel('Model')
plt.ylabel('Mean Squared Error')

# highlight significant differences
if p_value_mse12 < 0.05:
    plt.plot([0.15, 0.85], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse13 < 0.05:
    plt.plot([-0.15, 0.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(1.0, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse23 < 0.05:
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')

plt.show()
# create a dictionary to hold the results
compare = {
    'Model 1': ['Neural Network', 'Neural Network', 'Linear Regression'],
    'Model 2': ['Linear Regression', 'Random Forest', 'Random Forest'],
    't-statistic': [t_statistic_mse12, t_statistic_mse13, t_statistic_mse23],
    'p-value': [p_value_mse12, p_value_mse13, p_value_mse23]
}

# create a DataFrame from the dictionary
Diaresults_df = pd.DataFrame.from_dict(compare)

# format the p-values as scientific notation with two decimal places
Diaresults_df['p-value'] = Diaresults_df['p-value'].apply(lambda x: f'{x:.2e}')

# display the results
print(Diaresults_df)

Diaresults_df.to_csv('Diametercomparison.csv', index=False)
#%%====Bundle Weight
# Assuming your input data is stored in a pandas DataFrame called 'df'
le = LabelEncoder()
BuWt['Location'] = le.fit_transform(BuWt['Location'])  # replace 'column_name' with the name of the column containing the non-numeric value
#Split the dataframe into training and testing
X = BuWt.drop('Bu Wt', axis=1)
y = BuWt['Bu Wt']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)



#X_train = np.array(X_train).astype(np.float32)
#y_train = np.array(y_train).astype(np.float32)

# Build neural network model
model = Sequential()
model.add(Dense(1024, activation='relu', input_dim=X_train.shape[1]))
model.add(Dense(512, activation='relu'))
model.add(Dense(256, activation='relu'))
model.add(Dense(128, activation='relu'))
model.add(Dense(64, activation='relu'))
model.add(Dense(1, activation='linear'))



# Compile model
model.compile(loss='mse', optimizer='adam', metrics=['mse'])

# Define early stopping criteria
early_stopping = EarlyStopping(monitor='loss', patience=50)

# Train the model without a validation set
#model.fit(X_train, y_train, epochs=1000, batch_size=60, callbacks=[early_stopping])
# Train model
history = model.fit(X_train, y_train, epochs=300, batch_size=50)
buwtnnmse_values = history.history['mse']

# Evaluate model on test data
test_loss, test_mse = model.evaluate(X_test, y_test)

# Print test loss and mse
print("Test Loss:", test_loss)
print("Test MSE:", test_mse)

# Make predictions on test data
y_pred = model.predict(X_test)

# Calculate R^2 score on test data
r2 = r2_score(y_test, y_pred)
print("R^2 Score:", r2)

# Get feature importances
feature_importances = model.get_weights()[0]


# Print results
import math

print("MSE:", test_mse)
print("Predictions R2:", r2)
test_rmse = math.sqrt(test_mse)
print("Test RMSE:", test_rmse)




# Linear Regression Model
buwtlin_reg_mses = []
for i in range(50):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    lin_reg = LinearRegression()
    lin_reg.fit(X_train, y_train)
    lin_reg_preds = lin_reg.predict(X_test)
    lin_reg_mse = mean_squared_error(y_test, lin_reg_preds)
    buwtlin_reg_mses.append(lin_reg_mse)

lrbuwtmeanmse = sum(buwtlin_reg_mses) / 50
lin_reg_rmse = math.sqrt(lrbuwtmeanmse)
lin_reg_r2 = r2_score(y_test, lin_reg_preds)
print(f"Linear Regression Mean Squared Error: {lin_reg_mse}")
print(f"Linear Regression Root Mean Squared Error: {lin_reg_rmse}")
print(f"Linear Regression R2 Score: {lin_reg_r2}")





# Random Forest Model

# Define the number of iterations
num_iterations = 5

# Initialize an empty list to store the MSEs for each iteration
rand_forest_mses = []
r2_scores = []

for i in range(num_iterations):
    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=i)
    
    # Initialize the random forest model
    rand_forest = RandomForestRegressor(n_estimators=250, random_state=42)
    
    # Train the model on the training data
    rand_forest.fit(X_train, y_train)
    
    # Make predictions on the test data and calculate the MSE
    rand_forest_preds = rand_forest.predict(X_test)
    rand_forest_mse = mean_squared_error(y_test, rand_forest_preds)
    rand_forest_r2 = r2_score(y_test, rand_forest_preds)

    
    # Append the MSE to the list of MSEs
    rand_forest_mses.append(rand_forest_mse)
    r2_scores.append(rand_forest_r2)

# Print the list of MSEs
print(rand_forest_mses)

# Calculate mean MSE of all iterations
mean_rf_mse = sum(rand_forest_mses) / num_iterations
rf_rmse = math.sqrt(mean_rf_mse)

print(f"Random Forest R2 Score: {np.mean(r2_scores)}")

# Print mean MSE
print(f"Mean Random Forest MSE: {mean_rf_mse}")

#fill in Bundle Weight results for NN
results.loc['NNMSE', 'BundleWeight'] = test_mse
results.loc['NNRMSE', 'BundleWeight'] = test_rmse
results.loc['NNR2', 'BundleWeight'] = r2

#fill in TRS results for LR
results.loc['LRMSE', 'BundleWeight'] = lrbuwtmeanmse
results.loc['LRRMSE', 'BundleWeight'] = lin_reg_rmse
results.loc['LRR2', 'BundleWeight'] = lin_reg_r2

#fill in TRS results for RF
results.loc['RFMSE', 'BundleWeight'] = mean_rf_mse
results.loc['RFRMSE', 'BundleWeight'] = rf_rmse
results.loc['RFR2', 'BundleWeight'] = np.mean(r2_scores)


#perform t-tests to determine significant differences
from scipy.stats import ttest_rel, wilcoxon
from scipy.stats import ttest_ind

# Assume `mse`, `rmse`, and `r2` are arrays of the performance metrics for each model
# Compute Welch's t-test between nn and linreg
t_statistic_mse12, p_value_mse12 = ttest_ind(buwtnnmse_values, buwtlin_reg_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse12}")
print(f"p-value: {p_value_mse12}")

# Compute Welch's t-test between nn and rf
t_statistic_mse13, p_value_mse13 = ttest_ind(buwtnnmse_values, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse13}")
print(f"p-value: {p_value_mse13}")
# Compute Welch's t-test between linreg and rf
t_statistic_mse23, p_value_mse23 = ttest_ind(buwtlin_reg_mses, rand_forest_mses, equal_var=False)
# Print results
print("Welch's t-test")
print(f"t-statistic: {t_statistic_mse23}")
print(f"p-value: {p_value_mse23}")

#plotting the results
import seaborn as sns
import matplotlib.pyplot as plt

# combine the three arrays into one DataFrame
data = {'Model': ['Neural Network'] * 300 + ['Linear Regression'] * 50 + ['Random Forest'] * 5,
        'MSE': list(buwtnnmse_values) + list(buwtlin_reg_mses) + list(rand_forest_mses)}

# create box plot
sns.boxplot(x='Model', y='MSE', data=data)
plt.title('Comparison of MSEs for Three Models')
plt.xlabel('Model')
plt.ylabel('Mean Squared Error')

# highlight significant differences
if p_value_mse12 < 0.05:
    plt.plot([0.15, 0.85], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse13 < 0.05:
    plt.plot([-0.15, 0.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(1.0, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')
if p_value_mse23 < 0.05:
    plt.plot([0.85, 1.15], [max(data['MSE']) + 1, max(data['MSE']) + 1], lw=2, c='k')
    plt.text(0.5, max(data['MSE']) + 1.5, 'p < 0.05', ha='center')

plt.show()
# create a dictionary to hold the results
compare = {
    'Model 1': ['Neural Network', 'Neural Network', 'Linear Regression'],
    'Model 2': ['Linear Regression', 'Random Forest', 'Random Forest'],
    't-statistic': [t_statistic_mse12, t_statistic_mse13, t_statistic_mse23],
    'p-value': [p_value_mse12, p_value_mse13, p_value_mse23]
}

# create a DataFrame from the dictionary
Buwtresults_df = pd.DataFrame.from_dict(compare)

# format the p-values as scientific notation with two decimal places
Buwtresults_df['p-value'] = Buwtresults_df['p-value'].apply(lambda x: f'{x:.2e}')

# display the results
print(Buwtresults_df)

Buwtresults_df.to_csv('BundleWeightcomparison.csv', index=False)


#%%====Results graph


# Define colors for each bar
colors = [['#1f77b4', '#ff7f0e', '#2ca02c'],
          ['#ff7f0e', '#1f77b4', '#2ca02c', 'gray', 'purple'],
          ['#1f77b4', '#1f77b4', '#2ca02c', 'red', 'purple'],
          ['#2ca02c', '#2ca02c', '#2ca02c', 'gray', 'gray'],
          ['#ff7f0e', '#1f77b4', '#2ca02c', 'gray', 'gray']]

# Filter the dataframe to include only NNMSE, LRMSE, and RFMSE
nnmse_lrmse_rfms = results.loc[['NNMSE', 'LRMSE', 'RFMSE']]
nnmse_lrmse_rfms.loc['LRMSE', 'Diameter'] = np.nan

# Set up the subplots
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 5), sharey=False)

x_labels = ["NN", "LR", "RF"]

# Create a bar graph for each column
for i, col in enumerate(nnmse_lrmse_rfms.columns):
    ax = axs[i]
    for j, row in enumerate(nnmse_lrmse_rfms.index):
        color = colors[i][j]
        value = nnmse_lrmse_rfms.loc[row, col]
        ax.bar(x_labels[j], value, color=color)
        ax.text(x_labels[j], value, round(value, 2), ha='center', va='bottom', fontsize=8)
    ax.set_title(col)

# Add some overall plot titles and axis labels
fig.suptitle('MSE by Model for each Trait')
fig.text(0.5, 0.04, 'Models', ha='center')
fig.text(0.04, 0.5, 'Mean Squared Error', va='center', rotation='vertical')

# Display the plot
plt.show()


# Filter the dataframe to include only NNMSE, LRMSE, and RFMSE
nnr2_lrr2_rfr2 = results.loc[['NNR2', 'LRR2', 'RFR2']]
nnr2_lrr2_rfr2.loc['LRR2', 'Diameter'] = np.nan


# Define colors for each bar
colors = [['#1f77b4', '#1f77b4', '#1f77b4'],
          ['#1f77b4', '#1f77b4', '#1f77b4', 'gray', 'purple'],
          ['#1f77b4', '#1f77b4', '#1f77b4', 'red', 'purple'],
          ['#1f77b4', '#1f77b4', '#1f77b4', 'gray', 'gray'],
          ['#1f77b4', '#1f77b4', '#1f77b4', 'gray', 'gray']]



# Set up the subplots
fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15, 5), sharey=False)
x_labels = ["NN", "LR", "RF"]

# Create a bar graph for each column
for i, col in enumerate(nnr2_lrr2_rfr2.columns):
    ax = axs[i]
    for j, row in enumerate(nnr2_lrr2_rfr2.index):
        color = colors[i][j]
        value = nnr2_lrr2_rfr2.loc[row, col]
        ax.bar(x_labels[j], value, color=color)
        ax.text(x_labels[j], value, round(value, 2), ha='center', va='bottom', fontsize=8)
    ax.set_title(col)

# Add some overall plot titles and axis labels
fig.suptitle('R2 by Model for each Trait')
fig.text(0.5, 0.04, 'Models', ha='center')
fig.text(0.04, 0.5, 'R2', va='center', rotation='vertical')

# Display the plot
plt.show()




