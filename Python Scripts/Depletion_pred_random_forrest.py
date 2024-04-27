#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 14:08:29 2024

@author: zisan
"""

import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

# Read the data from the CSV file
data = pd.read_csv("plus_minus_mumax_trial.csv")

# Check for missing values in the columns of interest
if data[['Depletion_per', 'mumax_plusRha', 'Lag Time']].isnull().values.any():
    print("Dataset contains missing values in the columns of interest.")
else:
    # Define the independent variables (features) and the dependent variable (target)
    X = data[['mumax_plusRha', 'Lag Time']]
    y = data['Depletion_per']

    # Perform feature engineering
    poly = PolynomialFeatures(degree=2, include_bias=False)
    X_poly = poly.fit_transform(X)
    X_poly = pd.DataFrame(X_poly, columns=poly.get_feature_names_out(X.columns))

    # Fit a more sophisticated model
    model_rf = RandomForestRegressor(n_estimators=100, random_state=42)
    model_rf.fit(X_poly, y)

    # Display model parameters
    print("Random Forest Regression Model Parameters:")
    print("Number of Estimators:", model_rf.n_estimators)
    print("Random State:", model_rf.random_state)
    print("Criterion:", model_rf.criterion)
    # Add more parameters as needed

    # Calculate predicted values
    y_pred_rf = model_rf.predict(X_poly)

    # Calculate the mean squared error
    mse_rf = mean_squared_error(y, y_pred_rf)

    # Plot actual vs. predicted values
    plt.figure(figsize=(8, 6))
    plt.scatter(y, y_pred_rf, color='blue', label='Actual vs. Predicted')
    plt.plot(y, y, color='red', label='Perfect Prediction')
    plt.title('Actual vs. Predicted Depletion_per (Random Forest Regression)')
    plt.xlabel('Actual Depletion_per')
    plt.ylabel('Predicted Depletion_per')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Add predicted values to the DataFrame
    data['Predicted_Depletion_per'] = y_pred_rf

   # Get the coefficients of the polynomial features
    poly_coefficients = model_rf.feature_importances_


# Get the feature names
feature_names = X_poly.columns

# Define alphabet coefficients
alphabets = ['c0'] + [f'c{i}' for i in range(1, len(feature_names) + 1)]

# Print the polynomial equation with variables and alphabets
print("\nPolynomial Equation:")
equation = "Depletion_per = "
for i, coef in enumerate(poly_coefficients):
    if i == 0:
        equation += alphabets[i]  # Intercept term
    else:
        variable = feature_names[i-1]
        equation += f" + {alphabets[i]} * {variable}"
print(equation)

# Print the coefficients with corresponding variables
print("Coefficients with Variables:")
for i, coef in enumerate(poly_coefficients):
    if i == 0:
        print(f"{alphabets[i]} (Intercept): {coef:.2f}")
    else:
        variable = feature_names[i-1]
        print(f"{alphabets[i]} (Coefficient for {variable}): {coef:.2f}")



    # Save the DataFrame to a CSV file
    # data.to_csv("predicted_data.csv", index=False)

    # Print model performance
print("\nRandom Forest Regression MSE:", mse_rf)