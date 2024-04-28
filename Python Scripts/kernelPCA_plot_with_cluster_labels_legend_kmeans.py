#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 11:45:27 2024

@author: zisan
"""

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.decomposition import KernelPCA
import numpy as np
import seaborn as sns

# Load the data
data = pd.read_csv("plus_minus_mumax_trial.csv")

# Remove rows with #N/A values in COG column
data = data.dropna(subset=['COG'])

# Select relevant columns for clustering
X = data[['Depletion_per', 'mumax_plusRha', 'Lag Time','mumax_minusRha','SpikeN']]

# Standardize the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Elbow Method to find the optimal number of clusters
distortions = []
K = range(1, 11)
for k in K:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(X_scaled)
    distortions.append(kmeans.inertia_)

# Plot the Elbow curve
plt.figure(figsize=(4, 4))
plt.plot(K, distortions, 'bx-')
plt.xlabel('Number of Clusters (K)', fontsize=12, fontname='Arial')
plt.ylabel('Distortion', fontsize=12, fontname='Arial')
#plt.title('Elbow Method for Optimal K', fontsize=14, fontname='Arial')
plt.tick_params(axis='both', which='major', labelsize=10)
plt.grid(False)
#plt.savefig("elbow_curve.png", dpi=600)
plt.show()

# Based on the elbow method, let's say the optimal number of clusters is 3
optimal_k = 3

# Perform Kernel PCA
kpca = KernelPCA(n_components=2, kernel='rbf', random_state=12)
X_kpca = kpca.fit_transform(X_scaled)

# Perform KMeans clustering with the optimal number of clusters
kmeans = KMeans(n_clusters=optimal_k, random_state=12)
cluster_labels = kmeans.fit_predict(X_kpca)

# Plot the clustered data points with specified colors
plt.figure(figsize=(4, 4))
colors = {0: '#68217A', 1: 'gold', 2: '#008B8B',3:'darkred',4:'darkblue'}
for cluster in range(optimal_k):
    plt.scatter(X_kpca[cluster_labels == cluster, 0], X_kpca[cluster_labels == cluster, 1],
                label=f'Cluster {cluster+1}', color=colors[cluster])

#plt.title('Clustered Data using Kernel PCA with Elbow Method', fontsize=14, fontname='Arial')
plt.xlabel('Component 1', fontsize=12, fontname='Arial')
plt.ylabel('Component 2', fontsize=12, fontname='Arial')
plt.legend()
plt.grid(False)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.axis('off')

#plt.savefig("clustered_data_kpca.png", dpi=600)

plt.show()

# Add cluster labels to the DataFrame
data['Cluster'] = cluster_labels

# Save the clustered data along with the cluster labels to a CSV file
clustered_data_with_labels = data.copy()
#clustered_data_with_labels.to_csv("clustered_data_with_labels_kernelPCA.csv", index=False)

# Plot the proportion of different COGs in each cluster using circular bar plots
plt.figure(figsize=(15, 6))  # Increase figure size for larger circular bar plots
for i, cluster in enumerate(range(optimal_k)):
    plt.subplot(1, optimal_k, i+1, polar=True)
    cog_counts = data[data['Cluster'] == cluster]['COG'].value_counts(normalize=True).sort_index()
    angles = np.linspace(0, 2 * np.pi, len(cog_counts), endpoint=False).tolist()
    angles += angles[:1]
    cog_counts = cog_counts.tolist()
    cog_counts += cog_counts[:1]
    bars = plt.bar(angles, cog_counts, color=colors[cluster], width=0.4)  # Use color based on cluster
    plt.title(f'Cluster {cluster+1}')
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()

    # Add radial scale
    max_value = max(cog_counts)
    plt.ylim(0, max_value + 0.1)  # Add some padding for better visualization
    plt.yticks(np.linspace(0, max_value, 5), [f'{i:.0%}' for i in np.linspace(0, max_value, 5)])
    plt.gca().spines['polar'].set_visible(False)  # Hide the outer circle

# Adjust spacing between subplots
plt.subplots_adjust(wspace=0)

# Save and show the plot
#plt.savefig("proportion_cogs_cluster.png", dpi=600)
plt.show()

# Print the proportion values for the largest three bars in each cluster
for cluster in range(optimal_k):
    cog_counts = data[data['Cluster'] == cluster]['COG'].value_counts(normalize=True).sort_values(ascending=False)
    largest_three = cog_counts.head(3)
    print(f"Cluster {cluster+1}:")
    for cog, proportion in largest_three.items():
        print(f"   {cog}: {proportion:.2%}")
    print()

# Additional figure with clusters colored based on COG values
plt.figure(figsize=(8, 8))
for cog, color in zip(data['COG'].unique(), ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']):
    plt.scatter(X_kpca[data['COG'] == cog, 0], X_kpca[data['COG'] == cog, 1], label=cog, color=color)
#plt.title('Clustered Data based on Kernel PCA with COG Coloring', fontsize=14, fontname='Arial')
plt.xlabel('Component 1', fontsize=12, fontname='Arial')
plt.ylabel('Component 2', fontsize=12, fontname='Arial')
plt.legend(title='COG', fontsize=8)
plt.grid(False)
plt.tick_params(axis='both', which='major', labelsize=10)
#plt.savefig("clustered_data_cog_color.png", dpi=600)
plt.show()
