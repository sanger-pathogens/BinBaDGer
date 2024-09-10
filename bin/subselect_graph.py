#!/usr/bin/env python3

import argparse
import networkx as nx
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import umap

from sklearn.manifold import MDS
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from networkx.algorithms.community import greedy_modularity_communities

def read_phylip_distance(phylip_file):
    with open(phylip_file, 'r') as file:
        lines = file.readlines()

        #might be useful
        num_samples = int(lines[0].strip())

        accessions = []
        matrix = []

        for line in lines[1:]:
            parts = line.strip().split()
            accessions.append(parts[0])
            matrix.append([float(x) for x in parts[1:]])

        matrix = np.array(matrix)

    return num_samples, accessions, matrix

def hierarchy_cluster(matrix, accessions, num_clusters, output_file, plot_clusters, method='average'):
    condensed_matrix = squareform(matrix)

    Z = linkage(condensed_matrix, method=method)

    clusters = fcluster(Z, t=num_clusters, criterion='maxclust')

    if plot_clusters:
        plot_dendogram(Z, accessions, 'dendogram.png')
    
    return clusters

def optimal_number_of_clusters(matrix, max_clusters=30):
    """
    Determine the optimal number of clusters using silhouette analysis.
    Parameters:
        D (np.array): Pairwise distance matrix (n x n).
        max_clusters (int): Maximum number of clusters to evaluate.
    Returns:
        optimal_k (int): Optimal number of clusters.
    """
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
    X_transformed = mds.fit_transform(matrix)
    silhouette_scores = []
    
    for n_clusters in range(2, max_clusters + 1):
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        labels = kmeans.fit_predict(X_transformed)
        silhouette_avg = silhouette_score(X_transformed, labels)
        silhouette_scores.append(silhouette_avg)
        print(f"For n_clusters = {n_clusters}, the silhouette score is {silhouette_avg:.4f}")
    
    
    # Plot silhouette scores to visualize the optimal number of clusters
    plt.figure(figsize=(15, 9))
    plt.plot(range(2, max_clusters + 1), silhouette_scores, marker='o')
    plt.title('Silhouette Analysis for Optimal Number of Clusters')
    plt.xlabel('Number of clusters')
    plt.ylabel('Silhouette score')
    plt.grid(True)
    plt.savefig('optimal_k.png')
    optimal_k = np.argmax(silhouette_scores) + 2  # +2 because range starts from 2
    print(f"The optimal number of clusters is {optimal_k} with a silhouette score of {silhouette_scores[optimal_k-2]:.4f}")
    plt.close()

    return X_transformed, optimal_k


def kmeans_cluster(matrix):

    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)

    transformed_matrix, optimal_k = optimal_number_of_clusters(matrix)

    kmeans = KMeans(n_clusters=optimal_k, random_state=42)

    kmeans.fit(transformed_matrix)

    labels = kmeans.labels_

    plt.figure(figsize=(15, 9))
    plt.scatter(transformed_matrix[:, 0], transformed_matrix[:, 1], c=labels, cmap='tab20', marker='o', s=100, edgecolor='k')
    plt.title("Kmeans Clustering")
    plt.ylabel("MD1")
    plt.xlabel("MD2")
    plt.colorbar(label='Cluster label')

    plt.savefig('Kmeans.png', dpi=300, bbox_inches='tight')
    plt.close()

    save_clusters(transformed_matrix, labels, 'kmeans_clusters.csv')

    return labels

def save_clusters(transformed_matrix, labels, filename='cluster_assignments.csv'):
    df = pd.DataFrame({
        'MDS Dimension 1': transformed_matrix[:, 0],
        'MDS Dimension 2': transformed_matrix[:, 1],
        'Cluster Label': labels
    })

    df.to_csv(filename, index=False)


def plot_dendogram(Z, accessions, output_file):
    plt.figure(figsize=(15, 9))
    dendrogram(Z, labels=accessions, leaf_font_size=5, orientation='left')
    
    plt.title("Dendrogram of Hierarchical Clustering")
    plt.ylabel("Accessions")
    plt.xlabel("Distance")
    plt.tight_layout()
    
    # Save the plot as a PNG file
    plt.savefig(output_file)
    plt.close()

def edge_based_cluster(matrix, min, max, accessions):
    G = nx.Graph()

    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if min <= matrix[i, j] <= max:
                G.add_edge(i, j, weight=matrix[i, j])

    clusters = greedy_modularity_communities(G)
    
    return clusters
                

def relate_id_to_accession(clusters, accessions, output_file):
    cluster_dict = {}
    for idx, cluster in enumerate(clusters):
        if cluster not in cluster_dict:
            cluster_dict[cluster] = []
        cluster_dict[cluster].append(accessions[idx])
        print(f'Cluster {idx + 1}: {sorted(cluster)}')
    
    with open(output_file, 'w') as file:
        for cluster, members in cluster_dict.items():
            file.write(f"cluster {cluster}: {', '.join(members)}\n")

def plot_umap(matrix, labels, output_file):
    """
    Perform UMAP and plot clusters.

    Parameters:
        matrix (np.array): Pairwise distance matrix (n x n).
        labels (np.array): Cluster labels.

    Returns:
        X_umap (np.array): Data transformed by UMAP.
    """
    X_umap = umap.UMAP(metric='precomputed', random_state=42).fit_transform(matrix)
   
    plt.figure(figsize=(15, 9))
    scatter = plt.scatter(X_umap[:, 0], X_umap[:, 1], c=labels, cmap='tab20', marker='o', s=100, edgecolor='k')
    plt.title('UMAP Clustering')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.colorbar(label='Cluster Label')
    plt.savefig(output_file)
    plt.close()
   
    return X_umap

def umap_clustering(matrix):
    clusterable_embedding= umap.UMAP(metric='precomputed', min_dist=0.0, n_components=2, n_neighbors=30, random_state=42).fit_transform(matrix)

    labels = HDBSCAN(min_cluster_size=5).fit_predict(clusterable_embedding)

    clustered = (labels >= 0)

    plt.figure(figsize=(15, 9))
    scatter = plt.scatter(clusterable_embedding[:, 0], clusterable_embedding[:, 1], c=clustered, cmap='tab20', marker='o', s=100, edgecolor='k')
    plt.title('UMAP Clustering')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.colorbar(label='Cluster Label')
    plt.savefig('cluster_embedding.png')
    plt.close()

    return clustered

def select_closest_representatives(matrix, clusters, names, n_representatives=3):
    representatives = []
    unique_clusters = np.unique(clusters)
   
    for cluster in unique_clusters:
        indices = np.where(clusters == cluster)[0]
        sub_matrix = matrix[indices][:, indices]
        avg_distances = np.mean(sub_matrix, axis=1)
       
        # Get indices of the closest n representatives
        closest_indices = indices[np.argsort(avg_distances)[:n_representatives]]
        representatives.extend([names[i] for i in closest_indices])
   
    return representatives

def main():
    parser = argparse.ArgumentParser(description='subsample from a matrix')
    parser.add_argument("filename", type=str, help='Path to the Phylip file')

    args = parser.parse_args()

    _, accessions, matrix = read_phylip_distance(args.filename)

    #method = 'kmeans'
    
    #if method == 'kmeans':
    kmeans_labels = kmeans_cluster(matrix)
    plot_umap(matrix, kmeans_labels, 'kmeans_labels_umap.png')

    #if method == 'fcluster':
    hier_labels = hierarchy_cluster(matrix, accessions, 3, "clusters.txt", True, 'average')
    plot_umap(matrix, hier_labels, 'hier_labels_umap.png')
    

    hdbscan_labels = umap_clustering(matrix)
    plot_umap(matrix, hdbscan_labels, 'hdbscan_labels_umap.png')

    clusters = edge_based_cluster(matrix, 0.0001, 0.01, accessions)

    relate_id_to_accession(clusters, accessions, "clusters.txt")

    reps = select_closest_representatives(matrix, kmeans_labels, accessions, n_representatives=3)
    print("Selected Representatives:", reps)

if __name__ == "__main__":
    main()
