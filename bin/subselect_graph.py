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

##### input parsing functions #####

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

##### clustering functions #####

### hierarchy ###

def hierarchy_cluster(matrix, accessions, num_clusters=3, plot_clusters=True, method='average'):
    condensed_matrix = squareform(matrix)

    Z = linkage(condensed_matrix, method=method)

    clusters = fcluster(Z, t=num_clusters, criterion='maxclust')

    if plot_clusters:
        plot_dendogram(Z, accessions, 'dendogram.png')
    
    return clusters

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

### kmeans ###

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
    plt.close()

    return X_transformed, optimal_k

### edge_based ###

def edge_based_cluster(matrix, accessions, threshold):
    G = nx.Graph()
    
    num_nodes = len(matrix)
    
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if matrix[i, j] < threshold:
                G.add_edge(i, j, weight=matrix[i, j])
    
    clusters = greedy_modularity_communities(G)
    
    plot_network_subclusters(clusters, G)

    return clusters

### umap/HDBSCAN_based ###

def umap_clustering(matrix):
    clusterable_embedding = umap.UMAP(
        metric='precomputed', 
        min_dist=0.0, 
        n_components=2, 
        n_neighbors=30, 
        random_state=42
        ).fit_transform(matrix)

    labels = HDBSCAN(min_cluster_size=5, min_samples=10,).fit_predict(clusterable_embedding)

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

###### result plotting functions ######

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

def plot_network_subclusters(clusters, G):
    num_clusters = len(clusters)
    
    # Determine grid size
    grid_size = int(np.ceil(np.sqrt(num_clusters)))
    
    # Create subplots
    fig, axes = plt.subplots(grid_size, grid_size, figsize=(15, 15))
    axes = axes.flatten()
    
    # Create a color map for clusters
    cluster_colors = plt.cm.get_cmap('tab10', num_clusters)
    
    for idx, (cluster, ax) in enumerate(zip(clusters, axes)):
        subgraph = G.subgraph(cluster)
        pos = nx.spring_layout(subgraph, seed=42 + idx)  # Different seed for each cluster for separation
        
        # Draw the subgraph
        nx.draw(subgraph, pos, with_labels=True, node_color=[cluster_colors(idx)] * len(cluster), 
                edge_color=[cluster_colors(idx)] * len(subgraph.edges), node_size=500, font_size=10, 
                font_color='white', ax=ax)
        
        ax.set_title(f'Cluster {idx + 1}')
        ax.axis('off')  # Turn off axis
    
    # Hide any unused subplots
    for j in range(num_clusters, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.savefig('edge_network.png')
    plt.close()


###### cluster emission functions ######

def save_clusters(transformed_matrix, labels, filename='cluster_assignments.csv'):
    df = pd.DataFrame({
        'MDS Dimension 1': transformed_matrix[:, 0],
        'MDS Dimension 2': transformed_matrix[:, 1],
        'Cluster Label': labels
    })

    df.to_csv(filename, index=False)

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

def main():
    #methods dict key = name of arg to trigger
    # value = tuple (method_function, [required_args], [optional_args (set by default)])
    methods_dict = {
    'kmeans': (kmeans_cluster, ['matrix'], []),
    'heirarchy': (hierarchy_cluster, ['matrix', 'accessions'], []),
    'hdbscan': (umap_clustering, ['matrix'], []),
    'edge_based': (edge_based_cluster, ['matrix', 'accessions'], ['minimum_edge']),
    }

    parser = argparse.ArgumentParser(description='subsample from a matrix')
    parser.add_argument("--phylip", 
    type=str, 
    help='Path to the Phylip file'
    ),
    
    parser.add_argument('--methods', 
        nargs='+',
        choices=list(methods_dict.keys()) + ['all'],  # only allow valid methods and all trigger
        required=True,
        help="method or methods to use for clustering"
    ),

    args = parser.parse_args()

    _, accessions, matrix = read_phylip_distance(args.phylip)

    context = {
    'matrix': matrix,
    'accessions': accessions,
    'minimum_edge': 0.0005,
    }


    if 'all' in args.methods:
        selected_methods = methods_dict.keys()
    else:
        selected_methods = args.methods

    for method_name in selected_methods:
        method, required_args, optional_args = methods_dict[method_name]

        method_args = [context[arg] for arg in required_args if arg in context]

        #unused currently but can take args from context where needed
        method_optional_args = [context[arg] for arg in optional_args if arg in context and context[arg] is not None]

        result = method(*method_args, *method_optional_args)  # Unpack the arguments

        if method_name in ['kmeans', 'heirarchy', 'hdbscan']:
            plot_umap(matrix, result, f'{method_name}_umap.png')
        else:
            relate_id_to_accession(result, accessions, "clusters.txt")


    #reps = select_closest_representatives(matrix, kmeans_labels, accessions, n_representatives=3)

if __name__ == "__main__":
    main()
