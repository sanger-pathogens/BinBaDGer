#!/usr/bin/env python3

import argparse
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import umap
from collections import defaultdict

from sklearn.manifold import MDS
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from networkx.algorithms.community import greedy_modularity_communities

import imageio
import os

##### methods handling class #####

class ClusteringMethods:
    def __init__(self, **kwargs):
        self.context = kwargs

        self.methods = {
            'kmeans': (kmeans_cluster, ['matrix'], []),
            'heirarchy': (hierarchy_cluster, ['matrix', 'accessions'], []),
            'hdbscan': (umap_clustering, ['matrix'], []),
            'edge_based': (edge_based_cluster, ['matrix', 'accessions'], ['minimum_edge']),
            'network_based': (trim_network_to_n_nodes, ['matrix', 'accessions'], []),
        }

    def run_method(self, method_name):
        method, required_args, optional_args = self.methods[method_name]
        
        missing_args = [arg for arg in required_args if arg not in self.context]
        if missing_args:
            raise ValueError(f"Missing required arguments for '{method_name}': {', '.join(missing_args)}")

        method_args = [self.context[arg] for arg in required_args]
        method_optional_args = [self.context[arg] for arg in optional_args if arg in self.context]

        return method(*method_args, *method_optional_args)

##### input parsing functions #####

def read_phylip_distance(phylip_file):
    with open(phylip_file, 'r') as phylip:
        lines = phylip.readlines()

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

def edge_based_cluster(matrix, accessions, threshold, N=3, dissimilarity=True):
    G = nx.Graph()
    
    num_nodes = len(matrix)

    accession_map = {i: accessions[i] for i in range(num_nodes)}
    
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if matrix[i, j] < threshold:
                G.add_edge(accession_map[i], accession_map[j], weight=matrix[i, j])

    #as singletons don't get an edge above just save them at this stage
    all_nodes = set(accession_map.values())
    non_singletons = set(G.nodes())
    singletons = all_nodes - non_singletons
    
    clusters = greedy_modularity_communities(G)
    
    if not clusters and not singletons: #don't expect this really
        raise Exception('No samples identified')

    combined_clusters = [list(cluster) for cluster in clusters]
    
    for singleton in singletons:
        combined_clusters.append([singleton])
    
    representatives = {}
    
    # Loop through each cluster
    for cluster in combined_clusters:
        cluster = list(cluster)
        
        if len(cluster) >= N + 1:
            # Compute degree centrality for large clusters
            centrality = nx.degree_centrality(G.subgraph(cluster))
            sorted_nodes = sorted(centrality, key=centrality.get, reverse=dissimilarity)
            representatives[tuple(cluster)] = sorted_nodes[:N]
        else:
            # For smaller clusters, store all members
            representatives[tuple(cluster)] = cluster
    
    plot_network_subclusters(representatives.keys(), G, representatives)
    
    return representatives

### Network trimming ###

def trim_network_to_n_nodes(matrix, accessions, N=10, plot_seed=123, plot_iterations=True):
    G = nx.Graph()

    num_nodes = len(matrix)

    accession_map = {i: accessions[i] for i in range(num_nodes)}
    
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            G.add_edge(accession_map[i], accession_map[j], weight=matrix[i, j])

    if len(G.nodes) <= N:
        return {tuple(G.nodes): list(G.nodes)} #if there isn't enough to do the selection just return the nodes in the same structure as below

    # Sort edges by weight (shortest first)
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]['weight'])
    
    trimmed_graph = G.copy()
    iteration = 0
    filenames = []
    
    while len(trimmed_graph.nodes) > N:
        if plot_iterations:
            filename = plot_current_graph(trimmed_graph, iteration, plot_seed)
            filenames.append(filename)
        
        # Get the shortest edge and remove one of its nodes
        shortest_edge = sorted_edges.pop(0)
        node_to_remove = shortest_edge[0] if trimmed_graph.degree(shortest_edge[0]) <= trimmed_graph.degree(shortest_edge[1]) else shortest_edge[1]
        
        # Remove the chosen node and its edges
        trimmed_graph.remove_node(node_to_remove)
        
        # Remove edges related to the removed node from the sorted edge list
        sorted_edges = [edge for edge in sorted_edges if node_to_remove not in edge[:2]]
        
        iteration += 1
    
    plot_current_graph(trimmed_graph, iteration, plot_seed)

    if plot_iterations:
        # Final plot
        filename = plot_current_graph(trimmed_graph, iteration, plot_seed)
        filenames.append(filename)
        create_gif(filenames)
        # Clean up image files as we are saving gif
        for filename in filenames:
            os.remove(filename)

    clusters = list(nx.connected_components(trimmed_graph))
    representatives = {tuple(cluster): list(cluster) for cluster in clusters}
    
    return representatives

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

def plot_network_subclusters(clusters, G, representatives=None, plot_seed=123):
    num_clusters = len(clusters)
    
    # Determine grid size for subplots
    grid_size = int(np.ceil(np.sqrt(num_clusters)))
    
    # Create subplots
    fig, axes = plt.subplots(grid_size, grid_size, figsize=(15, 15))
    axes = np.array(axes).flatten()
    
    # Create a color map for clusters
    cluster_colors = plt.cm.get_cmap('tab10', num_clusters)
    
    for idx, (cluster, ax) in enumerate(zip(clusters, axes)):
        subgraph = G.subgraph(cluster)
        pos = nx.spring_layout(subgraph, seed=plot_seed) # Might want to change the seed sometimes?
        
        node_colors = []
        node_sizes = []
        
        for node in cluster:
            if representatives and node in representatives[tuple(cluster)]:
                node_colors.append('red')  # representatives in red
                node_sizes.append(1000)
            else:
                node_colors.append(cluster_colors(idx))
                node_sizes.append(500)
        
        # Draw the subgraph
        nx.draw(subgraph, pos, with_labels=True, node_color=node_colors, 
                edge_color=[cluster_colors(idx)] * len(subgraph.edges), 
                node_size=node_sizes, font_size=10, font_color='black', ax=ax)
        
        ax.set_title(f'Cluster {idx + 1}')
        ax.axis('off')
    
    # Hide any unused subplots
    for j in range(num_clusters, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()
    plt.savefig('edge_network_with_highlighted_representatives.png')
    plt.close()

def plot_current_graph(G, iteration, plot_seed):
    pos = nx.spring_layout(G, seed=plot_seed)  # Layout for consistent graph drawing
    plt.figure(figsize=(8, 8))
    
    nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', 
            node_size=700, font_size=10, font_color='black')

    filename=(f'network_iteration_{iteration}.png')
    
    plt.title(f"Iteration {iteration}: {len(G.nodes)} nodes remaining")
    plt.savefig(filename)
    plt.close()

    return filename #we getting fancy


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

def save_representatives_to_file(combined_clusters, output_file="representatives.txt"):
    with open(output_file, 'w') as file:
        for cluster, members in combined_clusters.items():
            for member in members:
                file.write(f"{member}\n")

def relate_id_to_accession(clusters, accessions, output_file="representatives.txt"):
    cluster_dict = defaultdict(list)
    for idx, cluster in enumerate(clusters):
        cluster_dict[cluster].append(accessions[idx])
    
    with open(output_file, 'w') as out:
        for cluster, members in cluster_dict.items():
            out.write(f"cluster {cluster}: {', '.join(members)}\n")

def create_gif(filenames, gif_filename='trimming_process.gif', duration=0.5):
    """Creates a GIF from a list of image filenames."""
    with imageio.get_writer(gif_filename, mode='I', duration=duration) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)

def main():
    parser = argparse.ArgumentParser(description='subsample from a matrix')
    parser.add_argument("--phylip", 
    type=str,
    required=True,
    help='Path to the Phylip file'
    ),
    
    parser.add_argument('--methods', 
        nargs='+',
        choices=['kmeans', 'heirarchy', 'hdbscan', 'edge_based', 'network_based', 'all'],
        required=True,
        help="method or methods to use for clustering"
    ),

    parser.add_argument('--minimum_edge', 
        type=float,
        default=0.01,
        help="minimum_edge for bringing forward to network"
    ),

    args = parser.parse_args()

    _, accessions, matrix = read_phylip_distance(args.phylip)

    clustering_methods = ClusteringMethods(
        matrix=matrix,
        accessions=accessions,
        minimum_edge=args.minimum_edge
    )

    selected_methods = clustering_methods.methods.keys() if 'all' in args.methods else args.methods

    for method_name in selected_methods:
        result = clustering_methods.run_method(method_name)

        if method_name in ['kmeans', 'heirarchy', 'hdbscan']:
            plot_umap(matrix, result, f'{method_name}_umap.png')
        else:
            save_representatives_to_file(result)


if __name__ == "__main__":
    main()