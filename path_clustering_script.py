import json
import os
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import euclidean_distances
import time

def create_id_to_node_dict(nodes):
    id_to_node = {}
    for node in nodes:
        node_id = node["id"]
        id_to_node[node_id] = node
    return id_to_node


def write_clusters_to_pdb(paths, labels, nodes, output_dir):
    node_dict = create_id_to_node_dict(nodes)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    clusters = {}
    for path, label in zip(paths, labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(path)

    for cluster_label, paths in clusters.items():
        filename = os.path.join(output_dir, f'cluster_{cluster_label}.pdb')
        with open(filename, 'w') as f:
            f.write('HEADER Cluster {}\n'.format(cluster_label))
            chain_id = 'A'
            atom_id = 1
            for path in paths:
                if not path:
                    continue
                for node_id in path:
                    node = node_dict[node_id]
                    f.write('ATOM  {:>5}  C   ALA 1111    {:>8.3f}{:>8.3f}{:>8.3f}  1.00  0.00           C\n'.format(atom_id, node['x'], node['y'], node['z']))
                    atom_id += 1
                f.write('TER\n')
                chain_id = chr(ord(chain_id) + 1)
            f.write('END\n')


def write_clusters_to_pdb_alt(paths, labels, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    clusters = {}
    for path, label in zip(paths, labels):
        if label not in clusters:
            clusters[label] = []
        clusters[label].append(path)

    for cluster_label, paths in clusters.items():
        filename = os.path.join(output_dir, f'cluster_{cluster_label}.pdb')
        with open(filename, 'w') as f:
            f.write('HEADER Cluster {}\n'.format(cluster_label))
            chain_id = '!'
            atom_id = 1
            for path in paths:
                for node_id in path:
                    x, y, z, radius = node_id
                    f.write('ATOM  {:>5}  C   ALA 1111    {:>8.3f}{:>8.3f}{:>8.3f}  1.00{:>5.2f}           C\n'.format(atom_id, x, y, z, radius))
                    atom_id += 1
                f.write('TER\n')
                chain_id = chr(ord(chain_id) + 1)
            f.write('END\n')


def path_similarity(path1, path2):
    path1, path2 = np.asarray(path1), np.asarray(path2)
    distance_matrix = euclidean_distances(path1, path2)
    # Compute the directed Hausdorff distances between the two paths
    hausdorff_distance_1 = np.max(np.min(distance_matrix, axis=1))
    hausdorff_distance_2 = np.max(np.min(distance_matrix, axis=0))
    # Return the maximum of the two directed Hausdorff distances
    return max(hausdorff_distance_1, hausdorff_distance_2)


def build_similarity_matrix(matrix, paths):
    for i in range(len(matrix)):
        for j in range(i, len(matrix)):
            if(i == j):
                matrix[i][j] = 0
            else:
                path1 = paths[i]
                path2 = paths[j]

                distance = path_similarity(path1, path2)

                matrix[i][j] = distance
                matrix[j][i] = distance


def get_representative_path(paths_as_nodes, labels):
    representative_paths = []

    unique_labels = np.unique(labels)
    for label in unique_labels:
        if label == -1:
            paths = [path for i, path in enumerate(paths_as_nodes) if labels[i] == label]
        else:
            paths = [path for i, path in enumerate(paths_as_nodes) if labels[i] == label]
            distances = np.zeros((len(paths), len(paths)))
            for i in range(len(paths)):
                for j in range(len(paths)):
                    distances[i][j] = path_similarity(paths[i], paths[j])
            medoid_index = np.argmin(distances.sum(axis=0))
            medoid_path = paths[medoid_index]
            representative_paths.append(medoid_path)
    print("Paths after clustering: ", len(representative_paths))
    return representative_paths


def write_representative_path_to_pdb(representative_paths, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cluster_label = 0
    for path in representative_paths:
        filename = os.path.join(output_dir, f'cluster_{cluster_label}_medioid.pdb')
        with open(filename, 'w') as f:
            f.write('HEADER Cluster {}\n'.format(cluster_label))
            atom_id = 1
            for node_id in path:
                x, y, z, radius = node_id
                f.write('ATOM  {:>5}  C   ALA 1111    {:>8.3f}{:>8.3f}{:>8.3f}  1.00{:>5.2f}           C\n'.format(atom_id, x, y, z, radius))
                atom_id += 1
            f.write('TER\n')
        cluster_label+=1


with open('output_alt.json', 'r') as f:
    data = json.load(f)

snapshots = data['snapshots']
paths_from_all_snapshots = []

start_total = time.time()

for snapshot in snapshots:
    snapshot_id = snapshot['snapshot_id']
    print("Snapshot: ", snapshot_id)
    nodes = snapshot['voronoi_nodes']
    paths = snapshot['paths']

    paths_as_nodes = []
    if (not len(paths)==0):
        if(len(paths)<5000):
            step = 1
            start_total_snapshot = time.time()
            start_building_dict = time.time()
            for path in paths:
                step+=1
                path_nodes = []
                for node_id in path:
                    node = next((n for n in nodes if n["id"] == node_id), None)
                    if node:
                        path_nodes.append([node["x"], node["y"], node["z"], node["radius"]])
                paths_as_nodes.append(np.array(path_nodes))
            end_building_dict = time.time()

            print("Paths before clustering: ", len(paths_as_nodes))
            print("Building path dictionary: ", end_building_dict - start_building_dict)
            start_snapshot_build_matrix = time.time()

            similarity_matrix = np.zeros(shape=(len(paths_as_nodes), len(paths_as_nodes)))
            build_similarity_matrix(similarity_matrix, paths_as_nodes)
            end_snapshot_build_matrix = time.time()
            snapshot_build_matrix_elapsed_time = end_snapshot_build_matrix - start_snapshot_build_matrix
            print("Building snapshot matrix: ", snapshot_build_matrix_elapsed_time)

            start_snapshot_clustering = time.time()
            dbscan = DBSCAN(eps=5, min_samples=5, metric='precomputed').fit(similarity_matrix)
            end_snapshot_clustering = time.time()
            snapshot_clustering_elapsed_time = end_snapshot_clustering - start_snapshot_clustering
            print("Clustering snapshot: ", snapshot_clustering_elapsed_time)

            labels = dbscan.labels_
            start_get_medioid = time.time()
            paths_from_all_snapshots.extend(get_representative_path(paths_as_nodes, labels))
            end_get_medioid = time.time()
            get_medioid_elapsed_time = end_get_medioid - start_get_medioid
            print("Calculating medioid and extending common matrix: ", get_medioid_elapsed_time)
            end_total_snapshot = time.time()
            print('Processing snapshot ', snapshot_id, ' calc time: ', end_total_snapshot - start_total_snapshot)
    print("-----------------------------")

print("Finished clustering each snapshot. Total paths: ", len(paths_from_all_snapshots))

start_final_clustering_build_matrix = time.time()

similarity_matrix_all = np.zeros(shape=(len(paths_from_all_snapshots), len(paths_from_all_snapshots)))
build_similarity_matrix(similarity_matrix_all, paths_from_all_snapshots)

end_final_clustering_build_matrix = time.time()
final_clustering_build_matrix_elapsed_time = end_final_clustering_build_matrix - start_final_clustering_build_matrix
print("Time of building matrix final clustering: ", final_clustering_build_matrix_elapsed_time)

start_final_clustering = time.time()
dbscan = DBSCAN(eps=5, min_samples=2, metric='precomputed').fit(similarity_matrix_all)
end_final_clustering = time.time()
final_clustering_elapsed_time = end_final_clustering - start_final_clustering
print("Time of final clustering: ", final_clustering_elapsed_time)

labels = dbscan.labels_
print("save to pdb")
write_clusters_to_pdb_alt(paths_from_all_snapshots, labels, "output_directory_clustering")
end_total = time.time()
total_elapsed_time = end_total - start_total
print('Total processing time: ',total_elapsed_time)

final_clustering_single_representatives = get_representative_path(paths_from_all_snapshots, labels)
print("represent: ", len(final_clustering_single_representatives))
write_representative_path_to_pdb(final_clustering_single_representatives,"output_single_medioids")

print("Results saved")
