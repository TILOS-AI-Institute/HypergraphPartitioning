import sys
import numpy as np
from sklearn import metrics
from sklearn.metrics import pairwise_distances
from sklearn.metrics import davies_bouldin_score

def read_files(embedding_file, labels_file):
    labels = np.loadtxt(labels_file, dtype=np.int)
    embedding = np.loadtxt(embedding_file, dtype=np.float32)
    return embedding, labels

def db_index(embedding, labels):
    return davies_bouldin_score(embedding, labels)

def vrc_index(embedding, labels):
    return metrics.calinski_harabasz_score(embedding, labels)

def sc_index(embedding, labels):
    return metrics.silhouette_score(embedding, labels)

embedding_file = sys.argv[1]
labels_file = sys.argv[2]

embedding, labels = read_files(embedding_file, labels_file)
dbi_score = db_index(embedding, labels)
print("[cluster-metric] Davies-Bouldin score ", dbi_score)
vrc_score = vrc_index(embedding, labels)
print("[cluster-metric] VRC/Calinski-Harabasz score ", vrc_score)
sc_score = sc_index(embedding, labels)
print("[cluster-metric] Silhouette Coefficient score ", sc_score)