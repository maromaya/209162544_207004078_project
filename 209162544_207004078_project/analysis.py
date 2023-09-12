import sys
import symnmf
import math
import pandas as pd
import numpy as np
from sklearn.metrics import silhouette_score



def my_main(argv):
    if len(argv) == 2:
        try:
            k = int(argv[0])
            file_path = argv[1]
            data = pd.read_csv(file_path, header= None, sep=',')
            dataP_lst = [list(row) for row in data.values]
            n = len(dataP_lst)
            if k >= n:
                print("An Error Has Occurred")
                return
            start_analysis(k,dataP_lst)
        except Exception as ex:
            print("An Error Has Occurred")
            return

    else:
        print("An Error Has Occurred")
        return


def symnmf_clus(matrix):
    np_mat = np.array(matrix)
    cluster_elem = np.argmax(np_mat, axis=1)
    return cluster_elem


def kmeans_clus(datapoints, k_mean_matrix):
    kmean_lst = []
    for elem in k_mean_matrix:
        kmean_lst.append(elem[0])
    kmean_mat_np = np.array(kmean_lst)
    cluster_elem = []
    for row in datapoints:
        distances = np.linalg.norm(kmean_mat_np - row, axis=1)  
        index = np.argmin(distances)  
        cluster_elem.append(index)
    return cluster_elem


def ED(lst1, lst2):
    dis = 0
    for i in range(len(lst1)):
        dis += math.pow(lst1[i]-lst2[i], 2)
    return math.sqrt(dis)



def conv_check(lst1, lst2):
    for i in range(len(lst1)):
        if 0.00001 < ED(lst1[i][0], lst2[i][0]):
            return False
    return True

def init_kmeans(distance):
    lst_vectors = []
    for i in range(len(distance)):
        lst = []
        kmeans_lst = calc_kmeans_lst(distance[i])
        lst.append(kmeans_lst)
        lst_vectors.append(lst)
    return lst_vectors


def start_analysis(k,dataP_lst):
    H = symnmf.start_symnmf(dataP_lst,k)
    kmeans = kmeans_func(k,dataP_lst,300)
    cluster_H_index = symnmf_clus(H)
    cluster_kmeans_index = kmeans_clus(dataP_lst,kmeans)
    kmeans_score = silhouette_score(dataP_lst, cluster_kmeans_index)
    H_score = silhouette_score(dataP_lst, cluster_H_index)
    print(f"nmf: {H_score:.4f}")
    print(f"kmeans: {kmeans_score:.4f}")




def calc_kmeans_lst(list_of_ki_including_the_first_param):
    kmeans_lst = []
    lst_of_vectors = list_of_ki_including_the_first_param[1:]
    for i in range(len(lst_of_vectors[0])):
        average = 0
        for j in range(len(lst_of_vectors)):
            average += lst_of_vectors[j][i]
        kmeans_lst.append(float(format(round(average/len(lst_of_vectors),4),'.4f')))
    return kmeans_lst


def init(k,vectors):
    lst_vectors = []
    for i in range(k):
        lst = []
        lst.append(vectors[i])
        lst_vectors.append(lst)
    return lst_vectors


def add_new_vec(distance, vectors):
    for vector in vectors:
        place = find_k(vector, distance)
        distance[place].append(vector)


def find_k(vector, K_lst):
    min = math.inf
    place = 0
    for i in range(len(K_lst)):
        temp = ED(vector, K_lst[i][0])
        if temp < min:
            min = temp
            place = i
    return place


def kmeans_func(k, data, iter1 = 300):
    vectors = data
    distance = init(k,vectors)
    converge = False
    i = 0
    while not converge and i < iter1:
        add_new_vec(distance, vectors)
        newDist = init_kmeans(distance)
        converge = conv_check(distance, newDist)
        distance = newDist
        i += 1
    return distance


if __name__ == "__main__":
    my_main(sys.argv[1:])
