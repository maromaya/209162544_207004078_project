import sys
import pandas as pd
import numpy as np
import symnmfmodule as symnmf


def print1(lst):
    for vector in lst:
        str_to_print = ""
        for num in vector:
            str_to_print += str(format(num,'.4f')) +','
        print(str_to_print[:-1])


def my_main(argv):

    if len(argv) == 3:
        try:
            k = int(argv[0])
            goal = argv[1]
            file_path = argv[2]
            data = pd.read_csv(file_path, header= None, sep=',')
            dataP_lst = [list(row) for row in data.values]
            n = len(dataP_lst)
            if k >= n:
                print("An Error Has Occurred")
                return
            start(k, goal, dataP_lst)
        except Exception as ex:
            print("An Error Has Occurred")
            return

    else:
        print("An Error Has Occurred")
        return

def start_symnmf (dataP_lst,k):
    W = symnmf.fit_sym_ddg_norm(dataP_lst, len(dataP_lst), 3, len(dataP_lst[0]))
    H_rand = init_H(W, k)
    res = symnmf.fit_symnmf(W, H_rand, len(W), len(H_rand), k)

    return res

    
def start(k, goal, dataP_lst):

    if(goal =="sym"):
        res = symnmf.fit_sym_ddg_norm(dataP_lst, len(dataP_lst), 1, len(dataP_lst[0]))
        print1(res)
    elif (goal == "ddg"):
        res = symnmf.fit_sym_ddg_norm(dataP_lst, len(dataP_lst), 2,  len(dataP_lst[0]))
        print1(res)
    elif (goal == "norm"):
        res = symnmf.fit_sym_ddg_norm(dataP_lst, len(dataP_lst), 3,  len(dataP_lst[0]))
        print1(res)
    elif (goal == "symnmf"):
        W = symnmf.fit_sym_ddg_norm(dataP_lst, len(dataP_lst), 3, len(dataP_lst[0]))
        H_rand = init_H(W, k)
        res = symnmf.fit_symnmf(W, H_rand, len(W), len(H_rand), k)
        print1(res)
    else:
        print("An Error Has Occurred")


def init_H(W_matrix, k):
    arr = np.array(W_matrix)
    avg = np.mean(arr)
    np.random.seed(0)
    upper_limit = 2 * (avg/k) ** 0.5
    H = [[np.random.uniform(0,upper_limit) for _ in range(k)] for _ in range(len(W_matrix))]
    return H


if __name__ == '__main__':
    my_main(sys.argv[1:])

