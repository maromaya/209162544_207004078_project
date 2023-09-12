# define PY_SSIZE_T_CLEAN
#include "symnmf.h"
#include <Python.h>

/* MUST include <Python.h>, this implies inclusion of the following standard headers:
     *  <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
/* include <Python.h> has to be before any standard headers are included */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static PyObject* myProject (PyObject *self, PyObject *args) {
    PyObject  *vectors;
    PyObject *vector;
    PyObject *number;
    PyObject *py_matrix;
    PyObject *py_vector;
    int temp;
    int n;
    int n_goal;
    double num;

    if (!PyArg_ParseTuple(args, "Oiii", &vectors,&n,&n_goal,&temp)) {
        return NULL;
    }

    double **all_vectors = malloc(n * sizeof(double *));
    if (all_vectors == NULL) {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    int e;
    for (e = 0; e < n; e++) {
        all_vectors[e] = calloc(temp,sizeof(double));
    }
    int i,j;

    for (i = 0; i < n; i++) {
        vector =  PyList_GetItem(vectors, i);

        if (!PyList_Check(vector)) {
            printf("Invalid vector object at index %d. Exiting.\n", i);
            return NULL;
        }
        for (j = 0; j<temp;j++){
            number =  PyList_GetItem(vector, j);

            if (!PyFloat_Check(number)) {
                printf("Invalid number object at index (%d, %d). Exiting.\n", i, j);
                return NULL;
            }
            num = PyFloat_AsDouble(number);
            all_vectors[i][j] = num;


        }
    }




    int g;
    double **result = malloc(n * sizeof(double *));
    for (g = 0; g < n; g++) {
        result[g] = calloc(n, sizeof(double));
    }
    if (n_goal ==1){
        result= sym(all_vectors,n, temp);
    }
    if(n_goal==2){
        result= ddg(all_vectors,n, temp);
    }
    if(n_goal==3){
        result= norm(all_vectors,n, temp);

    }
    // Create a Python list to hold the results
    py_matrix= PyList_New(n);

// Convert the results array to Python lists
    for (g = 0; g < n; g++) {
        py_vector = PyList_New(n);
        for (int w = 0; w < n; w++) {
            PyList_SetItem(py_vector, w, PyFloat_FromDouble(result[g][w]));
        }
        PyList_SetItem(py_matrix, g, py_vector);
    }

    int o;
    for (o= 0; o< n; o++){
        free(all_vectors[o]);
        free(result[o]);

    }
    free(all_vectors);
    free(result);
    return py_matrix;
}

static PyObject* my_symnmf (PyObject *self, PyObject *args) {
    PyObject  *vectors_W;
    PyObject  *vectors_H;
    PyObject *vector_W;
    PyObject *vector_H;
    PyObject *number_W;
    PyObject *number_H;
    PyObject *py_matrix;
    PyObject *py_vector;
    int W_num_of_rows;
    int H_num_of_rows;
    int H_num_of_cols;
    double num_W;
    double num_H;
    if (!PyArg_ParseTuple(args, "OOiii", &vectors_W,&vectors_H,&W_num_of_rows,&H_num_of_rows,& H_num_of_cols)) {
        return NULL;
    }

    double **all_vectors_W = malloc(W_num_of_rows * sizeof(double *));
    if (all_vectors_W == NULL) {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    int e;
    for (e = 0; e < W_num_of_rows; e++) {
        all_vectors_W[e] = calloc(W_num_of_rows,sizeof(double));
    }

    double **all_vectors_H = malloc(H_num_of_rows* sizeof(double *));
    if (all_vectors_H== NULL) {
        printf("Memory allocation failed. Exiting.\n");
        return NULL;
    }
    for (e = 0; e < H_num_of_rows; e++) {
        all_vectors_H[e] = calloc(H_num_of_cols,sizeof(double));
    }
    int i,j;

    for (i = 0; i < W_num_of_rows; i++) {
        vector_W =  PyList_GetItem(vectors_W, i);

        if (!PyList_Check(vector_W)) {
            printf("Invalid vector object at index %d. Exiting.\n", i);
            return NULL;
        }
        for (j = 0; j<W_num_of_rows;j++){
            number_W =  PyList_GetItem(vector_W, j);

            if (!PyFloat_Check(number_W)) {
                printf("Invalid num object at index (%d, %d). Exiting.\n", i, j);
                return NULL;
            }
            num_W = PyFloat_AsDouble(number_W);
            all_vectors_W[i][j] = num_W;


        }
    }

    for (i = 0; i < H_num_of_rows; i++) {
        vector_H =  PyList_GetItem(vectors_H, i);

        if (!PyList_Check(vector_H)) {
            printf("Invalid vector object at index %d. Exiting.\n", i);
            return NULL;
        }
        for (j = 0; j<H_num_of_cols;j++){
            number_H =  PyList_GetItem(vector_H, j);

            if (!PyFloat_Check(number_H)) {
                printf("Invalid num object at index (%d, %d). Exiting.\n", i, j);
                return NULL;
            }
            num_H = PyFloat_AsDouble(number_H);
            all_vectors_H[i][j] = num_H;


        }
    }


    int g;
    double **result = malloc(H_num_of_rows * sizeof(double *));
    for (g = 0; g < H_num_of_rows; g++) {
        result[g] = calloc(H_num_of_cols, sizeof(double));
    }
 

    result = symnmf(all_vectors_H,all_vectors_W,W_num_of_rows,H_num_of_rows,H_num_of_cols);


    // Create a Python list to hold the results
    py_matrix = PyList_New(H_num_of_rows);


    // Convert the results array to Python lists
    for (g = 0; g < H_num_of_rows; g++) {
    py_vector = PyList_New(H_num_of_cols);
    for (int w = 0; w < H_num_of_cols; w++) {
    PyList_SetItem(py_vector, w, PyFloat_FromDouble(result[g][w]));
    }
    PyList_SetItem(py_matrix, g, py_vector);
    }

    int o;
    for (o= 0; o< W_num_of_rows; o++){
    free(all_vectors_W[o]);
    }
    for (o= 0; o< H_num_of_rows; o++){
    free(all_vectors_H[o]);
    free(result[o]);

    }
    free(all_vectors_H);
    free(all_vectors_W);
    free(result);

    return py_matrix;
}
static PyMethodDef My_Methods[] = { 
        {
                "fit_sym_ddg_norm" ,     
                (PyCFunction) myProject, 
                            METH_VARARGS,         
                PyDoc_STR("return sym ddg or norm matrix") 
        },

        {
                "fit_symnmf" ,      
                (PyCFunction) my_symnmf, 
                            METH_VARARGS,         
                PyDoc_STR("return final H") 
        },
        {
                NULL, NULL, 0, NULL
        }
};

static struct PyModuleDef Moduledef = {

        PyModuleDef_HEAD_INIT,
        "symnmfmodule",     
        "module for somthing that you told us to do" , 
        -1,
        My_Methods
};


PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    return PyModule_Create(&Moduledef);
}