#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"
int num_of_vectors;
int vector_size;
double b = 0.5;
double eps = 1e-4;
int iter = 300;



double **load_data(const char* path) {
    double **data = calloc(1, sizeof(double *));
    int k = 0;
    size_t read;
    size_t k_one = -1;
    char *line = NULL;
    size_t row_len = 0;
    FILE *file;
    char *ptr;
    double val;
    size_t count;


    file = fopen(path, "r");
    if (file == NULL) {
        printf("Error opening file");
        exit(1);
    }
    read = (size_t)getline(&line, &row_len, file);
    while (read != k_one) {
        double *row = calloc(1, sizeof(double));
        row_len = 0;
        count = 0;
        ptr = line;
        while (*ptr != '\n') {
            int scanned = sscanf(ptr, "%lf", &val);
            if(scanned==-1){
                break;
            }
            if (scanned == 1) {
                if (count >= row_len) {
                    row_len = count + 1;
                    row = realloc(row, row_len * sizeof(double));
                    if (row == NULL) {
                        free(row);
                    }
                }
                row[count] = val;
                count++;
            }


            while ((*ptr >= '0' && *ptr <= '9') || *ptr == '.' || *ptr == '-') {
                ptr++;
            }
            if (*ptr == ',') {
                ptr++;
            }
        }


        if (count > 0) {
            data = realloc(data, (k + 1) * sizeof(double *));
            data[k] = row;
            k++;


        } else {
            free(row);
        }
    }
    free(line);
    num_of_vectors = k;
    vector_size = row_len;
    fclose(file);
    return data;
}




double** sym (double** data_points, int rowNum, int row_size){
    double **sym_matrix;
    int i,k,j;
    double new_ed;
    sym_matrix = malloc(rowNum * sizeof(double *));
    for (i = 0;i<rowNum;i++) {
        sym_matrix[i] = calloc(rowNum, sizeof(double));
    }
    for(j = 0 ;j<rowNum;j++){
        for( k=0;k<rowNum;k++){
            if(j!=k) {
                if (j < k) {
                double ed = ED(data_points[j], data_points[k], row_size);
                new_ed = exp((ed / 2) * -1);
                sym_matrix[j][k] = new_ed;
                sym_matrix[k][j] = new_ed;
              }

            }
        }
    }
    return sym_matrix;
}

double ED(double *lst1,double *lst2,int n) {
    double distance = 0;
    int i;
    for (i = 0; i < n; i++) {
        distance += pow(lst1[i] - lst2[i], 2);
    }

    return distance;
}

double sum_row(double *lst1,int n) {
    double sum = 0;
    int i;
    
    for (i = 0; i < n; i++) {
        sum += lst1[i];
    }

    return sum;
}


double sum_vectors(double *lst1,double *lst2,int n) {
    double sum = 0;
    int i;
    for (i = 0; i < n; i++) {
        sum += lst1[i] * lst2[i];
    }
    return sum;
}

double** ddg (double ** datapoints, int rowNum, int row_size ){
    double **result;
    double **sym_ddg_mat;
    int i;
    sym_ddg_mat = sym(datapoints,rowNum,row_size);
    result = malloc(rowNum * sizeof(double *));
    for (i = 0;i<rowNum;i++) {
        result[i] = calloc(rowNum, sizeof(double));
    }
    for (i = 0;i<rowNum;i++) {
        result[i][i] = sum_row(sym_ddg_mat[i],rowNum);

    }


    for (i= 0; i< rowNum; i++){

        free(sym_ddg_mat[i]);
    }

    free(sym_ddg_mat);
    return result;
}

double** mul_diag_norm_left(double** matrix1, double** matrix2, int n){
    int i,j;
    double **res;

    res = malloc(n * sizeof(double *));
    for (i = 0;i<n;i++) {
        res[i] = calloc(n, sizeof(double));
    }
    for (i=0;i<n;i++){
        double x = pow(matrix1[i][i],-0.5);
        for(j=0;j<n;j++){
            double temp = x*matrix2[i][j];
            res[i][j] = temp;
        }
    }
    return res;
}


double** mul_diag_norm_right(double** matrix1, double** matrix2, int n){
    int i,j;
    double **res;

    res = malloc(n * sizeof(double *));
    for (i = 0;i<n;i++) {
        res[i] = calloc(n, sizeof(double));
    }
    for (i=0;i<n;i++){
        double x = pow(matrix1[i][i],-0.5);
        for(j=0;j<n;j++){
            double temp = x*matrix2[j][i];
            res[i][j] = temp;
        }
    }
    return res;
}


double** norm (double** data_points, int rowNum, int row_size ){
    double **sym_matrix_norm;
    double **dig_matrix_norm;
    double **norm_matrix;
    double **DA_norm_matrix;
    int i;

    sym_matrix_norm = sym(data_points,rowNum,row_size);
    dig_matrix_norm = ddg(data_points,rowNum,row_size);

    DA_norm_matrix =mul_diag_norm_left(dig_matrix_norm, sym_matrix_norm, rowNum);
    norm_matrix =mul_diag_norm_right(dig_matrix_norm, DA_norm_matrix, rowNum);
    for (i= 0; i< rowNum; i++){

        free(sym_matrix_norm[i]);
        free(dig_matrix_norm[i]);
        free(DA_norm_matrix[i]);
    }

    free(sym_matrix_norm);
    free(dig_matrix_norm);
    free(DA_norm_matrix);

    return norm_matrix;
}



double** transpose(double ** matrix, int rowNum, int colNum){
    double** transposed;
    int i,j,k;
    transposed = malloc(colNum * sizeof(double *));
    for (i = 0;i<colNum;i++) {
        transposed[i] = calloc(rowNum, sizeof(double));
    }
    for (j = 0; j<rowNum; j++){
        for (k = 0; k<colNum; k++){
            transposed[k][j] = matrix[j][k];
        }
    }
    return transposed;
}

double** updateH(double** W, double**H, int rowsW, int rowsH, int colsH){
    double** matrix;
    double** WH;
    double** transH;
    double** H_transH_H;
    double** H_transH;
    double** divMat;
    double** divMat_new;
    int i,j;
    matrix = malloc(rowsH * sizeof(double *));
    for (i = 0;i<rowsH;i++) {
        matrix[i] = calloc(colsH, sizeof(double));
    }
    WH = multi_matrix(W, H, rowsW, rowsW, colsH);
    transH =transpose(H, rowsH, colsH);
    H_transH = multi_matrix(H, transH, rowsH, colsH, rowsH);
    H_transH_H = multi_matrix(H_transH, H, rowsH, rowsH, colsH);

    divMat =divide_mat(WH, H_transH_H, rowsW, colsH);
    divMat_new =calcB(divMat, rowsW, colsH);
    matrix = mul_onebyone(H, divMat_new,rowsH, colsH);
    for (i= 0; i< rowsW; i++){
        free(WH[i]);
        free(divMat[i]);
        free(divMat_new[i]);
        free(H_transH_H[i]);
        free(H_transH[i]);
    }
    for(j=0;j<colsH; j++){
        free(transH[j]);
    }
    free(WH);
    free(divMat);
    free(divMat_new);
    free(H_transH_H);
    free(H_transH);
    free(transH);
    return matrix;
}

double converge(double** cur_H, double** old_H, int rowNum, int colNum){
    double conv_val = 0;
    int k,j;
    for (j = 0; j<rowNum; j++){
        for (k = 0; k<colNum; k++){
            conv_val += pow(cur_H[j][k] - old_H[j][k], 2);
        }
    }
    return conv_val;
}

double** calcB(double** matrix, int rowNum, int colNum){
    double** res;
    int i,j,k;

    res = malloc(rowNum * sizeof(double *));
    for (i = 0;i<rowNum;i++) {
        res[i] = calloc(colNum, sizeof(double));
    }
    for (j = 0; j<rowNum; j++){
        for (k = 0; k<colNum; k++){
            res[j][k] = (1-b) + (b * matrix[j][k]);
        }
    }
    return res;
}

double** divide_mat(double** matrix1, double** matrix2, int rowNum, int colNum){
    double** res;
    int i,j,k;
    res = malloc(rowNum * sizeof(double *));
    for (i = 0;i<rowNum;i++) {
        res[i] = calloc(colNum, sizeof(double));
    }
    for (j = 0; j<rowNum; j++){
        for (k = 0; k<colNum; k++){
            res[j][k] = matrix1[j][k] / matrix2[j][k];
        }
    }
    return res;
}

double** multi_matrix(double** matrix1, double** matrix2, int rowNum_matrix1, int colNum_matrix1, int colNum_matrix2){
    int i,j,k;
    double **res;


    res = malloc(rowNum_matrix1 * sizeof(double *));
    for (i = 0;i<rowNum_matrix1;i++) {
        res[i] = calloc(colNum_matrix2, sizeof(double));
    }
    
    for (i=0;i<rowNum_matrix1;i++){
        for(j=0;j<colNum_matrix2;j++){
            for(k = 0;k<colNum_matrix1;k++){
                res[i][j]+= matrix1[i][k]*matrix2[k][j];
                
                  }
        }
    }

    return res;
}

double** mul_onebyone(double** matrix1, double** matrix2, int rowNum, int colNum){
    double** res;
    int i,j,k;
    res = malloc(rowNum * sizeof(double *));
    for (i = 0;i<rowNum;i++) {
        res[i] = calloc(colNum, sizeof(double));
    }
    for (j = 0; j<rowNum; j++){
        for (k = 0; k<colNum; k++){
            res[j][k] = matrix1[j][k] * matrix2[j][k];
        }
    }
    return res;
}

double** symnmf(double** H, double** W, int rowsW,int rowsH, int colsH){
    double** cur_H = NULL;
    int  j=0;
    int check = 0;
    double conv;
    while (check == 0 && j < iter){
        cur_H = updateH(W, H, rowsW, rowsH, colsH);
        conv = converge(cur_H, H, rowsH, colsH);
        H = cur_H;
        if (conv < eps){
            check = 1;
        }
        j++;

    }
    return cur_H;

}

void print1(double **mat, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%.4f", mat[i][j]);
            }
            if (j < n - 1) {
                printf(",");
            }
            else {
                printf("\n");
            }
        }
    }


void print2(double **mat, int n,int a) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < a; j++) {
            printf("%.4f", mat[i][j]);
            }
            if (j < a - 1) {
                printf(",");
            }
            else {
                printf("\n");
            }
        }
    }



int main(int argc, char *argv[]){
    int i,j;
    char *goal;
    double **data;
    double **res;
    if (argc==3){
        goal = argv[1];
        data = load_data(argv[2]);
        res = malloc(num_of_vectors* sizeof(double *));
        for (i = 0;i<num_of_vectors;i++) {
            res[i] = calloc(num_of_vectors, sizeof(double));
        }
        if (strcmp(goal, "sym") == 0){
            res = sym(data,num_of_vectors,vector_size);

        }
        else if (strcmp(goal, "ddg") == 0){
            res = ddg(data,num_of_vectors,vector_size);
        }
        else if (strcmp(goal, "norm") == 0){
            res = norm(data,num_of_vectors,vector_size);

        }
        print1(res,num_of_vectors);
        for (j= 0; j< num_of_vectors; j++){
            free(data[j]);
            free(res[j]);
        }
        free(data);
        free(res);
        exit(1);
    }
    else{
        printf("An error has occurred");
        exit(1);
    }
}

