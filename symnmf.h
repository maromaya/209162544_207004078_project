#ifndef UNTITLED_SYMNMF_H
#define UNTITLED_SYMNMF_H

#endif 

double** norm (double ** data_points, int num_of_rows, int row_size );
double** ddg (double ** data_points, int num_of_rows, int row_size );
double** sym (double ** data_points, int num_of_rows,int row_size );
double** symnmf(double ** H, double** W, int rowsW,int rowsH, int colsH);
double sum_row(double *lst,int n);
double sum_vectors(double *lst1,double *lst2,int n);
double** updateH(double** W, double**H, int rowsW, int rowsH, int colsH);
double converge(double** cur_H, double** old_H, int num_of_rows, int num_of_cols);
double** transpose(double ** matrix, int num_of_rows, int num_of_cols);
double** divide_mat(double** matrix1, double** matrix2, int num_of_rows, int num_of_cols);
double** mul_diag_norm_left(double** matrix1, double** matrix2, int n);
double** mul_diag_norm_right(double** matrix1, double** matrix2, int n);
double** multi_matrix(double** matrix1, double** matrix2, int num_of_rows_matrix1, int num_of_rows_matrix2, int num_of_cols_matrix2);
double** calcB(double** matrix, int num_of_rows, int num_of_cols);
double **load_data(const char* path) ;
double ED(double *lst1,double *lst2,int n);
double** mul_onebyone(double** matrix1, double** matrix2, int num_of_rows, int num_of_cols);
void print1(double **matrix, int n);
void print2(double **matrix, int n,int a);

