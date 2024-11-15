#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include "../impl/matrix.h"
#include "../impl/hess.h"
#include "../impl/givens.h"
#include "../impl/schmidt.h"

double complex cnum(double real, double imag){
    return (double complex) (real + imag*I);
}

int main(){
    int m;
    scanf("%d", &m);

    compl** mat = mzeroes(m, m);

    for(int i = 0; i < m; i++){
        for(int j = 0; j < m; j++){
            double a, b;
            scanf("%lf %lf", &a, &b);
            mat[i][j] = cnum(a, b);
        }
    }

    mat = hess(mat, m);
    //mprint(mat, m, m);

    int n = 200;
    for(int i = 0; i < n; i++) mat = schmidt(mat, m, m);

    for(int i = 0; i < m; i++) printf("%.6lf %.6lf\n", creal(mat[i][i]), cimag(mat[i][i]));
}

int main2(){
    int m = 4;
    compl** mat = mzeroes(m, m);
    compl** A = mzeroes(m, m);

    mat[0][0] = cnum(5, 0), mat[0][1] = cnum(9, 0), mat[0][2] = cnum(3, 0), mat[0][3] = cnum(5, 0);
    mat[1][0] = cnum(3, 0), mat[1][1] = cnum(3, 0), mat[1][2] = cnum(0, 0), mat[1][3] = cnum(2, 0);
    mat[2][0] = cnum(6, 0), mat[2][1] = cnum(2, 0), mat[2][2] = cnum(8, 0), mat[2][3] = cnum(9, 0);
    mat[3][0] = cnum(8, 0), mat[3][1] = cnum(7, 0), mat[3][2] = cnum(1, 0), mat[3][3] = cnum(8, 0);
    /*mat[0][0] = cnum(1, -1), mat[0][1] = cnum(2, 3), mat[0][2] = cnum(-1, 4), mat[0][3] = cnum(5, -2);
    mat[1][0] = cnum(3, 2), mat[1][1] = cnum(4, -1), mat[1][2] = cnum(2, 1), mat[1][3] = cnum(-3, 3);
    mat[2][0] = cnum(-2, 1), mat[2][1] = cnum(1, -3), mat[2][2] = cnum(3, 0), mat[2][3] = cnum(4, 2);
    mat[3][0] = cnum(5, 0), mat[3][1] = cnum(-1, -2), mat[3][2] = cnum(4, -1), mat[3][3] = cnum(3, 1);*/
    /*mat[0][0] = cnum(1, 0), mat[0][1] = cnum(2, 3), mat[0][2] = cnum(-1, 4);
    mat[1][0] = cnum(3, 2), mat[1][1] = cnum(4, -1), mat[1][2] = cnum(2, 1);
    mat[2][0] = cnum(-2, 1), mat[2][1] = cnum(1, -3), mat[2][2] = cnum(3, 0);*/
    /*A[0][0] = 1 + 1 * I;   A[0][1] = 2 + 5 * I;   A[0][2] = 3 + 9 * I;    A[0][3] = 4 + 13 * I;
    A[1][0] = 5 + 2 * I;   A[1][1] = 6 + 6 * I;   A[1][2] = 7 + 10 * I;   A[1][3] = 8 + 14 * I;
    A[2][0] = 9 + 3 * I;   A[2][1] = 10 + 7 * I;  A[2][2] = 11 + 11 * I;  A[2][3] = 12 + 15 * I;
    A[3][0] = 13 + 4 * I;   A[3][1] = 14 + 8 * I;  A[3][2] = 15 - 12 * I;  A[3][3] = 16 + 16 * I;*/

    mat = hess(mat, m);
    //mprint(mat, m, m);
    //mprint(mat, m, m);

    /*A = hess(A, m);
    mprint(A, m, m);*/
    //mprint(mmul(mat, mat, m, m, m), m, m);
    
    int n = 200;
    for(int i = 0; i < n; i++){
        //compl*** temp = schmidt(mat, m, m);
        //compl** Q = temp[0];
        //compl** R = temp[1];

        //mprint(Q, m, m);
        //mprint(R, m, m);
        //mprint(mmul(mT(Q, m, m), Q, m, m, m), m, m);
        //mprint(mmul(Q, R, m, m, m), m, m);
        
        //mat = mmul(R, Q, m, m, m);
        mat = schmidt(mat, m, m);
        //mat = givens(mat, m);
    }
    /*
    compl** vec = mzeroes(2, 1);
    vec[0][0] = 1 + 1*I;
    vec[1][0] = 2 + 3*I;
    printf("%lf", vnorm(vec, 2));*/
    /*int n = 100;
    for(int i = 0; i < n; i++){
        compl*** temp = givens(mat, m);
        compl** Q = temp[0];
        compl** R = temp[1];
        

        mprint(mat, m, m);
        mprint(mmul(Q, R, m, m, m), m, m);
        mat = mmul(R, Q, m, m, m);
    }*/

    //mprint(mat, m, m);
    for(int i = 0; i < m; i++) printf("%lf %lf ", creal(mat[i][i]), cimag(mat[i][i]));
}