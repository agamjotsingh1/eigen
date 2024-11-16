#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

int is_triangular(compl** mat, int m, double tolerance){
    for(int i = 0; i < m - 1; i++){
        if(cabs(mat[i + 1][i]) > tolerance) return 0;
    }

    return 1;
}

void shift(compl** mat, int m, compl s){
    for(int i = 0; i < m; i++){
        mat[i][i] += s;
    }
}

compl* eigen_givens(compl** mat, int m, int max_iterations, double tolerance, int* no_iterations){
    mat = hess(mat, m, tolerance);
    int i;

    for(i = 0; i < max_iterations; i++) {
        if(is_triangular(mat, m, tolerance)) break;
        compl temp = mat[m - 1][m - 1];

        shift(mat, m, -temp);
        mat = schmidt(mat, m, m, tolerance);
        shift(mat, m, temp);
    }

    *no_iterations = i;

    compl* eigen_values = malloc(sizeof(compl)*m);
    for(int i = 0; i < m; i++) {
        if(i + 1 < m && cabs(mat[i + 1][i]) > tolerance){
            compl a = mat[i][i], b = mat[i][i + 1], c = mat[i + 1][i], d = mat[i + 1][i + 1];
            compl D = (a + d)*(a + d) - 4*(a*d - b*c);
            
            eigen_values[i] = ((a + d) + csqrt(D))/2;
            eigen_values[i + 1] = ((a + d) - csqrt(D))/2;
            i++;
        }
        else {
            eigen_values[i] = mat[i][i];
        }
    }

    return eigen_values;
}