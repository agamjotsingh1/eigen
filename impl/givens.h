#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

compl** gmat(int m, int i, int j, compl** vec){
    compl** mat = meye(m);
    compl xi = vec[i][0];
    compl xj = vec[j][0];

    compl c = conj(xi)/sqrt(pow(cabs(xi), 2) + pow(cabs(xj), 2));
    compl s = conj(xj)/sqrt(pow(cabs(xi), 2) + pow(cabs(xj), 2));

    mat[i][i] = c;
    mat[i][j] = s;
    mat[j][i] = -conj(s);
    mat[j][j] = conj(c);
    return mat;
}

// H is hessenberg form
compl*** givens(compl** H, int m){
    compl** Q = meye(m);

    for(int i = 0; i < m - 1; i++){
        compl** vec = mgetcol(H, m, m, i);
        compl** G = gmat(m, i, i + 1, vec);
        Q = mmul(Q, mT(G, m, m), m, m, m);
        H = mmul(G, H, m, m, m); 
    }

    compl*** ret = (compl***) malloc(sizeof(compl**)*2);
    ret[0] = Q;
    ret[1] = H;
    return ret;
}