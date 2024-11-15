#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

// Gives the givens matrix associated with vector 'vec' at position (i, j)
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

// Gives QR decomposition of the hessenberg matrix H
compl** givens(compl** H, int m){
    for(int i = 0; i < m - 1; i++){
        compl** vec = mgetcol(H, m, m, i);
        compl** G = gmat(m, i, i + 1, vec);
        H = mmul(G, H, m, m, m); 
        H = mmul(H, mT(G, m, m), m, m, m); 
    }

    return H;
}