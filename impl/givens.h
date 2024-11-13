#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

void cprint(compl a){
    printf("%lf + %lfi", creal(a), cimag(a));
}

compl** gmat(int m, int i, int j, compl** vec){
    compl** mat = meye(m);
    /*cprint(vec[i][0]);
    printf(",");
    cprint(vec[j][0]);
    printf("\n");*/
    compl xi = vec[i][0];
    compl xj = vec[j][0];

    //compl c = vec[i][0]/(sqrt(pow((double) cabs(vec[i][0]), 2) + pow((double) cabs(vec[j][0]), 2)));
    //compl s = -vec[j][0]/(sqrt(pow((double) cabs(vec[i][0]), 2) + pow((double) cabs(vec[j][0]), 2)));
    compl c = xi/sqrt(pow(cabs(xi), 2) + pow(cabs(xj), 2));
    compl s = xj/sqrt(pow(cabs(xi), 2) + pow(cabs(xj), 2));
    /*cprint(c);
    printf(", ");
    cprint(s);
    printf("\n");*/

    mat[i][i] = c;
    mat[i][j] = s;
    mat[j][i] = -s;
    mat[j][j] = c;
    return mat;
}

// H is hessenberg form
compl*** givens(compl** H, int m){
    compl** Q = mzeroes(m, m);

    for(int i = 0; i < m - 1; i++){
        compl** vec = mgetcol(H, m, m, i);
        //mprint(vec, m, 1);
        compl** G = gmat(m, i, i + 1, vec);
        //mprint(mT(G, m, m), m, m);
        //mprint(G, m, m);
        //mprint(mmul(G, H, m, m, m), m, m);
        //mprint(mmul(mT(G, m, m), H, m, m, m), m, m);
        //mprint(mmul(G, vec, m, m, 1), m, 1);
        Q = mmul(Q, G, m, m, m);
    }

    compl** R = mmul(mT(Q, m, m), H, m, m, m);

    compl*** ret = (compl***) malloc(sizeof(compl**)*2);
    ret[0] = Q;
    ret[1] = R;
    return ret;
}