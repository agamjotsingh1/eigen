#include <complex.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

compl*** QR(compl** A, int m, int n){
    compl** Q = mzeroes(m, n);
    compl** R = mzeroes(n, n);
    for(int i = 0; i < n; i++){
        compl** a = mgetcol(A, m, n, i);
        compl** q = mdup(a, m, 1);

        for(int j = 0; j < i; j++){
            compl x = mmul(mT(a, m, 1), mgetcol(Q, m, n, j), 1, m, 1)[0][0];
            q = madd(q, mscale(mgetcol(Q, m, n, j), m, 1, -x), m, 1);
            R[j][i] = x;
        }

        R[i][i] = vnorm(q, m) + 0*I;
        q = mscale(q, m, 1, (compl) 1/vnorm(q, m));

        for(int j = 0; j < m; j++){
            Q[j][i] = q[j][0];
        }
    }

    compl*** ret = (compl***) malloc(sizeof(compl**)*2);
    ret[0] = Q;
    ret[1] = R;
    return ret;
}