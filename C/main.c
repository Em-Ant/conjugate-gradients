#include <stdio.h>
#include <stdlib.h>
#include "CG.h"

int main()
{
    ssparse_mtrx *M;
    double *B,*X;
    unsigned long n;
    int i;
    n=6;
    ssparse_init(&M,&B,&X,n);
    ssparse_set(1.0,0,0,M,1);
    ssparse_set(1.0,1,1,M,1);
    ssparse_set(1.0,2,2,M,1);
    ssparse_set(2.0,0,3,M,1);
    ssparse_set(3.0,1,4,M,1);
    ssparse_set(4.0,2,5,M,1);

    printf("\n\n");
    for(i=0;i<6;i++)
    {
        B[i]=1.0;
    }
    ssparse_CG(M,B,X,6);
    printf("\n");
    for(i=0;i<6;i++)
    {
        printf("%8.6f ",X[i]);
    }
    printf("\n");
    ssparse_multAV(M,X,B,6);
        for(i=0;i<6;i++)
    {
        printf("%8.6f ",B[i]);
    }
    printf("\n\n");

    return 0;
}
