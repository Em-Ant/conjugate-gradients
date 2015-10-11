
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "CG.h"

int ssparse_init(ssparse_mtrx **A,
                 double **B,
                 double **X,
                 unsigned long n)
{
    *A = calloc(n,sizeof(ssparse_mtrx));
    *B = calloc(n,sizeof(double));
    *X = calloc(n,sizeof(double));
    return 0;
}

int ssparse_set(double value,
                unsigned long i,
                unsigned long j,
                ssparse_mtrx *A,
                char flag)
{
    unsigned long app;
    r_el *el,*prev;

    if(i>j)
    {
        app = i;
        i = j;
        j = app;
    }

    if(i==j)
    {
        switch (flag)
        {
            case E_REP:
            {
                A[i].val = value;
                break;
            }
            case E_ADD:
                A[i].val += value;
                break;
        }
    }

    else if(A[i].r == 0)
    {
        A[i].r = calloc(1,sizeof(r_el));
        A[i].r->ind = j;
        A[i].r->val = value;
    }
    else
    {
        el = A[i].r;
        prev = 0;
        while(el)
        {
            if(el->ind == j)
            {
                switch (flag)
                {
                    case E_ADD:{
                        el->val += value;
                        break;
                    }
                    case E_REP:{
                        el->val = value;
                        break;
                    }
                }
                return 0;
            }
            else if(el->ind > j)
                break;
            prev = el;
            el = el->next;
        }
        if (el)
        {
            if(prev)
            {
                prev->next = malloc(sizeof(r_el));
                prev->next->val = value;
                prev->next->ind = j;
                prev->next->next = el;
            }
            else
            {
                A[i].r = malloc(sizeof(r_el));
                A[i].r->val = value;
                A[i].r->ind = j;
                A[i].r->next = el;
            }
        }
        else
        {
            prev->next = malloc(sizeof(r_el));
            prev->next->val = value;
            prev->next->ind = j;
            prev->next->next = 0;
        }
    }
    return 0;
}

int ssparse_free(ssparse_mtrx *A,
                 unsigned long n)
{
    unsigned long i;
    r_el *el,*nxt;

    for(i=0;i<n;i++)
    {
        el = A[i].r;
        while(el)
        {
            nxt = el->next;
            free(el);
            el = nxt;
        }
    }
    free(A);
    return 0;
}

int ssparse_multAV(ssparse_mtrx *A,
                   double *V,
                   double *R,
                   unsigned long n)
{
    unsigned long i;
    r_el *el;

    for(i=0;i<n;i++)
    {
        R[i] = A[i].val*V[i];
    }
    for(i=0;i<n;i++)
    {
        el  = A[i].r;
        while(el)
        {
            R[i] += el->val*V[el->ind];
            R[el->ind] += el->val*V[i];
            el = el->next;
        }
    }
    return 0;
}

/* Simple CG*/
int ssparse_CG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n)
{
    unsigned long i,j,k;
    double *Apk,*P,*R,
            alpha,beta,	/*norm,*/
            num,den;

    
    #ifdef _OUT
		printf("\nCONJUGATE GRADIENT SOLVER for SYMMETRIC SPARSE MATRICES,\
		    v1.0 - Emant (R), 2010.\n\nAllocating Memory...\n");
	#endif
	
    R = malloc(n*sizeof(double));
    P = malloc(n*sizeof(double));
    Apk = malloc(n*sizeof(double));

	#ifdef _OUT
    	printf("Setting Starting Point...\n");
	#endif
    ssparse_multAV(M,X,R,n);

    for(i=0;i<n;i++)
    {
        R[i] = B[i] - R[i];
        P[i] = R[i];
    }
    j=0;
    k=0;
    #ifdef _OUT
    	printf("Starting Iterations. Please Wait...\n\n");
    #endif
    do
    {
        num = 0;
        den = 0;
        ssparse_multAV(M,P,Apk,n);

        for(i=0;i<n;i++)
        {
            num += R[i]*R[i];
            den += P[i]*Apk[i];
        }
        alpha = num/den;
        den = 0;
        /*norm = 0;*/
        if(j==ITER_R_UPDATE)
        {
            for(i=0;i<n;i++)
                X[i] += alpha*P[i];

            ssparse_multAV(M,X,R,n);
            for(i=0;i<n;i++)
            {
                R[i] = B[i]-R[i];
                den += R[i]*R[i];
            }
    	#ifdef _OUT
            printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(den),CG_TOLER);
        #endif
            j=0;
        }
        else
        {
            for(i=0;i<n;i++)
            {
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                den += R[i]*R[i];
            }
    	#ifdef _OUT
            if(j==ITER_R_UPDATE/2)
                printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(den),CG_TOLER);
        #endif
        }
        beta = den/num;
        for(i=0;i<n;i++)
        {
            P[i] = R[i] + beta*P[i];
        }
        j++;
        k++;
    }while(sqrt(den)>CG_TOLER);

    free(Apk);
    free(P);
    free(R);
    #ifdef _OUT
    printf("\nSolution Done in %li Iterations.\nCurrent NORM is %5.3e.\n\
         TOLER is Set to %5.3e.\n\n",k,sqrt(den),CG_TOLER);
    #endif

    return 0;
}

/* JACOBI Proconditioner*/
int ssparse_JCG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n)
{
    unsigned long i,j,k;
    double *Apk,*P,*R,*J,*Z,
            alpha,beta,norm,
            num,den;

    printf("\nJACOBI CONJUGATE GRADIENT SOLVER for SYMMETRIC SPARSE MATRICES,\n"
        "v1.0 - Emant (R), 2010.\n\nAllocating Memory...\n");

    R = malloc(n*sizeof(double));
    P = malloc(n*sizeof(double));
    Apk = malloc(n*sizeof(double));
    J = malloc(n*sizeof(double));
    Z = malloc(n*sizeof(double));

    printf("Setting Starting Point...\n");

    ssparse_multAV(M,X,R,n);

    for(i=0;i<n;i++)
    {
        J[i] = 1/M[1].val; /* Jacobi Preconditioner*/
        R[i] = B[i] - R[i];
        Z[i] = J[i]*R[i];
        P[i] = Z[i];
    }
    j=0;
    k=0;
    printf("Starting Iterations. Please Wait...\n\n");
    do
    {
        num = 0;
        den = 0;
        ssparse_multAV(M,P,Apk,n);

        for(i=0;i<n;i++)
        {
            num += R[i]*Z[i];
            den += P[i]*Apk[i];
        }
        alpha = num/den;
        den = 0;
        norm = 0;
        if(j==ITER_R_UPDATE)
        {
            for(i=0;i<n;i++)
                X[i] += alpha*P[i];

            ssparse_multAV(M,X,R,n);
            for(i=0;i<n;i++)
            {
                R[i] = B[i]-R[i];
                Z[i] = J[i]*R[i];
                den += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(norm),CG_TOLER);
            j=0;
        }
        else
        {
            for(i=0;i<n;i++)
            {
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                Z[i] = J[i]*R[i];
                den += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            if(j==ITER_R_UPDATE/2)
                printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(norm),CG_TOLER);
        }
        beta = den/num;
        for(i=0;i<n;i++)
        {
            P[i] = Z[i] + beta*P[i];
        }
        j++;
        k++;
    }while(sqrt(norm)>CG_TOLER);

    free(Apk);
    free(P);
    free(R);
    free(J);
    free(Z);

    printf("\nSolution Done in %li Iterations.\nCurrent NORM is %5.3e.\n"
        "TOLER is Set to %5.3e.\n\n",k,sqrt(norm),CG_TOLER);

    return 0;
}

int ssparse_SSOR(ssparse_mtrx *M,
              double *X,
              double *Y,
              unsigned long n,
              double omega)
{
    unsigned long i;
    double c;
    r_el *el;

    c= (2-omega)/omega;

    for(i=0;i<n;i++)
        Y[i] = c*X[i];

    /* Solve [D/omega + L]*z = X*/
    for(i=0;i<n;i++)
    {
        Y[i] /= M[i].val;
        el = M[i].r;
        while(el)
        {
            Y[el->ind] -= el->val*Y[i]*omega;
            el = el->next;
        }
    }
    /*Solve [D^-1]*w = z;*/
    for(i=0;i<n;i++)
        Y[i] *= M[i].val;

    /*Solve [D/omega + U]*Y=w;*/
    for(i=n;i>0;i--)
    {
        el=M[i-1].r;
        while(el)
        {
            Y[i-1] -= el->val*Y[el->ind]*omega;
            el = el->next;
        }
        Y[i-1] /= M[i-1].val;
    }
    return 0;
}

/* SSOR Proconditioner*/
int ssparse_PCG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n)
{
    unsigned long i,j,k;
    double *Apk,*P,*R,*Z,
            alpha,beta,norm,
            num,den,omega;

    printf("\nSSOR PRECONDITIONED CONJUGATE GRADIENT SOLVER\nfor SYMMETRIC SPARSE MATRICES,\n"
        "v1.0 - Emant (R), 2010.\n\nAllocating Memory...\n");

    R = malloc(n*sizeof(double));
    P = malloc(n*sizeof(double));
    Apk = malloc(n*sizeof(double));
    Z = malloc(n*sizeof(double));

    printf("Setting Starting Point...\n");

    omega = 0.000001;

    ssparse_multAV(M,X,R,n);

    for(i=0;i<n;i++)
    {
        R[i] = B[i] - R[i];
        ssparse_SSOR(M,R,Z,n,omega);
        P[i] = Z[i];
    }
    j=0;
    k=0;
    printf("Starting Iterations. Please Wait...\n\n");
    do
    {
        num = 0;
        den = 0;
        ssparse_multAV(M,P,Apk,n);

        for(i=0;i<n;i++)
        {
            num += R[i]*Z[i];
            den += P[i]*Apk[i];
        }
        alpha = num/den;
        den = 0;
        norm = 0;
        if(j==ITER_R_UPDATE)
        {
            for(i=0;i<n;i++)
                X[i] += alpha*P[i];

            ssparse_multAV(M,X,R,n);
            for(i=0;i<n;i++)
            {
                R[i] = B[i]-R[i];
                ssparse_SSOR(M,R,Z,n,omega);
                den += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(norm),CG_TOLER);
            j=0;
        }
        else
        {
            for(i=0;i<n;i++)
            {
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                ssparse_SSOR(M,R,Z,n,omega);
                den += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            if(j==ITER_R_UPDATE/2)
                printf("ITER: %li, NORM: %5.3e, TOLER: %5.3e\n",
                    k,sqrt(norm),CG_TOLER);
        }
        beta = den/num;
        for(i=0;i<n;i++)
        {
            P[i] = Z[i] + beta*P[i];
        }
        j++;
        k++;
    }while(sqrt(norm)>CG_TOLER);

    free(Apk);
    free(P);
    free(R);
    free(Z);

    printf("\nSolution Done in %li Iterations.\nCurrent NORM is %5.3e.\n"
        "TOLER is Set to %5.3e.\n\n",k,sqrt(norm),CG_TOLER);

    return 0;
}

int ssparse_GET(double *V,
                unsigned long i,
                unsigned long j,
                ssparse_mtrx *M,
                char remove)
{
    r_el *el,*prev;
    unsigned int a;

    if(j<i)
    {
        a=j; j=i; i=a;
    }
    else if(i==j)
    {
        *V = M[i].val;
        if(remove)
            M[i].val = 0.0;
    }
    else
    {
        el = M[i].r;
        prev = 0;
        while(el)
        {
            if(el->ind==j)
            {
                *V=el->val;
                if(remove)
                {
                    if(prev==0)
                    {
                        M[i].r = el->next;
                        free(el);
                    }
                    else
                    {
                        prev->next= el->next;
                        free(el);
                    }
                }
                break;
            }
        }
    }
    return 0;
}



int ssparse_DIRICHLET(ssparse_mtrx *M,
              double *B,
              double V,
              unsigned long i)
{
    r_el *el,*prev;
    unsigned long k;
    double val=0.;
    el = M[i].r;
    M[i].r=0;
    prev = 0;
    B[i] = M[i].val*V;
    while(el)
    {
        B[el->ind] -= el->val*V;
        prev = el;
        el = el->next;
        free(prev);
    }
    for(k=0;k<i;k++)
    {
        ssparse_GET(&val,k,i,M,1);
        B[k] -= val*V;
    }
    return 0;
}
