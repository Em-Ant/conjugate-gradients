#ifndef CG_H_INCLUDED
#define CG_H_INCLUDED

#define CG_TOLER 1e-12
#define ITER_R_UPDATE 2000

#define E_ADD 0
#define E_REP 1

typedef struct _rwe
{
    double val;         /*Element on Row,*/
    unsigned long ind;   /* in Column ind.*/
    struct _rwe *next;  /* Pointer to Next Element in Row*/
}r_el;

typedef struct
{
    double val;         /* Element on DIAGONAL*/
    r_el *r;            /* Pointer to Elements on that Row*/
}ssparse_mtrx;

int ssparse_multAV(ssparse_mtrx *A,
                   double *V,
                   double *R,
                   unsigned long n);

int ssparse_free(ssparse_mtrx *A,
                 unsigned long n);

int ssparse_GET(double *V,
                unsigned long i,
                unsigned long j,
                ssparse_mtrx *M,
                char remove);

int ssparse_set(double value,
                unsigned long i,
                unsigned long j,
                ssparse_mtrx *A,
                char flag);

int ssparse_init(ssparse_mtrx **A,
                 double **B,
                 double **X,
                 unsigned long n);

int ssparse_CG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n);

int ssparse_JCG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n);

int ssparse_PCG(ssparse_mtrx *M,
              double *B,
              double *X,
              unsigned long n);

int ssparse_SSOR(ssparse_mtrx *M,
              double *X,
              double *Y,
              unsigned long n,
              double omega);

int ssparse_DIRICHLET(ssparse_mtrx *M,
              double *B,
              double V,
              unsigned long i);
#endif 
/* CG_H_INCLUDED */
