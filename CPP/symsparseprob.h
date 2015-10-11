
/******************************************************************************

                                    CG SOLVER
                                    EMANT 2013

******************************************************************************/

#ifndef SYMSPARSEPROB_H_INCLUDED
#define SYMSPARSEPROB_H_INCLUDED

#include "sparse.h"

using namespace std;

class SymSparseProb
{
    public:
        unsigned long  N;
        double *X,*B;
        Sparse A;
        SymSparseProb(unsigned long nrow,unsigned int mxi = 50000,
                      double tlr = 1e-15, unsigned int upd = 2000);
        ~SymSparseProb();

        void setToler       (double tlr);
        void setMaxiter     (unsigned int mxi);
        void setMatrixValue (unsigned long i, unsigned long j,
                            double valu , bool replace);
        void setIterUpdate  (double upd);
        void printOut();
        double dotp(double *V1, double *V2);

        void CGSolve();
        void JCGSolve();
        void SSORCGSolve(double omega = 1.);
        void SSORprec(double *X, double*Y, double omega);

    private:
        unsigned int MAX_ITER;
        double CG_TOLER;
        unsigned int ITER_UPDATE;
        unsigned int iter;
        double toler;
        bool out_of_iter;
};

inline double SymSparseProb::dotp(double *V1, double *V2)
{
    double d = 0;
    for (unsigned long i = 0 ;i<N; i++)
        d += V1[i]*V2[i];
    return d;

};

inline void SymSparseProb::SSORprec(double *X, double *Y, double omega)
{
    //*************************************************************************
    //
    //      PRECONDIZIONATORE SSOR (Evans, Axelsson - 1974)
    //      v1.0 - Emant 2013
    //
    //  [G. Meurant - Computer Solution of Large Linear Systems, 1999 - Cap.8]
    //*************************************************************************
    //
    // Per una mattrice simmetrica il precondizionatore SSOR è :
    // P = 1/(omega*(2-omega))*(D+omega*L)*inv(D)*(D+omega*L'), omega = ]0..2[
    //
    // Il calcolo di Y = inv(P)*X avviene in 3 fasi:
    //
    // 1) [D/omega + L]*z = X*(2-omega)/omega
    //
    // 2) [D^-1]*w = z;
    //
    // 3) D/omega + U]*Y =w;
    //
    // N.B - Nella procedura descritta, dividere per omega nella 1) e 3)
    // equivale a moltiplicare Y per omega^2 ! Per tale ragione nella 1)
    // moltiplicando per (2-omega)/omega si ottiene il risultato corretto.
    //*************************************************************************

    unsigned long i;
    double c= (2-omega)/omega;
    Entry *el;

    // Risolve [D/omega + L]*z = X
    for(i=0;i<N;i++)
        Y[i] = c*X[i];                          // z[i] = (X(i] - z[0]*L[0,i] - z[i]*L[1,i]-...)/D[i].
                                                // Si inizia assegnando z = X ( Ovviamente z = Y).

    for(i=0;i<N;i++)                            // FORWARD ELIMINATION
    {
        Y[i] /= A.diag[i].val;                  // z[0..i-1] è noto. Si completa il calcolo di z[i] dvidendo per D[i]
        el = A.diag[i].first;
        while(el)                               // Si calcola il contributo di z[i] nella COLONNE sottostanti
        {
            Y[el->col] -= el->val*Y[i]*omega;
            el = el->next;
        }
    }
    //Risolve [D^-1]*w = z;
    for(i=0;i<N;i++)
        Y[i] *= A.diag[i].val;

    //Risolve [D/omega + U]*Y=w;
    for(i=N;i>0;i--)                            // BACK SUBSTITUTION
    {
        el=A.diag[i-1].first;
        while(el)
        {
            Y[i-1] -= el->val*Y[el->col]*omega;         // Nella riga i-esima du L' (cioè colonna i-esima di L)
            el = el->next;                              // si calcola il contributo di w[0..i] ad Y[i]
        }
        Y[i-1] /= A.diag[i-1].val;                      // Si conclude dividendo per D[i]
    }
}

#endif // SYMSPARSEPROB_H_INCLUDED
