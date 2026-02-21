
/******************************************************************************

                                    CG SOLVER
                                    EMANT 2013

******************************************************************************/

#include <cmath>
#include "symsparseprob.h"
#include <iostream>


SymSparseProb::SymSparseProb( unsigned long nrow,
                unsigned int mxi, double tlr,unsigned int upd) : A(nrow)
{
    N = nrow;
    X = new double[N]();
    B = new double[N]();

    CG_TOLER = tlr;
    MAX_ITER = mxi;
    ITER_UPDATE = upd;
    iter = 0;
    toler = 0.;
    out_of_iter = false;

}

SymSparseProb::~SymSparseProb()
{
    delete X;
    delete B;
}

void SymSparseProb::setToler( double tlr)
{
    CG_TOLER = tlr;
}

void SymSparseProb::setMaxiter(unsigned int mxi)
{
    MAX_ITER = mxi;
}

void SymSparseProb::setIterUpdate(double upd)
{
    ITER_UPDATE = upd;
}

void SymSparseProb::setMatrixValue(unsigned long i , unsigned long j,
                                   double valu , bool replace)
{
    A.setValue(i,j,valu,replace);
}

void SymSparseProb::printOut()
{
    cout << " - SymSparseProb Data : -"<<endl<<endl;
    cout << "[A]";
    A.printOut();
    cout << "[B]        [X]"<<endl<<endl;
    for(unsigned long i = 0; i<N ; i++)
        cout <<"["<<B[i]<<"]        ["<<X[i]<<"]"<<endl;
    cout <<endl<<endl;
}



void SymSparseProb::CGSolve()
{
    // GRADIENTI CONIUGATI SENZA PRECONDIZIONATORE

    unsigned long i,j,k;
    double *Apk,*P,*R,
            alpha,beta,
            d1,d2;

    R = new double[N]();
    P = new double[N]();
    Apk = new double[N]();
    double normB = 0.;
    out_of_iter = false;

    for(i=0;i<N;i++)
        normB += B[i]*B[i];                         //Calcola |B| per il criterio di convergenza.
    normB = sqrt(normB);

    A.gaxpy(X,R);                                   // R = B - M*X_start // Calcolo Residui da X iniziale
    for(i=0;i<N;i++)
    {
        R[i] = B[i] - R[i];
        P[i] = R[i];                                //La prima direzione di ricerca è P0 = R0

    }

    j=0;
    k=0;
    do
    {
        d1 = 0.;                                    // Calcola alpha : min(F(xk + aplha*pk))
        d2 = 0.;
        A.gaxpy(P,Apk);
        for(i=0;i<N;i++)
        {
            d1 += R[i]*R[i];
            d2 += P[i]*Apk[i];
        }
        alpha = d1/d2;

        d2 = 0;
        if(j==ITER_UPDATE)                          // Ogni ITER_UPDATE iterazioni calcola i residui come R = B - AX
        {                                           // per limitare l'accumulo di errore d1erico
            for(i=0;i<N;i++)
                X[i] += alpha*P[i];

            A.gaxpy(X,R);
            for(i=0;i<N;i++)
            {
                R[i] = B[i]-R[i];
                d2 += R[i]*R[i];
            }
            j=0;
        }
        else
        {
            for(i=0;i<N;i++)                        // Nelle iterazioni ordinarie i residui sono
            {                                       // calcolati come R = Xk-1 * alpha*pk-1
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                d2 += R[i]*R[i];
            }
        }
        beta = d2/d1;                               // beta = rk'*rk/rk-1'*rk-1 E il valore che minimizza F(x) nella direzione pk
        for(i=0;i<N;i++)                            // Calcolo della nuova direzione di riceerca pk = rk-1 + beta*pk-1
        {
            P[i] = R[i] + beta*P[i];
        }
        j++; k++;
        if (k > MAX_ITER) out_of_iter = true;
    }while(sqrt(d2)>CG_TOLER*normB && !out_of_iter); // Controllo Convergenza. d2 = (Rk*Rk) cioè |R|^2
    delete Apk ;
    delete P;
    delete R;
    iter = k;
    toler = sqrt(d2);
}

void SymSparseProb::JCGSolve()
{
    //***********************************
    // Il precondizionatore di Jacobi è:
    // J = diag(A)
    // Quindi si calcola z = inv(J)*r
    //***********************************

    unsigned long i,j,k;
    double *Apk,*P,*R,*J,*Z,
            alpha,beta,d1,d2,norm=0.;

    R = new double[N]();
    P = new double[N]();
    Apk = new double[N]();
    J = new double[N]();
    Z = new double[N]();
    double normB = 0.;
    out_of_iter = false;

    for(i=0;i<N;i++)
        normB += B[i]*B[i];                 //Calcola |B| per il criterio di convergenza.
    normB = sqrt(normB);

    A.gaxpy(X,R);                           // R = B - M*X_start // Calcolo Residui da X iniziale
    for(i=0;i<N;i++)
    {
        J[i] = 1/A.diag[i].val;
        R[i] = B[i] - R[i];
        Z[i] = J[i]*R[i];
        P[i] = Z[i];
    }
    j=0;
    k=0;
    do
    {
        d1 = 0.;                           // Calcola alpha : min(F(xk + aplha*pk))
        d2 = 0.;
        A.gaxpy(P,Apk);
        for(i=0;i<N;i++)
        {
            d1 += R[i]*Z[i];
            d2 += P[i]*Apk[i];
        }
        alpha = d1/d2;

        d2 = 0;
        norm = 0;
        if(j==ITER_UPDATE)                  // Ogni ITER_UPDATE iterazioni calcola i residui come R = B - AX
        {                                   // per limitare l'accumulo di errore numerico
            for(i=0;i<N;i++)
                X[i] += alpha*P[i];

            A.gaxpy(X,R);
            {
                R[i] = B[i]-R[i];
                Z[i] = J[i]*R[i];
                d2 += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            j=0;
        }
        else
        {
            for(i=0;i<N;i++)
            {
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                Z[i] = J[i]*R[i];
                d2 += Z[i]*R[i];
                norm += R[i]*R[i];
            }
        }
        beta = d2/d1;                       // beta = zk'*rk/zk-1'*rk-1. E il valore che minimizza F(x) nella direzione pk
        for(i=0;i<N;i++)                    // Calcolo della nuova direzione di riceerca pk = rk-1 + beta*pk-1
        {
            P[i] = Z[i] + beta*P[i];
        }
        j++; k++;
        if (k > MAX_ITER) out_of_iter = true;
    }while(sqrt(norm)>CG_TOLER*normB && !out_of_iter); // Controllo Convergenza. norm = (Rk*Rk) cioè |R|^2
    delete Apk ;
    delete P;
    delete R;
    iter = k;
    toler = sqrt(norm);
}


void SymSparseProb::SSORCGSolve(double omega)
{
    //*************************************************************************
    //
    //      PRECONDIZIONATORE SSOR (Evans, Axelsson - 1974)
    //      v1.0 - Emant 2013
    //
    //  [G. Meurant - Computer Solution of Large Linear Systems, 1999 - Cap.8]
    //*************************************************************************
    //
    // Per una matrice simmetrica il precondizionatore SSOR è :
    // P = 1/(omega*(2-omega))*(D+omega*L)*inv(D)*(D+omega*L'), omega = ]0..2[
    //
    // omega default = 1
    //*************************************************************************

    unsigned long i,j,k;
    double *Apk,*P,*R,*Z,
            alpha,beta,norm,
            d1,d2;

    R = new double[N]();
    P = new double[N]();
    Apk = new double[N]();
    Z = new double[N]();
    double normB = 0.;
    out_of_iter = false;

    for(i=0;i<N;i++)
        normB += B[i]*B[i];                         //Calcola |B| per il criterio di convergenza.
    normB = sqrt(normB);

    A.gaxpy(X,R);                                   // R = B - M*X_start // Calcolo Residui da X iniziale
    for(i=0;i<N;i++)
    {
        R[i] = B[i] - R[i];
        SSORprec(R,Z,omega);
        P[i] = Z[i];                                //La prima direzione di ricerca è P0 = R0
    }

    j=0;
    k=0;

    do
    {
        d1 = 0.;                                    // Calcola alpha : min(F(xk + aplha*pk))
        d2 = 0.;
        A.gaxpy(P,Apk);
        for(i=0;i<N;i++)
        {
            d1 += R[i]*Z[i];
            d2 += P[i]*Apk[i];
        }
        alpha = d1/d2;
        norm = 0.;
        d2 = 0;
        if(j==ITER_UPDATE)                          // Ogni ITER_UPDATE iterazioni calcola i residui come R = B - AX
        {                                           // per limitare l'accumulo di errore d1erico
            for(i=0;i<N;i++)
                X[i] += alpha*P[i];

            A.gaxpy(X,R);
            for(i=0;i<N;i++)
            {
                R[i] = B[i]-R[i];
                SSORprec(R,Z,omega);
                d2 += Z[i]*R[i];
                norm += R[i]*R[i];
            }
            j=0;
        }
        else
        {
            for(i=0;i<N;i++)                        // Nelle iterazioni ordinarie i residui sono
            {                                       // calcolati come R = Xk-1 * alpha*pk-1
                X[i] += alpha*P[i];
                R[i] -= alpha*Apk[i];
                SSORprec(R,Z,omega);
                d2 += Z[i]*R[i];
                norm += R[i]*R[i];
            }
        }
        beta = d2/d1;                               // beta = rk'*rk/rk-1'*rk-1 E il valore che minimizza F(x) nella direzione pk
        for(i=0;i<N;i++)                            // Calcolo della nuova direzione di riceerca pk = rk-1 + beta*pk-1
        {
            P[i] = R[i] + beta*P[i];
        }
        j++; k++;
        if (k > MAX_ITER) out_of_iter = true;
        cout << sqrt(norm) <<","<< k <<endl;
    }while(sqrt(norm)>CG_TOLER*normB && !out_of_iter); // Controllo Convergenza. d2 = (Rk*Rk) cioè |R|^2
    delete Apk ;
    delete P;
    delete R;
    iter = k;
    toler = sqrt(norm);
}




