#ifndef SPARSE_H_INCLUDED
#define SPARSE_H_INCLUDED

using namespace std;

class Entry                                             // Matrice in formato ListOfLists
{                                                       // Ogni riga è una lista dei valori non nulli SOPRA la diagonale;


    public:
        double val;                                     // Valore  non-zero
        unsigned long col;                              // indice COLONNA
        Entry *next;                                    // Puntatore al prossimo elemento della lista RIGA
        Entry( double valu, unsigned long  colu,
                Entry *nxt=0):val(valu),col(colu),next(nxt){};
};

class DiagEntry
{
    public:
        double val;
        Entry *first;
        DiagEntry():first(0){};
};

class Sparse                                                    // Matrice Simmetrica
{
    public:
        Sparse():N(0),diag(0){};
        Sparse(unsigned long nrow);
        ~Sparse();
        void setValue(unsigned long i,
                 unsigned long j,
                 double valu,
                 bool replace);
        void gaxpy( double *V,                                  // Prodotto Matrice*Vettore (nella terminologia LAPACK)
                    double *R);
        void printOut();

    //private:
        unsigned long N;
        DiagEntry *diag;
};

inline void Sparse::gaxpy( double *V , double*R)               // Moltiplicazione MATRICE*VETTORE (nella terminologia LAPACK) R = M*V
{

    for(unsigned long i=0 ;i < N ; i++)
        R[i] = diag[i].val*V[i];

    for(unsigned long i=0 ; i < N ; i++)
    {
        Entry *curr  = diag[i].first;
        while(curr)
        {
            R[i] += curr->val*V[curr->col];
            R[curr->col] += curr->val*V[i];
            curr = curr->next;
        }
    }
}

#endif // SPARSE_H_INCLUDED
