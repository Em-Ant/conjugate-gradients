#ifndef SPARSE_MATRIX_H_INCLUDED
#define SPARSE_MATRIX_H_INCLUDED

using namespace std;

class Entry                                             // Matrix in ListOfLists format
{                                                       // Each row is a list of non-null values ABOVE the diagonal


    public:
        double val;                                     // Non-zero value
        unsigned long col;                              // Column index
        Entry *next;                                    // Pointer to next element in ROW list
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

class Sparse                                                    // Symmetric Matrix
{
    public:
        Sparse():N(0),diag(0){};
        Sparse(unsigned long nrow);
        ~Sparse();
        void setValue(unsigned long i, unsigned long j, double valu);
        void addToValue(unsigned long i, unsigned long j, double valu);
        void gaxpy( double *V,                                  // Matrix*Vector product (LAPACK terminology)
                    double *R);
        void printOut();

    //private:
        unsigned long N;
        DiagEntry *diag;
};

inline void Sparse::gaxpy( double *V , double*R)               // Matrix*Vector multiplication (LAPACK terminology) R = M*V
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
