#ifndef SPARSE_MATRIX_COMPLEX_H_INCLUDED
#define SPARSE_MATRIX_COMPLEX_H_INCLUDED

#include <complex>
#include <vector>

using namespace std;

class ComplexEntry
{
    public:
        complex<double> val;
        unsigned long col;
        ComplexEntry *next;
        ComplexEntry(complex<double> valu, unsigned long colu, ComplexEntry *nxt=0)
            : val(valu), col(colu), next(nxt){};
};

class ComplexDiagEntry
{
    public:
        complex<double> val;
        ComplexEntry *first;
        ComplexDiagEntry():first(0){};
};

class ComplexSparseMatrix
{
    public:
        ComplexSparseMatrix();
        ComplexSparseMatrix(unsigned long number_of_rows);
        ~ComplexSparseMatrix();
        void setValue(unsigned long row_index, unsigned long column_index, complex<double> value);
        void addToValue(unsigned long row_index, unsigned long column_index, complex<double> value);
        void gaxpy(complex<double> *input_vector, complex<double> *result_vector);
        void printOut();

        unsigned long number_of_rows;
        ComplexDiagEntry *diagonal;
};

#endif
