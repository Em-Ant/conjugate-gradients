#ifndef SPARSE_MATRIX_H_INCLUDED
#define SPARSE_MATRIX_H_INCLUDED

#include "linear_solver_traits.h"
#include <iostream>
#include <algorithm>

template <typename Traits>
class SparseMatrix
{
public:
    using Scalar = typename Traits::Scalar;

    struct Entry
    {
        Scalar val;
        unsigned long col;
        Entry *next;
        Entry(const Scalar &v, unsigned long c, Entry *n = nullptr)
            : val(v), col(c), next(n) {}
    };

    struct DiagEntry
    {
        Scalar val;
        Entry *first;
        DiagEntry() : first(nullptr) {}
    };

    unsigned long N;
    DiagEntry *diag;

    SparseMatrix() : N(0), diag(nullptr) {}

    SparseMatrix(unsigned long n) : N(n)
    {
        diag = new DiagEntry[N]();
    }

    ~SparseMatrix()
    {
        for (unsigned long i = 0; i < N; i++)
        {
            while (diag[i].first)
            {
                Entry *del = diag[i].first;
                diag[i].first = del->next;
                delete del;
            }
        }
        delete[] diag;
    }

    void setValue(unsigned long i, unsigned long j, const Scalar &value)
    {
        if (Traits::enforces_symmetric_storage && i > j)
        {
            std::swap(i, j);
        }
        if (j >= N)
            return;

        if (i == j)
        {
            diag[i].val = value;
        }
        else if (diag[i].first == nullptr)
        {
            diag[i].first = new Entry(value, j);
        }
        else
        {
            Entry *curr = diag[i].first, *prev = nullptr;
            while (curr && curr->col < j)
            {
                prev = curr;
                curr = curr->next;
            }

            if (!curr)
            {
                prev->next = new Entry(value, j);
            }
            else if (curr->col == j)
            {
                curr->val = value;
            }
            else
            {
                if (prev)
                {
                    prev->next = new Entry(value, j, curr);
                }
                else
                {
                    diag[i].first = new Entry(value, j, curr);
                }
            }
        }
    }

    void addToValue(unsigned long i, unsigned long j, const Scalar &value)
    {
        if (Traits::enforces_symmetric_storage && i > j)
        {
            std::swap(i, j);
        }
        if (j >= N)
            return;

        if (i == j)
        {
            diag[i].val += value;
        }
        else if (diag[i].first == nullptr)
        {
            diag[i].first = new Entry(value, j);
        }
        else
        {
            Entry *curr = diag[i].first, *prev = nullptr;
            while (curr && curr->col < j)
            {
                prev = curr;
                curr = curr->next;
            }

            if (!curr)
            {
                prev->next = new Entry(value, j);
            }
            else if (curr->col == j)
            {
                curr->val += value;
            }
            else
            {
                if (prev)
                {
                    prev->next = new Entry(value, j, curr);
                }
                else
                {
                    diag[i].first = new Entry(value, j, curr);
                }
            }
        }
    }

    void gaxpy(Scalar *V, Scalar *R)
    {
        for (unsigned long i = 0; i < N; i++)
        {
            R[i] = diag[i].val * V[i];
        }

        for (unsigned long i = 0; i < N; i++)
        {
            Entry *curr = diag[i].first;
            while (curr)
            {
                R[i] += curr->val * V[curr->col];
                R[curr->col] += curr->val * V[i];
                curr = curr->next;
            }
        }
    }

    void printOut()
    {
        std::cout << std::endl
                  << std::endl;
        for (unsigned long i = 0; i < N; i++)
        {
            std::cout << i << " - <[" << diag[i].val << "]>";
            Entry *p = diag[i].first;
            while (p)
            {
                std::cout << "->[" << p->val << "]," << p->col;
                p = p->next;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl
                  << std::endl;
    }
};

using Sparse = SparseMatrix<RealTraits>;
using ComplexSparseMatrix = SparseMatrix<ComplexSymmetricTraits>;

#endif
