
#include "sparse_matrix.h"
#include <iostream>

using namespace std;

Sparse::Sparse(unsigned long nrow)
{
    N = nrow;
    diag = new DiagEntry[N]();
}

Sparse::~Sparse()
{

    for(unsigned long i = 0; i < N ; i++)
    {
        while(diag[i].first)
        {
            Entry *del = diag[i].first;
            diag[i].first = del->next;
            delete del;
        }
    }
    delete diag;
}


void Sparse::setValue(unsigned long i,unsigned long j,
                    double valu,bool replace)
{

    if(i>j){unsigned long k = i ; i = j ; j = k;}
    if (j>=N) return;

    if(i==j)                                            // Update diagonal
    {
        if(replace)
            diag[i].val = valu;                         // Replace or add coefficients
        else
            diag[i].val += valu;
    }

    else if(diag[i].first == 0)                         // Row is empty
        diag[i].first = new Entry(valu,j);             // Create new entry

    else
    {
        Entry *curr = diag[i].first, *prev = 0;        // Initialize search

        while(curr && curr->col < j)                    // Scan row list forward
        { prev = curr ; curr = curr->next; }

        if (!curr)                                      // Create new entry at end
            prev->next = new Entry(valu,j);

        else if(curr->col == j)                         // Non-null value at position i,j
        {
            if(replace)
                curr->val += valu;
            else
                curr->val = valu;
        }

        else                                            // Insert in middle of row list
        {
            if(prev)
                prev->next = new Entry(valu,j,curr);
            else                                        // Insert at beginning of list
                diag[i].first = new Entry(valu,j,curr);
        }
    }
}


void Sparse::printOut()
{
    cout << endl<<endl;
    for(unsigned long i=0 ;i < N; i++)
    {
        cout<<i<<" - <["<<diag[i].val<<"]>";
        Entry *p = diag[i].first;
        while(p)
        {
            cout<<"->["<<p->val<<"],"<<p->col;
            p = p->next;
        }
        cout <<endl;
    }
    cout << endl<<endl;
}
