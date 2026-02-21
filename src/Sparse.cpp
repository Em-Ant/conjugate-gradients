
#include "sparse.h"
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

    if(i==j)                                            // Aggiorna la diagonale
    {
        if(replace)
            diag[i].val = valu;                         // Sostituzione o Somma dei Coefficienti
        else
            diag[i].val += valu;
    }

    else if(diag[i].first == 0)                         // La riga è vuota
        diag[i].first = new Entry(valu,j);             // quindi crea un nuovo Entry

    else
    {
        Entry *curr = diag[i].first, *prev = 0;        // Inizializza la ricerca

        while(curr && curr->col < j)                    // Scorre la lista riga in avanti...
        { prev = curr ; curr = curr->next; }

        if (!curr)                                      // Crea un nuovo Entry in coda
            prev->next = new Entry(valu,j);

        else if(curr->col == j)                         // C'è un valore non nullo in posizione i,j
        {
            if(replace)
                curr->val += valu;
            else
                curr->val = valu;
        }

        else                                            // Inserimento in mezzo alla lista riga
        {
            if(prev)
                prev->next = new Entry(valu,j,curr);
            else                                        // Inserimento in testa alla lista
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
