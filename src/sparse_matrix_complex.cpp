#include "sparse_matrix_complex.h"
#include <iostream>

ComplexSparseMatrix::ComplexSparseMatrix() : number_of_rows(0), diagonal(0)
{
}

ComplexSparseMatrix::ComplexSparseMatrix(unsigned long n)
{
    number_of_rows = n;
    diagonal = new ComplexDiagEntry[number_of_rows]();
}

ComplexSparseMatrix::~ComplexSparseMatrix()
{
    for(unsigned long index = 0; index < number_of_rows; index++)
    {
        while(diagonal[index].first)
        {
            ComplexEntry *delete_element = diagonal[index].first;
            diagonal[index].first = delete_element->next;
            delete delete_element;
        }
    }
    delete[] diagonal;
}

void ComplexSparseMatrix::setValue(unsigned long row_index, unsigned long column_index, complex<double> value, bool replace)
{
    if (column_index >= number_of_rows) return;

    if(row_index == column_index)
    {
        if(replace)
            diagonal[row_index].val = value;
        else
            diagonal[row_index].val += value;
    }
    else if(diagonal[row_index].first == 0)
        diagonal[row_index].first = new ComplexEntry(value, column_index);
    else
    {
        ComplexEntry *current = diagonal[row_index].first;
        ComplexEntry *previous = 0;
        while(current && current->col < column_index)
        {
            previous = current;
            current = current->next;
        }

        if (!current)
            previous->next = new ComplexEntry(value, column_index);
        else if(current->col == column_index)
        {
            if(replace)
                current->val = value;
            else
                current->val += value;
        }
        else
        {
            if(previous)
                previous->next = new ComplexEntry(value, column_index, current);
            else
                diagonal[row_index].first = new ComplexEntry(value, column_index, current);
        }
    }
}

void ComplexSparseMatrix::gaxpy(complex<double> *input_vector, complex<double> *result_vector)
{
    for(unsigned long index = 0; index < number_of_rows; index++)
        result_vector[index] = diagonal[index].val * input_vector[index];

    for(unsigned long index = 0; index < number_of_rows; index++)
    {
        ComplexEntry *current = diagonal[index].first;
        while(current)
        {
            result_vector[index] += current->val * input_vector[current->col];
            result_vector[current->col] += current->val * input_vector[index];
            current = current->next;
        }
    }
}

void ComplexSparseMatrix::printOut()
{
    std::cout << std::endl << std::endl;
    for(unsigned long index = 0; index < number_of_rows; index++)
    {
        std::cout << index << " - <[" << diagonal[index].val << "]>";
        ComplexEntry *pointer = diagonal[index].first;
        while(pointer)
        {
            std::cout << "->[" << pointer->val << "]," << pointer->col;
            pointer = pointer->next;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
}
