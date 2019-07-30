#include <vector>
#include <iostream>

#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"

using std::vector;
using boost::numeric::ublas::matrix;
using std::cout, std::endl;


matrix<int> make_matrix(vector<vector<int>> values) {
    auto the_matrix = matrix<int>(values.size(), values[0].size());

    // Copy 'values' into 'the_matrix', one row at a time
    auto iter1 = the_matrix.begin1();
    for (auto values_iter = values.begin(); values_iter != values.end(); ++values_iter, ++iter1)
        std::copy(values_iter->begin(), values_iter->end(), iter1.begin());

    return the_matrix;
}

matrix<int> make_matrix(int size1, int size2, int value) {
    auto the_matrix = matrix<int>(size1, size2);

    for (auto iter1 = the_matrix.begin1(); iter1 != the_matrix.end1(); ++iter1)
        for (auto iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
            *iter2 = value;

    return the_matrix;
}

int main() {
    matrix<int> the_matrix = make_matrix({{1, 0, 2, 4, 3},
                                          {4, 6, 5, 2, 1},
                                          {4, 4, 5, 2, 1},
                                          {5, 6, 8, 5, 3}});

    cout << "Print the matrix by rows:" << endl;
    for (auto iter1 = the_matrix.begin1(); iter1 != the_matrix.end1(); ++iter1)
        for (auto iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
            cout << *iter2 << " ";

    cout << endl << endl << "Print the matrix by columns:" << endl;
    for (auto iter2 = the_matrix.begin2(); iter2 != the_matrix.end2(); ++iter2)
        for (auto iter1 = iter2.begin(); iter1 != iter2.end(); ++iter1)
            cout << *iter1 << " ";

    cout << endl;
}