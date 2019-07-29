#include <vector>
#include <limits>
#include <iostream>
#include <string>

#define CATCH_CONFIG_MAIN

#include "../include/catch.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"

using std::vector;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::zero_matrix;
using boost::numeric::ublas::scalar_matrix;
using std::cout, std::endl;
using std::max;
using std::string;

matrix<int> make_matrix(vector<vector<int>> values) {
    auto size1 = values.size();
    auto size2 = values[0].size();
    auto the_matrix = matrix<int>(size1, size2);

    for (unsigned long row = 0; row < size1; ++row)
        for (unsigned long col = 0; col < size2; ++col)
            the_matrix(row, col) = values[row][col];

    return the_matrix;
}

matrix<int> make_matrix(int size1, int size2, int value) {
    auto the_matrix = matrix<int>(size1, size2);

    // Fancier, do it with an interator
    for (unsigned long row = 0; row < size1; ++row)
        for (unsigned long col = 0; col < size2; ++col)
            the_matrix(row, col) = value;

    return the_matrix;
}

int dp_change(int money, const vector<int> &coins) {
    vector<int> min_num_coins = vector<int>(money + 1);

    for (int m = 1; m <= money + 1; ++m) {
        min_num_coins[m] = std::numeric_limits<int>::max();
        for (auto denomination: coins) {
            if (m < denomination)
                continue;
            auto n_coins = min_num_coins[m - denomination] + 1;
            min_num_coins[m] = std::min(n_coins, min_num_coins[m]);
        }
    }
    return min_num_coins[money];
}

int manhattan_tourist(const matrix<int> &down, const matrix<int> &right) {
    auto n = down.size1();
    auto m = right.size2();
    matrix<int> dp(zero_matrix<int>(n + 1, m + 1));
    // 'i' is row index and 'j' is column index
    for (unsigned long i = 1; i <= n; ++i)
        dp(i, 0) = dp(i - 1, 0) + down(i - 1, 0);
    for (unsigned long j = 1; j <= m; ++j)
        dp(0, j) = dp(0, j - 1) + right(0, j - 1);
    for (unsigned long i = 1; i <= n; ++i)
        for (unsigned long j = 1; j <= m; ++j)
            dp(i, j) = max(dp(i - 1, j) + down(i - 1, j), dp(i, j - 1) + right(i, j - 1));

    return dp(n, m);
}

TEST_CASE("align") {
    REQUIRE(dp_change(40, {50, 25, 20, 10, 5, 1}) == 2);
    REQUIRE(dp_change(19415, {18, 16, 7, 5, 3, 1}) == 1080);
    REQUIRE(dp_change(16042, {16, 13, 11, 8, 7, 5, 3, 1}) == 1003);

}

TEST_CASE("manhattan_tourist") {
    matrix<int> down = make_matrix({{1, 0, 2, 4, 3},
                                    {4, 6, 5, 2, 1},
                                    {4, 4, 5, 2, 1},
                                    {5, 6, 8, 5, 3}});
    matrix<int> right = make_matrix({{3, 2, 4, 0},
                                     {3, 2, 4, 2},
                                     {0, 7, 3, 3},
                                     {3, 3, 0, 2},
                                     {1, 3, 2, 2}});
    REQUIRE(manhattan_tourist(down, right) == 34);

    down = make_matrix({{4, 4, 2, 1, 3, 1, 0, 0, 1},
                        {4, 3, 0, 2, 0, 4, 3, 4, 4},
                        {2, 3, 3, 1, 2, 1, 2, 2, 0},
                        {3, 0, 3, 3, 2, 1, 1, 3, 4},
                        {2, 3, 1, 2, 2, 0, 2, 3, 2},
                        {2, 2, 2, 0, 4, 2, 1, 0, 3},
                        {4, 3, 1, 1, 0, 1, 1, 4, 2},
                        {0, 0, 2, 2, 2, 1, 2, 4, 2},
                        {4, 3, 0, 3, 1, 3, 2, 3, 1},
                        {1, 4, 1, 0, 3, 4, 1, 2, 1},
                        {4, 4, 0, 4, 1, 4, 3, 1, 2},
                        {4, 1, 2, 3, 1, 3, 3, 3, 0},
                        {3, 1, 0, 2, 2, 0, 4, 4, 0},
                        {2, 0, 1, 0, 0, 3, 1, 1, 1},
                        {0, 1, 3, 2, 2, 2, 1, 2, 1},
                        {0, 2, 0, 3, 1, 2, 2, 4, 2},
                        {2, 0, 4, 1, 3, 3, 2, 4, 0},
                        {2, 3, 1, 3, 4, 2, 1, 4, 4}});


    right = make_matrix({{3, 3, 1, 1, 3, 4, 4, 4},
                         {4, 0, 3, 1, 0, 3, 4, 4},
                         {2, 2, 2, 3, 3, 1, 1, 4},
                         {1, 3, 1, 4, 4, 2, 0, 1},
                         {0, 2, 0, 3, 3, 3, 1, 0},
                         {3, 2, 0, 4, 1, 4, 4, 3},
                         {3, 0, 1, 1, 0, 3, 3, 0},
                         {3, 1, 1, 0, 2, 3, 4, 0},
                         {2, 4, 2, 1, 1, 3, 1, 2},
                         {1, 0, 4, 3, 0, 3, 3, 0},
                         {2, 3, 2, 4, 4, 3, 3, 0},
                         {3, 1, 2, 0, 3, 4, 3, 2},
                         {0, 0, 4, 1, 4, 0, 3, 4},
                         {3, 2, 3, 2, 0, 1, 2, 1},
                         {4, 3, 3, 2, 0, 1, 1, 2},
                         {0, 0, 4, 1, 2, 4, 0, 3},
                         {3, 4, 0, 1, 2, 3, 0, 1},
                         {4, 0, 2, 4, 2, 2, 4, 0},
                         {4, 3, 4, 2, 2, 3, 2, 3}});
    REQUIRE(manhattan_tourist(down, right) == 80);
}

int argmax(int a, int b, int c) {
    if (a >= b && a >= c)
        return 0;
    if (b >= a && b >= c)
        return 1;
    if (c >= a && c >= b)
        return 2;
    assert(false);
}

string longest_common_string(const string &string1, const string &string2) {
    enum Direction {
        undefined, right, down, diagonal
    };

    int size1 = size(string1);
    int size2 = size(string2);

    matrix<Direction> path(size1 + 1, size2 + 1);
    auto dp = make_matrix(size1 + 1, size2 + 1, std::numeric_limits<int>::max());
    // dp(0,0) = 0;
    // 'i' is row index (corresponding to string1) and 'j' is column index (corresponding to string2).

    dp(0, 0) = 0;
    path(0, 0) = undefined;

    for (int i = 1; i <= size1; ++i) {
        dp(i, 0) = 0;
        path(i, 0) = down;
    }

    for (int j = 1; j <= size2; ++j) {
        dp(0, j) = 0;
        path(0, j) = right;
    }

    for (int i = 1; i <= size1; ++i)
        for (int j = 1; j <= size2; ++j) {
            int the_argmax = argmax(dp(i - 1, j - 1), dp(i - 1, j), dp(i, j - 1));
            switch (the_argmax) {
                case 0: // Diagonal
                    dp(i, j) = (string1.at(i - 1) == string2.at(j - 1)) ? dp(i - 1, j - 1) + 1 : dp(i - 1, j - 1);
                    path(i, j) = diagonal;
                    break;
                case 1: // Going down
                    dp(i, j) = dp(i - 1, j);
                    path(i, j) = down;
                    break;
                case 2: // Going right
                    dp(i, j) = dp(i, j - 1);
                    path(i, j) = right;
                    break;
                default:
                    assert(false);
            }
        }

    int i = size1;
    int j = size2;

    string result;
    while (i != 0 or j != 0) {
        switch (path(i, j)) {
            case down:
                // result.insert(begin(result), string1.at(i-1));
                --i;
                break;
            case right:
                // result.insert(begin(result), string2.at(j-1));
                --j;
                break;
            case diagonal:
                if (string1.at(i - 1) == string2.at(j - 1))
                    result.insert(begin(result), string1.at(i - 1));
                --i;
                --j;
                break;
            case undefined:
                assert(false);
                break;
            default:
                assert(false);
        }
    }
    return result;
}

TEST_CASE("longest_common_string") {
    REQUIRE(longest_common_string("ACGATACGT", "CCCATTAAGT") == "CATAGT");
    REQUIRE(longest_common_string("ACGATACGT", "GACTATAGAA") == "ACATAG");
    REQUIRE(longest_common_string("CCCATTAAGT", "GACTATAGAA") == "CATAA");
    REQUIRE(longest_common_string("CTCGAT", "TACGTC") == "TCGT");
    REQUIRE(longest_common_string("ATGTTATA", "ATCGTCC") == "ATGT");
    REQUIRE(longest_common_string("AACCTTGG", "ACACTGTGA") == "AACTTG");

    string s1 = "AGCAGTTCCCTGATTGTTTAGTATTTGACTCCGTAGTTGAGCCTATATCGTAATTCTGCCAAGGAA";
    string s2 = "ATTATAATCCCCCGGACAAGAACCTAGTGGGCGCGTGGGACGGGACAAGAGCGACGTTCCGTAGGTGTGAGAGCCGTACTCAATTTTGGTTATTACCAT";
    REQUIRE(longest_common_string(s1, s2) == "AATTCCCGATTGTAGAGACTCCGTAGTTGAGCCTATATGTATTCCA");
}
