#include <vector>
#include <limits>
#include <string>
#include <utility>
#include <unordered_map>
#include <set>
#include <boost/functional/hash.hpp>

#define CATCH_CONFIG_MAIN

#include "../include/catch.hpp"
#include "boost/numeric/ublas/matrix.hpp"

using std::vector;
using boost::numeric::ublas::matrix;
using boost::numeric::ublas::zero_matrix;
using boost::numeric::ublas::scalar_matrix;
using std::cout, std::endl;
using std::max;
using std::string;
using std::pair, std::make_pair;
using std::unordered_map;
using std::numeric_limits;
using std::set;

enum Direction {
    undefined, right, down, diagonal
};

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

int dp_change(int money, const vector<int> &coins) {
    vector<int> min_num_coins = vector<int>(money + 1);

    for (int m = 1; m <= money + 1; ++m) {
        min_num_coins[m] = numeric_limits<int>::max();
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
    int size1 = size(string1);
    int size2 = size(string2);

    matrix<Direction> path(size1 + 1, size2 + 1);
    auto dp = make_matrix(size1 + 1, size2 + 1, numeric_limits<int>::max());

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
                --i;
                break;
            case right:
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

pair<string, matrix<int>> get_blosum62() {
    string alphabet = "ACDEFGHIKLMNPQRSTVWY";
    auto the_matrix = make_matrix({{4,  0,  -2, -1, -2, 0,  -2, -1, -1, -1, -1, -2, -1, -1, -1, 1,  0,  0,  -3, -2},
                                   {0,  9,  -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2},
                                   {-2, -3, 6,  2,  -3, -1, -1, -3, -1, -4, -3, 1,  -1, 0,  -2, 0,  -1, -3, -4, -3},
                                   {-1, -4, 2,  5,  -3, -2, 0,  -3, 1,  -3, -2, 0,  -1, 2,  0,  0,  -1, -2, -3, -2},
                                   {-2, -2, -3, -3, 6,  -3, -1, 0,  -3, 0,  0,  -3, -4, -3, -3, -2, -2, -1, 1,  3},
                                   {0,  -3, -1, -2, -3, 6,  -2, -4, -2, -4, -3, 0,  -2, -2, -2, 0,  -2, -3, -2, -3},
                                   {-2, -3, -1, 0,  -1, -2, 8,  -3, -1, -3, -2, 1,  -2, 0,  0,  -1, -2, -3, -2, 2},
                                   {-1, -1, -3, -3, 0,  -4, -3, 4,  -3, 2,  1,  -3, -3, -3, -3, -2, -1, 3,  -3, -1},
                                   {-1, -3, -1, 1,  -3, -2, -1, -3, 5,  -2, -1, 0,  -1, 1,  2,  0,  -1, -2, -3, -2},
                                   {-1, -1, -4, -3, 0,  -4, -3, 2,  -2, 4,  2,  -3, -3, -2, -2, -2, -1, 1,  -2, -1},
                                   {-1, -1, -3, -2, 0,  -3, -2, 1,  -1, 2,  5,  -2, -2, 0,  -1, -1, -1, 1,  -1, -1},
                                   {-2, -3, 1,  0,  -3, 0,  1,  -3, 0,  -3, -2, 6,  -2, 0,  0,  1,  0,  -3, -4, -2},
                                   {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7,  -1, -2, -1, -1, -2, -4, -3},
                                   {-1, -3, 0,  2,  -3, -2, 0,  -3, 1,  -2, 0,  0,  -1, 5,  1,  0,  -1, -2, -2, -1},
                                   {-1, -3, -2, 0,  -3, -2, 0,  -3, 2,  -2, -1, 0,  -2, 1,  5,  -1, -1, -3, -3, -2},
                                   {1,  -1, 0,  0,  -2, 0,  -1, -2, 0,  -2, -1, 1,  -1, 0,  -1, 4,  1,  -2, -3, -2},
                                   {0,  -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0,  -1, -1, -1, 1,  5,  0,  -2, -2},
                                   {0,  -1, -3, -2, -1, -3, -3, 3,  -2, 1,  1,  -3, -2, -2, -3, -2, 0,  4,  -3, -1},
                                   {-3, -2, -4, -3, 1,  -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2},
                                   {-2, -2, -3, -2, 3,  -3, 2,  -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2,  7}});

    return make_pair<string, matrix<int>>("ACDEFGHIKLMNPQRSTVWY", std::move(the_matrix));
}

pair<string, matrix<int>> get_pam250() {
    auto the_matrix = make_matrix({{2,  -2, 0,  0,  -3, 1,  -1, -1, -1, -2, -1, 0,  1,  0,  -2, 1,  1,  0,  -6, -3},
                                   {-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4, 0,  -2, -2, -8, 0},
                                   {0,  -5, 4,  3,  -6, 1,  1,  -2, 0,  -4, -3, 2,  -1, 2,  -1, 0,  0,  -2, -7, -4},
                                   {0,  -5, 3,  4,  -5, 0,  1,  -2, 0,  -3, -2, 1,  -1, 2,  -1, 0,  0,  -2, -7, -4},
                                   {-3, -4, -6, -5, 9,  -5, -2, 1,  -5, 2,  0,  -3, -5, -5, -4, -3, -3, -1, 0,  7},
                                   {1,  -3, 1,  0,  -5, 5,  -2, -3, -2, -4, -3, 0,  0,  -1, -3, 1,  0,  -1, -7, -5},
                                   {-1, -3, 1,  1,  -2, -2, 6,  -2, 0,  -2, -2, 2,  0,  3,  2,  -1, -1, -2, -3, 0},
                                   {-1, -2, -2, -2, 1,  -3, -2, 5,  -2, 2,  2,  -2, -2, -2, -2, -1, 0,  4,  -5, -1},
                                   {-1, -5, 0,  0,  -5, -2, 0,  -2, 5,  -3, 0,  1,  -1, 1,  3,  0,  0,  -2, -3, -4},
                                   {-2, -6, -4, -3, 2,  -4, -2, 2,  -3, 6,  4,  -3, -3, -2, -3, -3, -2, 2,  -2, -1},
                                   {-1, -5, -3, -2, 0,  -3, -2, 2,  0,  4,  6,  -2, -2, -1, 0,  -2, -1, 2,  -4, -2},
                                   {0,  -4, 2,  1,  -3, 0,  2,  -2, 1,  -3, -2, 2,  0,  1,  0,  1,  0,  -2, -4, -2},
                                   {1,  -3, -1, -1, -5, 0,  0,  -2, -1, -3, -2, 0,  6,  0,  0,  1,  0,  -1, -6, -5},
                                   {0,  -5, 2,  2,  -5, -1, 3,  -2, 1,  -2, -1, 1,  0,  4,  1,  -1, -1, -2, -5, -4},
                                   {-2, -4, -1, -1, -4, -3, 2,  -2, 3,  -3, 0,  0,  0,  1,  6,  0,  -1, -2, 2,  -4},
                                   {1,  0,  0,  0,  -3, 1,  -1, -1, 0,  -3, -2, 1,  1,  -1, 0,  2,  1,  -1, -2, -3},
                                   {1,  -2, 0,  0,  -3, 0,  -1, 0,  0,  -2, -1, 0,  0,  -1, -1, 1,  3,  0,  -5, -3},
                                   {0,  -2, -2, -2, -1, -1, -2, 4,  -2, 2,  2,  -2, -1, -2, -2, -1, 0,  4,  -6, -2},
                                   {-6, -8, -7, -7, 0,  -7, -3, -5, -3, -2, -4, -4, -6, -5, 2,  -2, -5, -6, 17, 0},
                                   {-3, 0,  -4, -4, 7,  -5, 0,  -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2, 0,  10}});

    return make_pair<string, matrix<int>>("ACDEFGHIKLMNPQRSTVWY", std::move(the_matrix));
}

struct ScoredAlignment {
    int score;
    string alignment1;
    string alignment2;

    ScoredAlignment(int score,
                    const string &alignment1,
                    const string &alignment2) : score(score), alignment1(alignment1), alignment2(alignment2) {}
};

typedef unordered_map<pair<char, char>, int, boost::hash<pair<char, char>>> ScoringMatrix;

ScoringMatrix scoring_matrix_as_map(const matrix<int> &scoring_matrix, const string &alphabet) {
    ScoringMatrix res;
    for (unsigned long i = 0; i < alphabet.size(); ++i)
        for (unsigned long j = 0; j < alphabet.size(); ++j)
            res[pair<char, char>(alphabet.at(i), alphabet.at(j))] = scoring_matrix(i, j);

    return res;
}

ScoredAlignment best_scored_alignment(const string &string1,
                                      const string &string2,
                                      const matrix<int> &scoring_matrix,
                                      const string &alphabet,
                                      int sigma,
                                      bool local = false) {

    auto scores = scoring_matrix_as_map(scoring_matrix, alphabet);
    int size1 = size(string1);
    int size2 = size(string2);

    auto dp = make_matrix(size1 + 1, size2 + 1, numeric_limits<int>::max());
    matrix<pair<int, int>> previous(size1 + 1, size2 + 1);

    // Initialize all previous coordinates to (-1, -1), useful for debugging
    for (auto iter1 = previous.begin1(); iter1 != previous.end1(); ++iter1)
        for (auto iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
            *iter2 = make_pair(-1, -1);

    // 'i' is row index (corresponding to string1) and 'j' is column index (corresponding to string2).
    dp(0, 0) = 0;

    for (int i = 1; i <= size1; ++i) {
        // These are all free rides in case of local alignment, otherwise deletions
        dp(i, 0) = local ? 0 : dp(i - 1, 0) - sigma;
        previous(i, 0) = local ? make_pair(0, 0) : make_pair(i - 1, 0);
    }

    for (int j = 1; j <= size2; ++j) {
        // These are all free rides in case of local alignment, otherwise insertions
        dp(0, j) = local ? 0 : dp(0, j - 1) - sigma;
        previous(0, j) = local ? make_pair(0, 0) : make_pair(0, j - 1);
    }

    for (int i = 1; i <= size1; ++i)
        for (int j = 1; j <= size2; ++j) {
            // Cost of match/mismatch, based on the scoring matrix
            int cost_from_upper_left =
                    dp(i - 1, j - 1) + scores[pair<char, char>(string1.at(i - 1), string2.at(j - 1))];
            // Cost of deletion
            int cost_from_up = dp(i - 1, j) - sigma;
            // Cost of insertion
            int cost_from_left = dp(i, j - 1) - sigma;
            // Maximize the cost
            int the_argmax = argmax(cost_from_upper_left, cost_from_up, cost_from_left);
            // Check if a free-ride from the source would be preferable
            if (local && cost_from_left < 0 && cost_from_up < 0 && cost_from_upper_left < 0) {
                dp(i, j) = 0;
                previous(i, j) = make_pair(0, 0);
            } else
                switch (the_argmax) {
                    case 0: // Diagonal
                        dp(i, j) = cost_from_upper_left;
                        previous(i, j) = pair(i - 1, j - 1);
                        break;
                    case 1: // Going down
                        dp(i, j) = cost_from_up;
                        previous(i, j) = pair(i - 1, j);
                        break;
                    case 2: // Going right
                        dp(i, j) = cost_from_left;
                        previous(i, j) = pair(i, j - 1);
                        break;
                    default:
                        assert(false);
                }
        }

    // Allow a free-ride to the sink from any other vertex
    if (local) {
        for (unsigned long i = 0; i <= size1; ++i)
            for (unsigned long j = 0; j <= size2; ++j) {
                if (i == size1 && j == size2)
                    continue;
                if (dp(i, j) > dp(size1, size2)) {
                    dp(size1, size2) = dp(i, j);
                    previous(size1, size2) = make_pair(i, j);
                }
            }
    }

    // Backtracking of the solution, to find the two aligned strings.
    string alignment1, alignment2;
    int i = size1;
    int j = size2;
    while (i != 0 or j != 0) {
        int prev_i = previous(i, j).first;
        int prev_j = previous(i, j).second;
        if (prev_i < i && prev_j == j) {
            alignment1.insert(begin(alignment1), string1.at(i - 1));
            alignment2.insert(begin(alignment2), '-');
        } else if (prev_i == i && prev_j < j) {
            alignment2.insert(begin(alignment2), string2.at(j - 1));
            alignment1.insert(begin(alignment1), '-');
        } else if (prev_i < i && prev_j < j) {
            if (!local || (dp(prev_i, prev_j) != dp(i, j))) {
                alignment1.insert(begin(alignment1), string1.at(i - 1));
                alignment2.insert(begin(alignment2), string2.at(j - 1));
            }
        }
        i = prev_i;
        j = prev_j;
    }

    ScoredAlignment res(dp(size1, size2), alignment1, alignment2);
    return res;
}

int edit_distance(const string &string1, const string &string2) {
    const int sigma = 1;

    string all_characters = string1 + string2;
    set<char> characters(all_characters.begin(), all_characters.end());
    string alphabet(characters.begin(), characters.end());

    auto scoring_matrix = make_matrix(alphabet.size(), alphabet.size(), -1);
    for (unsigned long i = 0; i < alphabet.size(); ++i)
        for (unsigned long j = 0; j < alphabet.size(); ++j)
            if (alphabet.at(i) == alphabet.at(j))
                scoring_matrix(i, j) = 0;

    auto score_alignment = best_scored_alignment(string1,
                                                 string2,
                                                 scoring_matrix,
                                                 alphabet,
                                                 sigma,
                                                 false);

    return -score_alignment.score;
}

TEST_CASE("make_matrix") {
    matrix<int> down = make_matrix({{1, 0, 2, 4, 3},
                                    {4, 6, 5, 2, 1},
                                    {4, 4, 5, 2, 1},
                                    {5, 6, 8, 5, 3}});

    // By row...
    vector<int> by_row;
    for (auto iter1 = down.begin1(); iter1 != down.end1(); ++iter1)
        for (auto iter2 = iter1.begin(); iter2 != iter1.end(); ++iter2)
            by_row.emplace_back(*iter2);
    REQUIRE(by_row == vector<int>({1, 0, 2, 4, 3, 4, 6, 5, 2, 1, 4, 4, 5, 2, 1, 5, 6, 8, 5, 3}));

    // By column...
    vector<int> by_column;
    for (auto iter2 = down.begin2(); iter2 != down.end2(); ++iter2)
        for (auto iter1 = iter2.begin(); iter1 != iter2.end(); ++iter1)
            by_column.emplace_back(*iter1);
    REQUIRE(by_column == vector<int>({1, 4, 4, 5, 0, 6, 4, 6, 2, 5, 5, 8, 4, 2, 2, 5, 3, 1, 1, 3}));

    auto the_matrix = make_matrix(3, 4, -1);
    REQUIRE(the_matrix.size1() == 3);
    REQUIRE(the_matrix.size2() == 4);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            REQUIRE(the_matrix(i, j) == -1);

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


TEST_CASE("align") {
    REQUIRE(dp_change(40, {50, 25, 20, 10, 5, 1}) == 2);
    REQUIRE(dp_change(19415, {18, 16, 7, 5, 3, 1}) == 1080);
    REQUIRE(dp_change(16042, {16, 13, 11, 8, 7, 5, 3, 1}) == 1003);

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

TEST_CASE("best_scored_alignment") {
    // auto alphabet_scores = get_blosum62();
    auto[alphabet, scores] = get_blosum62();
    // auto alphabet = alphabet_scores.first;
    // auto scores = alphabet_scores.second;

    auto res = best_scored_alignment("PLEASANTLY", "MEANLY", scores, alphabet, 5);
    REQUIRE(res.score == 8);
    REQUIRE(res.alignment1 == "PLEASANTLY");
    REQUIRE(res.alignment2 == "-ME--AN-LY");

    res = best_scored_alignment(
            "KHLGRRPTYGFPFWYMVWDFQCQDDKEQKFFCKPRHVPCTWLGCEVTDEMWMDLHVEVQPQFCLVRQEFWHIFPPFSSIYWMYFDPSDVNRIMHDD",
            "KPTYGFPFWYMDWDFQCQDEWKKEIRFCKEQKFFCKPRHVPCWWLGCEVMDLHTQRYFWH", scores, alphabet, 5);
    REQUIRE(res.score == 43);
    REQUIRE(res.alignment1 ==
            "KHLGRRPTYGFPFWYMVWDFQCQD----D----KEQKFFCKPRHVPCTWLGCEVTDEMWMDLHVEVQPQFCLVRQEFWHIFPPFSSIYWMYFDPSDVNRIMHDD");
    REQUIRE(res.alignment2 ==
            "K-----PTYGFPFWYMDWDFQCQDEWKKEIRFCKEQKFFCKPRHVPCWWLGCEV-----MDLH--T--Q----R------Y--F----W------------H--");
}

int score_alignment(const string &s1,
                    const string &s2,
                    const pair<string, matrix<int>> &scoring_matrix,
                    const int sigma) {
    assert(s1.size() == s2.size());
    auto scores = scoring_matrix_as_map(scoring_matrix.second, scoring_matrix.first);
    int score = 0;
    for (unsigned long i = 0; i < s1.size(); ++i)
        score += (s1.at(i) == '-' || s2.at(i) == '-') ? -sigma : scores[make_pair(s1.at(i), s2.at(i))];
    return score;
}

TEST_CASE("best_scored_alignment (local)") {
    const auto alphabet_scores = get_pam250();
    const auto[alphabet, scores] = alphabet_scores;

    auto res = best_scored_alignment("AMMY", "PMMTN", scores, alphabet, 5, true);
    auto score = score_alignment(res.alignment1, res.alignment2, make_pair(alphabet, scores), 5);
    REQUIRE(res.score == 13);
    REQUIRE(score == res.score);
    REQUIRE(res.alignment1 == "AMM");
    REQUIRE(res.alignment2 == "PMM");

    res = best_scored_alignment("MEANLY", "PENALTY", scores, alphabet, 5, true);
    score = score_alignment(res.alignment1, res.alignment2, alphabet_scores, 5);
    REQUIRE(res.score == 15);
    REQUIRE(score == res.score);
    REQUIRE(res.alignment1 == "EL-Y");
    REQUIRE(res.alignment2 == "ELTY");

    res = best_scored_alignment(
            "AMTAFRYRQGNPRYVKHFAYEIRLSHIWLLTQMPWEFVMGIKMPEDVFQHWRVYSVCTAEPMRSDETYPCELFTVFDDIFTAEPVVCSCFYDDPM",
            "WQEKAVDGTVPSRHQYREKEDRQGNEIGKEFRRGPQVCEYSCNSHSCGWMPIFCIVCMSYVAFYCGLEYPMSRKTAKSQFIEWCDWFCFNHEFIPWVLRRYVVYDKIRYNYSYRNSASMEFV",
            scores, alphabet, 5, true);
    score = score_alignment(res.alignment1, res.alignment2, alphabet_scores, 5);
    REQUIRE(res.score == 56);
    REQUIRE(score == res.score);
    REQUIRE(res.alignment1 == "QGPYVKHFYEIRLHIWLLQMWEFVGIKMPE-VFQH---W-VYSVCEPMSDTYPCL-FVFDFAEPVV-C--CFYDD");
    REQUIRE(res.alignment2 == "DGV-SRH-Y--REEDRQGEIKEFRGPQVCEYCNSHSCGWMIF--C-CMYVFY-CLEYMSRASQFIEWCWFCFNHE");
}

TEST_CASE("edit_distance") {
    auto dist = edit_distance("PLEASANTLY", "MEANLY");
    REQUIRE(dist == 5);

    auto s1 = "GGACRNQMSEVNMWGCWWASVWVSWCEYIMPSGWRRMKDRHMWHWSVHQQSSPCAKSICFHETKNQWNQDACGPKVTQHECMRRRLVIAVKEE";
    auto s2 = "GMWGFVQVSTQSRFRHMWHWSVHQQSSECAKSICHHEWKNQWNQDACGPKVTQHECMANMPMHKCNNWFWRLVIAVKEEKVRETKMLDLIHRHWLVLNQGRMNEHNVTLRKSPCVKRIMHKWKSRTTFHR";
    dist = edit_distance(s1, s2);
    REQUIRE(dist == 97);

    dist = edit_distance("AC", "AC");
    REQUIRE(dist == 0);

    dist = edit_distance("AT", "G");
    REQUIRE(dist == 2);

    dist = edit_distance("CAGACCGAGTTAG", "CGG");
    REQUIRE(dist == 10);

    dist = edit_distance("CGT", "CAGACGGTGACG");
    REQUIRE(dist == 9);
}