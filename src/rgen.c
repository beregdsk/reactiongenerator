#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "include/row_echelon_integer.h"
#include "rgen.h"

bool contains_zero(int* x, int n) {
    for (int i = 0; i < n; ++i) {
        if (x[i] == 0) return true;
    }
    return false;
}

int gcd(int a, int b) { return b == 0 ? a : gcd(b, a % b); }

int lcm(int a, int b) { return a * b / gcd(a, b); }

int* minimal_nullspace(int* A, int m, int n, bool* found) {
    i4mat_rref(m, n, A);

    int n_pivots = 0;
    int j = 0;
    for (int i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            if (A[j * m + i] != 0) {
                ++n_pivots;
                break;
            }
        }

        if (j == n) break;
    }

    int n_null = n - n_pivots;
    *found = n_null > 0;

    if (n_null == 1) {
        int* null_vector = malloc(sizeof(int) * n);

        int l = 1;
        for (int i = 0; i < n - 1; ++i) {
            l = lcm(l, A[i * m + i]);
        }

        for (int i = 0; i < n - 1; ++i) {
            null_vector[i] = A[(n - 1) * m + i] * l / A[i * m + i];
        }
        null_vector[n - 1] = -l;

        if (!contains_zero(null_vector, n)) {
            return null_vector;
        }
        free(null_vector);
    }
    return 0;
}

DLL00_EXPORT_API void generate(int* S, int* len_S, int* ref, int dim, int n_vecs, int maxn) {
    memset(S, 0, sizeof S);

    int* vec_nums = malloc(sizeof(int) * n_vecs);
    vec_nums[0] = 0;
    vec_nums[1] = 1;
    int len_nums = 2;

    while (len_nums > 1 && *len_S < maxn) {
        int* A = malloc(sizeof(int) * len_nums * dim);

        for (int i = 0; i < len_nums; ++i) {
            for (int j = 0; j < dim; ++j) {
                A[i * dim + j] = ref[vec_nums[i] * dim + j];
            }
        }

        bool found = false;
        int* null_vector = minimal_nullspace(A, dim, len_nums, &found);
        
        free(A);

        if (found) {
            if (null_vector != 0) {
                for (int i = 0; i < len_nums; ++i) {
                    S[(*len_S) * n_vecs + vec_nums[i]] = null_vector[i];
                }
                ++(*len_S);
            }

            int last = vec_nums[--len_nums];
            if (last + 1 >= n_vecs) {
                last = vec_nums[--len_nums];
            }
            vec_nums[len_nums++] = last + 1;
        }
        else {
            int last = vec_nums[len_nums - 1];
            if (last + 1 < n_vecs) {
                vec_nums[len_nums++] = last + 1;
            }
            else {
                last = vec_nums[--len_nums];

                int j = 0;
                while (last == n_vecs - j && len_nums > 1) {
                    last = vec_nums[--len_nums];
                    ++j;
                }
                vec_nums[len_nums++] = last + 1;
            }
        }

        free(null_vector);
    }

    free(vec_nums);
}