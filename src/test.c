#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

#include "include/row_echelon_integer.h"
#include "rgen.h"

int main() {
    int m = 5;
    int n = 18;
    int a[5 * 18] = { //Ag C Cl H N     0 3 12 14 16 17 
       1, 44,  0, 28,  4 ,
       0,  6,  0,  6,  0,
       0,  1,  0,  4,  0,
       0,  0,  1,  1,  0,
       0,  8,  0, 10,  0,
       0,  4,  0,  5,  1,
       0,  1,  0,  5,  1,
       0,  0,  0,  3,  1,
       0,  2,  0,  7,  1,
       0,  7,  0,  8,  0,
       0,  3,  0,  6,  0,
       0,  2,  0,  4,  0,
       0,  2,  0,  6,  0,
       0,  3,  0,  8,  0,
       1,  0,  1,  0,  0,
       0,  5,  0,  5,  1,
       0,  6,  1,  5,  0,
       0,  6,  0,  7,  1 };

    i4mat_print(m, n, a, "  Input A:");

    int* S = calloc(1000 * n, sizeof(int));
    int len_S = 0;
    generate(S, &len_S, a, m, n, 1000);

    printf("%d\n", len_S);
    for (int i = 0; i < len_S; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%d  ", S[i * n + j]);
        }
        printf("\n");
    }
    free(S);

    getchar();

    return 0;
}