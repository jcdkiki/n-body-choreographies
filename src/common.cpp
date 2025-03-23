#include "common.hpp"
#include <assert.h>

void M_InverseMatrix(double *mat, int n, double *inv)
{
    double *aug = (double *)malloc(n * 2 * n * sizeof(double));
    assert(aug != NULL);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * 2 * n + j] = mat[i * n + j];
        }
        for (int j = n; j < 2 * n; j++) {
            aug[i * 2 * n + j] = (j == i + n) ? 1.0 : 0.0;
        }
    }

    for (int c = 0; c < n; c++) {
        int pivot_row = c;
        double max_val = fabs(aug[c * 2 * n + c]);
        for (int i = c + 1; i < n; i++) {
            double val = fabs(aug[i * 2 * n + c]);
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        if (fabs(max_val) < 1e-10) {
            free(aug);
            return;
        }

        if (pivot_row != c) {
            for (int j = 0; j < 2 * n; j++) {
                double temp = aug[c * 2 * n + j];
                aug[c * 2 * n + j] = aug[pivot_row * 2 * n + j];
                aug[pivot_row * 2 * n + j] = temp;
            }
        }

        double pivot = aug[c * 2 * n + c];
        for (int j = 0; j < 2 * n; j++) {
            aug[c * 2 * n + j] /= pivot;
        }

        for (int r = 0; r < n; r++) {
            if (r == c) continue;
            double factor = aug[r * 2 * n + c];
            for (int j = 0; j < 2 * n; j++) {
                aug[r * 2 * n + j] -= factor * aug[c * 2 * n + j];
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i * n + j] = aug[i * 2 * n + (n + j)];
        }
    }

    free(aug);
}
