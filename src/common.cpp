#include "common.hpp"
#include <assert.h>

/*
void M_GetCof(double *mat, double *cof, int p, int q, int n)
{
    int i = 0, j = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (row != p && col != q) {
                cof[i*(n-1) + j++] = mat[row*n + col];
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double M_GetDeterminant(double *mat, int n)
{
    if (n == 1)
        return mat[0];
    
    double det = 0;
    
    double *cof = (double*)alloca(sizeof(double) * (n-1) * (n-1));
    double sign = 1.0;
    for (int f = 0; f < n; f++) {
        M_GetCof(mat, cof, 0, f, n);
        det += sign * mat[f] * M_GetDeterminant(cof, n - 1);
        sign = -sign;
    }
    return det;
}

void M_CalculateAdjoint(double *mat, int n, double *adj)
{
    if (n == 1) {
        adj[0] = 1.0;
        return;
    }
    
    double *cof = (double*)alloca(sizeof(double) * (n-1)*(n-1));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            M_GetCof(mat, cof, i, j, n);
            double sign = ((i + j) % 2 == 0) ? 1 : -1;
            adj[j*n + i] = sign * M_GetDeterminant(cof, n - 1);
        }
    }
}

void M_InverseMatrix(double *mat, int n, double *inv)
{
    double det = M_GetDeterminant(mat, n);
    assert(det != 0);

    M_CalculateAdjoint(mat, n, inv);

    for (int i = 0; i < n*n; i++) {
        inv[i] /= det;
    }
}
*/


void M_InverseMatrix(double *mat, int n, double *inv)
{
    double *aug = (double *)malloc(n * 2 * n * sizeof(double));
    assert(aug != NULL);

    // Initialize augmented matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aug[i * 2 * n + j] = mat[i * n + j];
        }
        for (int j = n; j < 2 * n; j++) {
            aug[i * 2 * n + j] = (j == i + n) ? 1.0 : 0.0;
        }
    }

    for (int c = 0; c < n; c++) {
        // Find pivot row
        int pivot_row = c;
        double max_val = fabs(aug[c * 2 * n + c]);
        for (int i = c + 1; i < n; i++) {
            double val = fabs(aug[i * 2 * n + c]);
            if (val > max_val) {
                max_val = val;
                pivot_row = i;
            }
        }

        // Check for singularity
        if (fabs(max_val) < 1e-10) {
            free(aug);
            return;
        }

        // Swap rows
        if (pivot_row != c) {
            for (int j = 0; j < 2 * n; j++) {
                double temp = aug[c * 2 * n + j];
                aug[c * 2 * n + j] = aug[pivot_row * 2 * n + j];
                aug[pivot_row * 2 * n + j] = temp;
            }
        }

        // Normalize pivot row
        double pivot = aug[c * 2 * n + c];
        for (int j = 0; j < 2 * n; j++) {
            aug[c * 2 * n + j] /= pivot;
        }

        // Eliminate other rows
        for (int r = 0; r < n; r++) {
            if (r == c) continue;
            double factor = aug[r * 2 * n + c];
            for (int j = 0; j < 2 * n; j++) {
                aug[r * 2 * n + j] -= factor * aug[c * 2 * n + j];
            }
        }
    }

    // Copy the inverse part
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv[i * n + j] = aug[i * 2 * n + (n + j)];
        }
    }

    free(aug);
}
