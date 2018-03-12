#ifndef INVERSION_H
#define INVERSION_H

#include "xtensor/xarray.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xindex_view.hpp"

typedef xt::xarray<double> mat;
const int N = 5;
const double zero = pow(10, -9);

//Use LUP decomposition to perform inverse calculation
class Inversion {
private:
    mat L, U, inv_A, B, b, x, y, diag, A_copy;
    xt::xarray<int> P;

    bool LUP_decomposition(const mat &A);

    xt::xarray<double> LUP_solve();

public:
    Inversion();

    mat solve(const mat &A);

};

Inversion::Inversion() {//init
    A_copy = xt::zeros<double>({N, N});
    inv_A = xt::zeros<double>({N, N});

    x = xt::zeros<double>({N});
    y = xt::zeros<double>({N});
    diag = xt::zeros<double>({N});
    L = xt::zeros<double>({N, N});
    U = xt::zeros<double>({N, N});

    B = xt::eye<int>(N);
    P = xt::zeros<int>({N});
    b = xt::zeros<int>({N});
}

mat Inversion::solve(const mat &A) {
    bool isSingle = LUP_decomposition(A);
    if (isSingle) {//if we met single matrix ,we will return zero matrix.
        std::cout << "This is a single matrix" << std::endl;
        return xt::zeros<int>({N, N});
    }
    for (int i = 0; i < N; ++i) {
        b = xt::view(B, i);
        xt::view(inv_A, i) = LUP_solve();
    }
    inv_A = xt::transpose(inv_A);
    return inv_A;
}

bool Inversion::LUP_decomposition(const mat &A) {
    int row = 0;
    double p;
    P = xt::arange(N), A_copy = A;
    for (int i = 0; i < N - 1; i++) {
        p = 0.0d;
        for (int j = i; j < N; j++) {
            if (fabs(A_copy[j * N + i]) > p) {
                p = fabs(A_copy[j * N + i]);
                row = j;
            }
        }
        if (0 == p)//is single
            return true;

        std::swap(P[i], P[row]);
        for (int j = 0; j < N; j++)
            std::swap(A_copy[i * N + j], A_copy[row * N + j]);

        //It is the same with LU.
        double u = A_copy[i * N + i], l = 0.0d;
        for (int j = i + 1; j < N; j++) {
            l = A_copy[j * N + i] / u;
            A_copy[j * N + i] = l;
            for (int k = i + 1; k < N; k++)
                A_copy[j * N + k] = A_copy[j * N + k] - A_copy[i * N + k] * l;
        }
    }

    //Construct L and U.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            if (i != j) {
                L[i * N + j] = A_copy[i * N + j];
            } else {
                L[i * N + j] = 1;
            }
        }
        for (int k = i; k < N; k++) {
            U[i * N + k] = A_copy[i * N + k];
        }
    }
    //det
    double det = 1.0;
    diag = xt::diagonal(A_copy);
    for (int i = 0; i < N; ++i)
        det *= diag[i];
    return fabs(det) < zero;//is single

}

xt::xarray<double> Inversion::LUP_solve() {
    //Forward Substitution
    for (int i = 0; i < N; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] = y[i] - L[i * N + j] * y[j];
        }
    }
    //Back Substitution
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = N - 1; j > i; j--) {
            x[i] = x[i] - U[i * N + j] * x[j];
        }
        x[i] /= U[i * N + i];
    }
    return x;
}
#endif
