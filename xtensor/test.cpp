#include "inversion.hpp"
#include <fstream>
#include <time.h>
#include <unistd.h>
#include "xtensor/xrandom.hpp"


//Multiply
mat dot(const mat &A, const mat &B) {
    mat res = xt::zeros<double>({N, N});
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                res[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
    for (int i = 0; i < N * N; i++) {
        if (fabs(res[i]) < zero) {
            res[i] = 0.0;
        }
    }
    return res;
}

void outStream(std::ostream &out, const mat &A, const mat &inv_A) {
    out << "A:" << std::endl;
    out << A << std::endl;
    out << "inv_A:" << std::endl;
    out << inv_A << std::endl;
    out << "A*inv_A:" << std::endl;
    out << dot(A, inv_A) << std::endl;
}

//Multiply the inverse of the matrix and matrix to determine whether it is a unit matrix
bool isCorrect(mat A, mat inv_A) {
    mat single = xt::zeros<int>({N, N});
    if (inv_A == single) {
        //In "Inversion::solve" function,if we met single matrix ,we will return zero matrix.
        return true;
    }
    mat res = dot(A, inv_A);
    mat eye = xt::eye<double>(N);
    //If the absolute value is less than 10^-10, it means equal
    for (int i = 0; i < N * N; i++) {
        if (fabs(res[i] - eye[i]) > zero) {
            return false;
        }
    }
    return true;
}

//run n times and record time
void test(int n) {
    Inversion inv;
    mat A_rand, inv_A;
    mat eye = xt::eye(N);
    clock_t start, totalTime;
    for (int i = 0; i < n; ++i) {
        A_rand = xt::random::randn<double>({N, N});
        start = clock();
        inv_A = inv.solve(A_rand);
        totalTime += clock() - start;
        if (!isCorrect(A_rand, inv_A)) {
            std::cout << "Wrong Answer" << std::endl;
            outStream(std::cout, A_rand, inv_A);
            return;
        }
    }
    std::cout << "Accept" << std::endl;
    std::cout << "Run " << n << " times." << std::endl;
    std::cout << "Total time:" << (double) totalTime / CLOCKS_PER_SEC << " s" << std::endl;
}

//Read a matrix from file and output its inverse to another file
void runFromFile(std::string input, std::string output) {
    std::ifstream inFile(input.c_str(), std::ios::in);
    if (!inFile)
        std::cout << "input error" << std::endl;

    double tmp;
    mat A = xt::zeros<double>({N, N});
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            inFile >> tmp;
            A[i * N + j] = tmp;
        }
    }

    Inversion inv;
    mat inv_A = inv.solve(A);
    std::ofstream outFile(output.c_str(), std::ios::out);
    outStream(outFile, A, inv_A);
    std::cout<<"Run successfully, please open the file to view"<<std::endl;
}

int main(int argc, char **argv) {
    int opt, times = -1;
    std::string input = "";
    std::string output = "";
    while ((opt = getopt(argc, argv, "i:o:t:")) != -1) {
        switch (opt) {
            case 'i':
                input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 't':
                times = atoi(optarg);
                break;
        }
    }
    if (times != -1) {
        test(times);
    } else if (input != "" && output != "") {
        runFromFile(input, output);
    } else {
        fprintf(stderr, "usage: %s -i [input path] -o [output path] -t [times] \n", argv[0]);
        fprintf(stderr, "usage: %s -t [times] \n", argv[0]);
        exit(-1);
    }
    return 0;
}