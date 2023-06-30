#pragma once
#include <iostream>
#include <unistd.h>
#include "Matrix.h"
using namespace std;


struct NNwork {
    int L;  int *Siz;
    Matrix *Z, *A, *W, *B;
    Matrix *dZ, *dW, *dB;

    NNwork() {
        L = 1;
        Siz = new int[L];
        Z = new Matrix[L];
        A = new Matrix[L];
        W = new Matrix[L];
        B = new Matrix[L];
        dZ = new Matrix[L];
        dW = new Matrix[L];
        dB = new Matrix[L];
    }

    NNwork(const int &num, std::initializer_list <int> il, const double &l, const double &r) {
        L = num;
        Siz = new int[L];
        Z = new Matrix[L];
        A = new Matrix[L];
        W = new Matrix[L];
        B = new Matrix[L];
        dZ = new Matrix[L];
        dW = new Matrix[L];
        dB = new Matrix[L];

        int cnt = -1;
        for(auto itm : il)  Siz[++cnt] = itm;
        for (int i = 1;i < L;i++) {
            W[i] = MatrixAssignment(Siz[i], Siz[i-1], l, r);
            B[i] = MatrixAssignment(Siz[i], 1, l, r);
        }
    }

    NNwork(const int &num, int arr[], const double &l, const double &r) {
        L = num;
        Siz = new int[L];
        Z = new Matrix[L];
        A = new Matrix[L];
        W = new Matrix[L];
        B = new Matrix[L];
        dZ = new Matrix[L];
        dW = new Matrix[L];
        dB = new Matrix[L];

        for (int i = 0;i < L;i++)   Siz[i] = arr[i];

        for (int i = 1;i < L;i++) {
            W[i] = MatrixAssignment(Siz[i], Siz[i-1], l, r);
            B[i] = MatrixAssignment(Siz[i], 1, l, r);
        }
    }

    ~NNwork(){
        delete [] Z, A, W, B, dZ, dW, dB;
    }

    Matrix Act(const Matrix &x) {
        return logistic(x);
    }

    Matrix Act_D(const Matrix &x) {
        return logistic_D(x);
    }



    Matrix ForwardPropagation(const Matrix &x) {
        A[0] = Z[0] = x;
        for (int l = 1; l < L; l++) {
            Matrix B_pi = B[l], tmp = MatrixAssignment(W[l].row, A[l-1].col);
            BroadCast(tmp, B_pi);
            Z[l] = W[l] * A[l-1] + B_pi;
            A[l] = Act(Z[l]);
        }
        //std::cout << "F: ";  disp(A[L-1]);
        return A[L-1];
    }

    void BackwardPropagation(const Matrix &y) {
        dZ[L-1] = 2 % (log(1 + A[L-1]) - log(1 + y)) / (1 + A[L-1]);

        //dZ[L-1] = 2 % (A[L-1] - y) % Act_D(Z[L-1]);

        for (int l = L-1; l > 0; l--) {
            dW[l] = dZ[l] * T(A[l-1]);
            dB[l] = sum(dZ[l], 1);
            if(l-1) dZ[l-1] = (T(W[l]) * dZ[l]) % Act_D(Z[l-1]);
        }
    }

    double MSLE(const Matrix &y) {
        return sumAll((log(1 + y) - log(1 + A[L-1])) % (log(1 + y) - log(1 + A[L-1]))) / y.row / y.col;
    }

    double MSE(const Matrix &y) {
        return sumAll((y - A[L-1]) % (y - A[L-1])) / y.row / y.col;
    }

    void GradientDescent(const double &alpha) {
        for (int i = 1;i < L;i++) {
            W[i] = W[i] - (alpha % dW[i]);
            B[i] = B[i] - (alpha % dB[i]);
        }
    }
    
    void Output(const string &name, Matrix &Samples_X_ori, Matrix &Samples_Y_ori) {
        int o = dup(fileno(stdout));
        freopen(name.c_str(), "w", stdout);
        cout << L << endl;
        for(int i = 0;i < L;i++)    cout << Siz[i] << " ";
        cout << endl << endl;
        for(int i = 1; i < L;i++)   disp(W[i]);
        for(int i = 1; i < L;i++)   disp(B[i]);
        cout << Samples_X_ori.row << " " << Samples_X_ori.col << endl;  disp(Samples_X_ori);
        cout << endl;
        cout << Samples_Y_ori.row << " " << Samples_Y_ori.col << endl;  disp(Samples_Y_ori);
        dup2(o, fileno(stdout));
        close(o);
    }
};