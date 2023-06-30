#include <iostream>
#include <cmath>
#include "Matrix.h"

using namespace std;

Matrix TargetFunction(const Matrix &x) {
    Matrix Mat = MatrixAssignment(x.row, x.col);
    Iter(Mat) Mat(i, j) = x(i, j)*x(i, j);
    return Mat;
}

int main() {
    double l = 0;   double r = 6;
    int n = 5000;
    int Input_Size = 1;
    int Output_Size = 1;

    Matrix Samples_X = MatrixAssignment(1, n, l, r);
    Matrix Samples_Y = TargetFunction(Samples_X);

    freopen("./Output/Train-Samples.txt", "w", stdout);
    disp(Samples_X);    disp(Samples_Y);   
    fclose(stdout);
    
    return 0;
}