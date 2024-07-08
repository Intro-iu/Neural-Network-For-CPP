#include <iostream>
#include <cmath>
#include "Matrix.h"

using namespace std;

Matrix TargetFunction(const Matrix &x) {
    Matrix Mat = MatrixAssignment(x.row, x.col);
    Iter(Mat) Mat(i, j) = sin(x(i, j)+1);
    return Mat;
}

int main() {
    double l = -5;   double r = 5;
    int n = 750;
    int Input_Size = 1;
    int Output_Size = 1;

    Matrix Samples_X = MatrixAssignment(1, n, l, r);
    Matrix Samples_Y = TargetFunction(Samples_X);

    freopen("./Output/Train-Samples.txt", "w", stdout);
    disp(Samples_X);    disp(Samples_Y);   
    fclose(stdout);
    
    return 0;
}