#include <iostream>
#include "Matrix.h"
using namespace std;

int main() {
    Matrix Mat = MatrixAssignment(5, 1);
    int t = 0;
    Iter(Mat) {
        Mat(i, j) = t++;
    }
    disp(Mat);
    disp(Normalize_ZScore(Mat));
    disp(Restore_ZScore(Mat, Normalize_ZScore(Mat)));
}

