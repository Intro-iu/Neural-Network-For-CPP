#include <iostream>
#include <cstdlib>
#include "Matrix.h"
#include "NeuralNetwork.h"

using namespace std;

void ReadSamples(const int &n, const int &Input_Size, const int &Output_Size, Matrix &Samples_X, Matrix &Samples_Y) {
    Samples_X = MatrixAssignment(Input_Size, n);
    Samples_Y = MatrixAssignment(Output_Size, n);
    
    freopen("./Output/Train-Samples.txt", "r", stdin);
    Iter(Samples_X) cin >> Samples_X(i, j);
    Iter(Samples_Y) cin >> Samples_Y(i, j);
    fclose(stdin);

}


int main() {
    int n = 200;
    int Generation = 100000;

    Matrix Samples_X_ori, Samples_Y_ori;
    ReadSamples(n, 1, 1, Samples_X_ori, Samples_Y_ori);

    //Matrix Samples_X_nor = Samples_X_ori;
    //Matrix Samples_Y_nor = Samples_Y_ori;

    Matrix Samples_X_nor = Normalize_MinMax(Samples_X_ori);
    Matrix Samples_Y_nor = Normalize_MinMax(Samples_Y_ori);

    //cout << "Samples:\n";   disp(Samples_X_nor);    disp(Samples_Y_nor);
    
    NNwork Network(4, {1, 2, 3, 1}, -2, 2);
    
    long long t = 0;

    //while (Network.CostFunction(Samples_Y_nor) >= 1.0e-10) {
    while (t < Generation) {
        Network.ForwardPropagation(Samples_X_nor);
        Network.BackwardPropagation(Samples_Y_nor);
        Network.GradientDescent(0.001);
        if(t++ % 5000 == 0 || t == Generation) {
            Network.Output("./Output/Train-Para.txt", Samples_X_ori, Samples_Y_ori);
            cout << "Generation: " << t << "  MSLE: " << Network.MSLE(Samples_Y_nor) << "  MSE: " << Network.MSE(Samples_Y_nor) << endl;
        }
    }
}
