#include <iostream>
#include "Matrix.h"
#include "NeuralNetwork.h"

using namespace std;

int L, row, col;
int Sample_x_row, Sample_x_col, Sample_y_row, Sample_y_col;
Matrix Input, Sample_x, Sample_y;

int main() {
    freopen("./Output/Train-Para.txt", "r", stdin);
    cin >> L;
    int Siz[L];
    for(int i = 0;i < L;i++)    cin >> Siz[i];
    NNwork Network(L, Siz, 0, 0);
    for(int l = 1;l < L;l++)    Iter(Network.W[l]) cin >> Network.W[l](i, j);
    for(int l = 1;l < L;l++)    Iter(Network.B[l]) cin >> Network.B[l](i, j);

    cin >> Sample_x_row >> Sample_x_col;
    Sample_x = MatrixAssignment(Sample_x_row, Sample_x_col);
    Iter(Sample_x) cin >> Sample_x(i, j);

    cin >> Sample_y_row >> Sample_y_col;
    Sample_y = MatrixAssignment(Sample_y_row, Sample_y_col);
    Iter(Sample_y) cin >> Sample_y(i, j);

    fclose(stdin);
    freopen("./Output/Test-Input.txt", "r", stdin);
    cin >> row >> col;
    Input = MatrixAssignment(row, col);
    Iter(Input) cin >> Input(i, j);

    fclose(stdin);

    Network.ForwardPropagation(Normalize_MinMax(Input));
    freopen("./Output/Test-Output.txt", "w", stdout);
    disp(Restore_MinMax(Sample_y, Network.A[L-1]));
    fclose(stdout);

}