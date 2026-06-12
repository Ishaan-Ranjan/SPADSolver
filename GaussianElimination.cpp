#include <iostream>
#include <vector>
#include "GaussianElimination.h"

void rowEchelon(std::vector<double>& matrix, int const bins){
    int const columns = bins;
    double mult_factor; 

    //Makes Matrix Upper-Triangular
    for(int i = 0; i < columns-1; i++){
        mult_factor = matrix[(i + 1) * (columns + 1) + i]/matrix[i * (columns + 1) + i];

        matrix[(i + 1) * (columns + 1) + i] = 0;
        matrix[(i + 1) * (columns + 1 ) + bins-1] -= mult_factor * matrix[i * (columns + 1) + bins-1];
        matrix[(i + 1) * (columns + 1 ) + bins] -= mult_factor * matrix[i * (columns + 1) + bins];
    }

    //Diagonalizes

    for(int i = 0; i < columns - 1; i++){
        mult_factor = matrix[i * (columns + 1) + bins-1]/matrix[(bins-1) * (columns + 1) + bins-1];

        matrix[i * (columns + 1) + bins] -= mult_factor*matrix[(bins-1) * (columns + 1) + bins];
        matrix[i * (columns + 1) + bins-1] = 0;
    }

    //Sets Coefficients to 1

    for(int i = 0; i < columns; i++){
        matrix[i * (columns + 1) + bins] = matrix[i * (columns + 1) + bins]/matrix[i * (columns + 1) + i];
        matrix[i * (columns + 1) + i] = 1;
    }
   /*for(int i = 0; i < matrix.size(); i++){
    std::cout << matrix[i] << " ";
    }*/
}
std::vector<double> solVector(std::vector<double>& matrix, int const bins){
    std::vector<double> sol_vec(bins);

    for(int i = 0; i < bins; i++){
        sol_vec[i] = matrix[i * (bins + 1) + bins];
    }
    return sol_vec;
}