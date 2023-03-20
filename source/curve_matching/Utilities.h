/*-----------------------------------------------------------------------*\
|     ____            _  ______ __  __  ____  _  ________                 |
|    / __ \       _  (_)/  ___ |  \/  |/ __ \| |/ /  ____|                |
|   | |  | |_ __ | |_ _|  (___ | \  / | |  | | ' /| |__    _     _        |
|   | |  | | '_ \|  _| |\___  \| |\/| | |  | |  < |  __| _| |_ _| |_      |
|   | |__| | |_) | |_| |____)  | |  | | |__| | . \| |___|_   _|_   _|     |
|    \____/| .__/\___|_|______/|_|  |_|\____/|_|\_\______||_|   |_|       |
|          | |                                                            |
|          |_|                                                            |
|                                                                         |
|            Authors: Magnus Fürst <magnus.furst@ulb.ac.be>               |
|                     Andrea Bertolino <andrea.bertolino@ulb.be>          |
|-------------------------------------------------------------------------|
|   License                                                               |
|                                                                         |
|   This file is part of OptiSMOKE.                                       |
|   Copyright (C) 2020 by Magnus Fürst and Andrea Bertolino               |
|                                                                         |
|   OptiSMOKE is free software: you can redistribute it and/or modify     |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OptiSMOKE is distributed in the hope that it will be useful,          |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OptiSMOKE. If not, see <http://www.gnu.org/licenses/>.     |
|                                                                         |
\*-----------------------------------------------------------------------*/

/* Prints a vector */
auto printV = [](const auto& std::vector) {

    for (int i=0; i<std::vector.size(); ++i)
        cout << "(" << i << ") " << std::vector[i] << endl;
    cout << endl;

};



/* Prints a matrix */
auto printM = [](const auto& matrix) {

    for (int i=0; i<matrix.size(); ++i) {
        for (int j=0; j<matrix[i].size(); ++j)
            cout << "(" << i << "," << j << ") " << matrix[i][j] << "\t    ";
        cout << endl;
    }
    cout << endl;

};



/* Prints the matrix equal to the difference of the two input matrices */
void print(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {

    for (int i = 0; i<A.size(); ++i) {
        for (int j = 0; j<A[i].size(); ++j)
            cout << "(" << i << "," << j << ") " << A[i][j]-B[i][j] << "\t    ";
        cout << endl;
    }
    cout << endl;

}



/* Prints the vector equal to the difference of the two input vectors */
void print(const std::vector<double>& a, const std::vector<double>& b) {

    for (int i=0; i<a.size(); ++i)
        cout << "(" << i << ") " << a[i]-b[i] << endl;
    cout << endl;

}



/* Prints the GCV1 vector and the corresponding log10lambda values */
void printGCV1(const std::vector<double>& GCV1,
               double log10lambdaMin,
               double log10lambdaStep) {

    for (int i=0; i<GCV1.size(); ++i) {
        cout << "(" << log10lambdaMin+i*log10lambdaStep << ")\t " << GCV1[i];
        if (i != GCV1.size()-1) {
            if (GCV1[i] > GCV1[i+1]) cout << "   \t ++";
            else if (GCV1[i] < GCV1[i+1]) cout << "   \t--";
            else cout << "   \tequal!";
        }
        cout << endl;
    }
    cout << endl;

}



/* Finds the minimum and the maximum value of a matrix */
void minMax(const std::vector<std::vector<double>>& matrix) {

    double min = matrix[0][0];
    int rowMin = 0;
    int columnMin = 0;
    double max = matrix[0][0];
    int rowMax = 0;
    int columnMax = 0;
    for (int i=0; i<matrix.size(); ++i) {
        for (int j=0; j<matrix[i].size(); ++j) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                rowMin = i;
                columnMin = j;
            }
            if (matrix[i][j] > max) {
                max = matrix[i][j];
                rowMax = i;
                columnMax = j;
            }
        }
    }
    cout << "Min: (" << rowMin << "," << columnMin << ") " << min << endl
         << "Max: (" << rowMax << "," << columnMax << ") " << max << endl
         << endl;

}



/* Finds the minimum and the maximum value of the difference of the two input
matrices */
void minMax(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {

    auto matrix = std::vector<std::vector<double>>(A.size(),std::vector<double>(A[0].size()));
    for (int i=0; i<A.size(); ++i)
        for (int j=0; j<A[0].size(); ++j)
            matrix[i][j] = A[i][j] - B[i][j];

    double min = matrix[0][0];
    int rowMin = 0;
    int columnMin = 0;
    double max = matrix[0][0];
    int rowMax = 0;
    int columnMax = 0;
    for (int i=0; i<matrix.size(); ++i) {
        for (int j=0; j<matrix[i].size(); ++j) {
            if (matrix[i][j] < min) {
                min = matrix[i][j];
                rowMin = i;
                columnMin = j;
            }
            if (matrix[i][j] > max) {
                max = matrix[i][j];
                rowMax = i;
                columnMax = j;
            }
        }
    }
    cout << "Min: (" << rowMin << "," << columnMin << ") " << min << endl
         << "Max: (" << rowMax << "," << columnMax << ") " << max << endl
         << endl;

}



/* Calculates the inverse of matrix A using the Gauss-Jordan method. Modifies
matrix A while running */
void invertWithGaussJordan(std::vector<std::vector<double>>& A,
                           std::vector<std::vector<double>>& I) {

    int K = A.size();

    // Matrix I starts as an identity matrix with the same dimensions as A
    I = std::vector<std::vector<double>>(K,std::vector<double>(K));
    for (int i=0; i<K; ++i)
        I[i][i] = 1.;

    for (int i=0; i<K; ++i) {
        double alfa = A[i][i];
        for (int j=0; j<K; ++j) {
            A[i][j] /= alfa;
            I[i][j] /= alfa;
        }
        for (int k=0; k<K; ++k) {
            if (k != i && A[k][i] != 0) {
                double beta = A[k][i];
                for (int j=0; j<K; ++j) {
                    A[k][j] -= beta * A[i][j];
                    I[k][j] -= beta * I[i][j];
                }
            }
        }
    }

}
