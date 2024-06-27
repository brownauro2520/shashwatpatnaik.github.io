#ifndef VECTORALGEBRA_H_INCLUDED
#define VECTORALGEBRA_H_INCLUDED

vector<double> add_v(vector<double> vec1, vector<double> vec2) {
    int size = vec1.size();
    vector<double> result(size);
    for (int i = 0; i < size; i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

double sum_V(vector<vector<double>> v) {
  double sum = 0;
  for (int i = 0; i < v.size(); i++) {
    for (int j = 0; j < v[i].size(); j++) {
      sum += abs(v[i][j]);
    }
  }
  return sum;
}
vector<double> multiplyVS(double scalar, vector<double> v) {
    vector<double> result;

    for (int i = 0; i < v.size(); i++) {
        double product = scalar * v[i];
        result.push_back(product);
    }

    return result;
}


vector<double> multiplySV( vector<double> v, double scalar) {
    vector<double> result;

    for (int i = 0; i < v.size(); i++) {
        double product = scalar * v[i];
        result.push_back(product);
    }

    return result;
}
// function to add two matrices
vector<vector<double>> addMatrices(vector<vector<double>> A, vector<vector<double>> B) {
    int n = A.size();
    int m = A[0].size();
    vector<vector<double>> C(n, vector<double>(m, 0.0));

    // check if matrices have same dimensions
    if (n != B.size() || m != B[0].size()) {
        cout << "Matrices have different dimensions, cannot add." << endl;
        return C;
    }

    // add matrices element-wise
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }

    return C;
}


vector<vector<double>> multiplyVV(vector<double>& v1, vector<double>& v2) {
    // Check if the size of the input vectors is valid for matrix multiplication
    if (v1.size() == 0 || v2.size() == 0 || v1.size() != v2.size()) {
        cout << "Error: Invalid vector size for matrix multiplication!" << endl;
        exit(1);
    }

    // Initialize the result matrix
    vector<vector<double>> result(v1.size(), vector<double>(v1.size(), 0));

    // Perform matrix multiplication
    for (int i = 0; i < v1.size(); i++) {
        for (int j = 0; j < v2.size(); j++) {
            result[i][j] = v1[i] * v2[j];
        }
    }

    return result;
}


double determinant(vector<vector<double>>& matrix) {
    int n = matrix.size();
    double det = 1.0;
    if (matrix.size() != matrix[0].size()) {
        throw invalid_argument("Matrix is not square!");
    }
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            double factor = matrix[j][i] / matrix[i][i];
            for (int k = i; k < n; k++) {
                matrix[j][k] -= factor * matrix[i][k];
            }
        }
        det *= matrix[i][i];
    }
    return det;
}

double determinantJ(vector<vector<double>> &mat) {
    if(mat.size() != 2 || mat[0].size() != 2 || mat[1].size() != 2) {
        cerr << "Error: matrix is not 2x2" << endl;
        return 0;
    }
    return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
}

double matrix_sum(const vector<vector<double>>& matrix) {
    double sum = 0.0;
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            sum += element;
        }
    }
    return sum;
}

vector<vector<double>> scalar_matrix_mult(double scalar, vector<vector<double>> matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<double>> result(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = scalar * matrix[i][j];
        }
    }

    return result;
}
//mat2 needs to be the bigger matrix
vector<vector<double>> combinematrix(vector<vector<double>> mat1, vector<vector<double>> mat2) {
    for (int i = 0; i < mat1.size(); i++) {
        mat2.push_back(mat1[i]);
    }
return mat2;
}


vector<double> scalarV(vector<double> vec, double scalar) {
    vector<double> result;
    for(int i = 0; i < vec.size(); i++) {
        result.push_back(scalar * vec[i]);
    }
    return result;
}
#endif // VECTORALGEBRA_H_INCLUDED
