#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <fstream>
//#define DETERMINANT_PROCESS_SHOW
//#define INVERSE_PROCESS_SHOW
//#define LA_SOLVE_PROCESS_SHOW
//#define LA_JACOBI_METHOD_SHOW
//#define LA_SEIDEL_METHOD_SHOW
//#define LS_REGRESSION_PROCESS_SHOW
#define DOUBLE_PRECISION 4
#define GNUPLOT_NAME "gnuplot"

template<typename T>
T abs_(T value){
    return value >= 0 ? value : -value;
}

static bool random_init = false;

class DimensionalException : public std::exception{
public:
    DimensionalException() = default;
};

template<typename T>
class Matrix {
public:
    Matrix();
    Matrix(int rows, int columns, T value = 0);
    Matrix(const Matrix& matrix);
    ~Matrix();
    Matrix& transpose();
    static Matrix& randMatrix(int rows, int columns);

    int getRows() const;
    int getColumns() const;

    template<typename D>
    friend std::ostream& operator << (std::ostream& os, const Matrix<D>& matrix);

    template<typename D>
    friend std::istream& operator >> (std::istream& is, Matrix<D>& matrix);

    template<typename D>
    friend bool operator==(const Matrix<D>& l, const Matrix<D>& r);

    template<typename D>
    friend bool operator!=(const Matrix<D>& l, const Matrix<D>& r);

    Matrix& operator=(const Matrix& r);
    Matrix& operator+=(const Matrix& r);
    Matrix& operator-=(const Matrix& r);
    Matrix& operator*=(const Matrix& r);

    template<typename D>
    friend Matrix<D>& operator-(const Matrix<D>& r);

    template<typename D>
    friend Matrix<D>& operator-(const Matrix<D>& l, const Matrix<D>& r);

    template<typename D>
    friend Matrix<D>& operator+(const Matrix<D>& l, const Matrix<D>& r);

    template<typename D>
    friend Matrix<D>& operator*(const Matrix<D>& l, const Matrix<D>& r);

    template<typename D>
    operator Matrix<D>();

    template<typename D>
    Matrix& multiplyBy(D value) const;

    template<typename D>
    Matrix& divideBy(D value) const;

    T* operator[](int row) const;
protected:
    T* arr = nullptr;
    int rows;
    int columns;
};

template<typename T>
Matrix<T>::Matrix(): rows(3), columns(3) {
    arr = new T[rows * columns];
    std::fill(arr, arr + rows * columns, 0);
}

template<typename T>
Matrix<T>::Matrix(int rows, int columns, T value): rows(rows), columns(columns) {
    arr = new T[rows * columns];
    std::fill(arr, arr + rows * columns, value);
}

template<typename T>
Matrix<T>::Matrix(const Matrix &matrix): rows(matrix.rows), columns(matrix.columns)  {
    arr = new T[rows * columns];
    std::copy(matrix.arr, matrix.arr + rows * columns, arr);
}

template<typename T>
Matrix<T>::~Matrix() {
    if (arr != nullptr){
        delete[] arr;
        arr = nullptr;
        rows = 0;
        columns = 0;
    }
}

template<typename T>
Matrix<T>& Matrix<T>::randMatrix(int rows, int columns) {
    Matrix<T> *out = new Matrix<T>(rows, columns);

    if (!random_init){
        srand(time(nullptr));
        random_init = true;
    }

    for (int i = 0; i < rows * columns; i++) {
        out->arr[i] = (T)(rand() % 100) / 10;
    }

    return *out;
}

template<typename D>
std::ostream &operator<<(std::ostream &os, const Matrix<D> &matrix) {
    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.columns; j++) {
            os << matrix.arr[i * matrix.columns + j];

            if (j < matrix.columns - 1){
                os << " ";
            }
        }
        os << "\n";
    }

    return os;
}

template<>
std::ostream &operator<<(std::ostream &os, const Matrix<double> &matrix) {
    std::cout << std::fixed;

    for (int i = 0; i < matrix.rows; i++) {
        for (int j = 0; j < matrix.columns; j++) {
            if (abs_(matrix.arr[i * matrix.columns + j]) >= 1/std::pow(10, DOUBLE_PRECISION)){
                os << std::setprecision(DOUBLE_PRECISION) << matrix.arr[i * matrix.columns + j];
            } else{
                os << std::setprecision(DOUBLE_PRECISION) << 0.00;
            }

            if (j < matrix.columns - 1){
                os << " ";
            }
        }
        os << "\n";
    }

    return os;
}

template<typename D>
bool operator==(const Matrix<D> &l, const Matrix<D> &r) {
    if (l.columns != r.columns || l.rows != r.rows)
        return false;

    for (int i = 0; i < l.columns * l.rows; i++) {
        if (l.arr[i] != r.arr[i]){
            return false;
        }
    }

    return true;
}

template<typename D>
bool operator!=(const Matrix<D> &l, const Matrix<D> &r) {
    return !(l == r);
}

template<typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix &r) {
    if (this->arr != nullptr){
        delete[] arr;
        rows = 0;
        columns = 0;
    }

    rows = r.rows;
    columns = r.columns;
    arr = new T[rows * columns];
    std::copy(r.arr, r.arr + rows * columns, arr);

    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix &r) {
    if (columns != r.columns || rows != r.rows){
        throw DimensionalException();
    }

    for (int i = 0; i < rows * columns; i++) {
        arr[i] += r.arr[i];
    }

    return *this;
}

template<typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix &r) {
    if (columns != r.columns || rows != r.rows){
        throw DimensionalException();
    }

    for (int i = 0; i < rows * columns; i++) {
        arr[i] -= r.arr[i];
    }

    return *this;
}

template<typename D>
Matrix<D> &operator*(const Matrix<D> &l, const Matrix<D> &r) {
    if (l.columns != r.rows){
        throw DimensionalException();
    }

    Matrix<D> *out = new Matrix<D>(l.rows, r.columns);

    for (int i = 0; i < l.rows; i++) {
        D* l_row = l.arr + i * l.columns;
        D* out_row = out->arr + i * out->columns;

        for (int j = 0; j < r.rows; j++) {
            D* r_row = r.arr + j * r.columns;
            D l_element = l_row[j];

            for (int k = 0; k < r.columns; k++) {
                out_row[k] += l_element * r_row[k];
            }
        }
    }

    return *out;
}

template<typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix &r) {
    Matrix<T> *result = &((*this) * r);
    *this = *result;

    return *this;
}

template<typename T>
template<typename D>
Matrix<T> &Matrix<T>::multiplyBy(D value) const{
    Matrix<T> *out = new Matrix<T>(rows, columns);

    for (int i = 0; i < rows * columns; i++) {
        out->arr[i] = arr[i] * (T)value;
    }

    return *out;
}

template<typename D>
Matrix<D> &operator-(const Matrix<D> &r) {
    return r.multiplyBy(-1);
}

template<typename D>
Matrix<D> &operator-(const Matrix<D> &l, const Matrix<D> &r) {
    if (l.columns != r.columns || l.rows != r.rows){
        throw DimensionalException();
    }

    Matrix<D> *out = new Matrix<D>(l.rows, l.columns);

    for (int i = 0; i < l.rows * l.columns; i++) {
        out->arr[i] = l.arr[i] - r.arr[i];
    }

    return *out;
}

template<typename D>
Matrix<D> &operator+(const Matrix<D> &l, const Matrix<D> &r) {
    if (l.columns != r.columns || l.rows != r.rows){
        throw DimensionalException();
    }

    Matrix<D> *out = new Matrix<D>(l.rows, l.columns);

    for (int i = 0; i < l.rows * l.columns; i++) {
        out->arr[i] = l.arr[i] + r.arr[i];
    }

    return *out;
}

template<typename T>
template<typename D>
Matrix<T>::operator Matrix<D>() {
    Matrix<D> *out = new Matrix<D>(rows, columns);

    for (int i = 0; i < rows * columns; i++) {
        out->arr[i] = (D)arr[i];
    }

    return *out;
}

template<typename T>
template<typename D>
Matrix<T> &Matrix<T>::divideBy(D value) const {
    if (value == 0){
        throw std::exception();
    }

    Matrix<T> *out = new Matrix<T>(rows, columns);

    for (int i = 0; i < rows * columns; i++) {
        out->arr[i] = arr[i] / (T)value;
    }

    return *out;
}

template<typename T>
Matrix<T> &Matrix<T>::transpose() {
    Matrix<T> *out = new Matrix<T>(columns, rows);

    for (int i = 0; i < columns; i++) {
        for (int j = 0; j < rows; j++) {
            out->arr[i * rows + j] = arr[j * columns + i];
        }
    }

    return *out;
}

template<typename D>
std::istream &operator>>(std::istream &is, Matrix<D> &matrix) {
    for (int i = 0; i < matrix.rows * matrix.columns; i++) {
        is >> matrix.arr[i];
    }

    return is;
}

template<typename T>
int Matrix<T>::getRows() const {
    return rows;
}

template<typename T>
int Matrix<T>::getColumns() const {
    return columns;
}

template<typename T>
T* Matrix<T>::operator[](int row) const{
    return arr + row * columns;
}

template <typename T>
class SquareMatrix : public Matrix<T>{
public:
    SquareMatrix();
    SquareMatrix(int size);
    SquareMatrix(const SquareMatrix& squareMatrix);

    static SquareMatrix& randMatrix(int size);
    T determinant() const;
    SquareMatrix& inverse() const;
};

template<typename T>
SquareMatrix<T>::SquareMatrix() : Matrix<T>() {}

template<typename T>
SquareMatrix<T>::SquareMatrix(int size) : Matrix<T>(size, size) {}

template<typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix &squareMatrix) : Matrix<T>(squareMatrix){}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::randMatrix(int size) {
    return *(SquareMatrix<T>*)(&Matrix<T>::randMatrix(size, size));
}

template<typename T>
class IdentityMatrix : public SquareMatrix<T>{
public:
    IdentityMatrix() = delete;
    IdentityMatrix(int size) = delete;
    IdentityMatrix(const IdentityMatrix& identityMatrix) = delete;
    IdentityMatrix(IdentityMatrix&& identityMatrix) = delete;
    void * operator new(size_t size) = delete;

    template<typename D>
    friend const IdentityMatrix<D>& I(int size);
};

template<typename D>
const IdentityMatrix<D>& I(int size) {
    SquareMatrix<D> *out = new SquareMatrix<D>(size);

    for (int i = 0; i < size; i++) {
        (*out)[i][i] = 1;
    }

    return *(const IdentityMatrix<D>*)out;
}

template<typename T>
class EliminationMatrix : public SquareMatrix<T>{
public:
    EliminationMatrix() = delete;
    EliminationMatrix(const Matrix<T>& matrix, int row, int column);
    EliminationMatrix(const EliminationMatrix& eliminationMatrix);
};

template<typename T>
EliminationMatrix<T>::EliminationMatrix(const Matrix<T> &matrix, int row, int column) {
    this->rows = matrix.getRows();
    this->columns = matrix.getRows();

    this->arr = new T[this->rows * this->columns];
    std::fill(this->arr, this->arr + this->rows * this->columns, 0);

    for (int i = 0; i < this->rows * this->columns; i += (this->rows + 1)) {
        this->arr[i] = 1;
    }

    this->arr[row * this->columns + column] = -matrix[row][column] / matrix[column][column];
}

template<typename T>
EliminationMatrix<T>::EliminationMatrix(const EliminationMatrix &eliminationMatrix) : SquareMatrix<T>(eliminationMatrix) {}

template<typename T>
class PermutationMatrix : public SquareMatrix<T>{
public:
    PermutationMatrix();
    PermutationMatrix(int size, int row1, int row2);
    PermutationMatrix(const PermutationMatrix& permutationMatrix);
};

template<typename T>
PermutationMatrix<T>::PermutationMatrix() : SquareMatrix<T>(I<int>(3)){}

template<typename T>
PermutationMatrix<T>::PermutationMatrix(int size, int row1, int row2) : SquareMatrix<T>(size){
    for (int i = 0; i < size; i++) {
        if (i != row1 && i != row2){
            this->arr[i*(size + 1)] = 1;
        } else if ( i == row1){
            this->arr[i*size + row2] = 1;
        } else{
            this->arr[i*size + row1] = 1;
        }
    }
}

template<typename T>
PermutationMatrix<T>::PermutationMatrix(const PermutationMatrix &permutationMatrix) : SquareMatrix<T>(permutationMatrix){}

template<typename T>
T SquareMatrix<T>::determinant() const{
    SquareMatrix<T> copy(*this);
    int factor = 1;
    int step = 1;

    for (int i = 0; i < this->rows - 1; i++) {
        int swap_with = i;
        T max_value = abs_(this->arr[i * this->rows + i]);

        bool all_rows_zero = (max_value == 0);

        for (int j = i + 1; j < this->rows; j++) {
            if (abs_(this->arr[j * this->rows + i]) > max_value){
                swap_with = j;
                max_value = abs_(this->arr[j * this->rows + i]);
                break;
            }
        }

        if (!all_rows_zero){
            if (swap_with != i){
                PermutationMatrix<T> p(this->rows, swap_with, i);
                copy = *(SquareMatrix<T>*)&(p * copy);
                factor *= -1;
#ifdef DETERMINANT_PROCESS_SHOW
                printf("step #%d: permutation\n", step);
                std::cout << copy;
#endif
                step++;
            }

            for (int j = i + 1; j < this->rows; j++) {
                EliminationMatrix<T> e(copy, j, i);
                copy = *(SquareMatrix<T>*)&(e * copy);
#ifdef DETERMINANT_PROCESS_SHOW
                printf("step #%d: elimination\n", step);
                std::cout << copy;
#endif
                step++;
            }
        }
    }

    T result = 1;
    for (int i = 0; i < this->rows; i++) {
        result *= copy[i][i];
    }

    return result * (T)factor;
}

template<typename T>
SquareMatrix<T> &SquareMatrix<T>::inverse() const {
    SquareMatrix<T> *out = new SquareMatrix<T>(*this);
    Matrix<T> augMatrix(out->rows, out->rows * 2);

    for (int i = 0; i < out->rows; i++) {
        for (int j = 0; j < out->rows; j++) {
            augMatrix[i][j] = this->arr[i * out->rows + j];
        }
    }

    for (int i = 0; i < out->rows; i++) {
        augMatrix[i][i + out->rows] = 1;
    }

    int step = 0;
#ifdef INVERSE_PROCESS_SHOW
    printf("step #%d: Augmented Matrix\n", step);
    std::cout << augMatrix;
    printf("Direct way:\n");
#endif
    step++;

    for (int i = 0; i < this->rows; i++) {
        int swap_with = i;
        T max_value = abs_(augMatrix[i][i]);

        bool all_rows_zero = (max_value == 0);

        for (int j = i + 1; j < this->rows; j++) {
            if (abs_(augMatrix[j][i]) > max_value){
                swap_with = j;
                max_value = abs_(augMatrix[j][i]);
                all_rows_zero = false;
            }
        }

        if (!all_rows_zero && i != this->rows - 1){
            if (swap_with != i){
                PermutationMatrix<T> p(this->rows, swap_with, i);
                augMatrix = p * augMatrix;
#ifdef INVERSE_PROCESS_SHOW
                printf("step #%d: permutation\n", step);
                std::cout << augMatrix;
#endif
                step++;
            }

            for (int j = i + 1; j < this->rows; j++) {
                if (augMatrix[j][i] != 0){
                    EliminationMatrix<T> e(augMatrix, j, i);
                    augMatrix = e * augMatrix;
#ifdef INVERSE_PROCESS_SHOW
                    printf("step #%d: elimination\n", step);
                    std::cout << augMatrix;
#endif
                    step++;
                }
            }
        } else if (all_rows_zero){
            return *out;
        }
    }

#ifdef INVERSE_PROCESS_SHOW
    printf("Way back:\n");
#endif

    for (int i = this->rows - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (augMatrix[j][i] != 0){
                EliminationMatrix<T> e(augMatrix, j, i);
                augMatrix = e * augMatrix;
#ifdef INVERSE_PROCESS_SHOW
                printf("step #%d: elimination\n", step);
                std::cout << augMatrix;
#endif
                step++;
            }
        }
    }

    for (int i = 0; i < this->rows; i++) {
        T value = augMatrix[i][i];
        for (int j = 0; j < 2 * this->rows; j++) {
            augMatrix[i][j] /= value;
        }
    }

#ifdef INVERSE_PROCESS_SHOW
    printf("Diagonal normalization:\n");
    std::cout << augMatrix;
#endif

    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->rows; j++) {
            out->arr[i * this->rows + j] = augMatrix[i][j + this->rows];
        }
    }

    return *out;
}

template<typename T>
class ColumnVector : public Matrix<T>{
public:
    ColumnVector();
    ColumnVector(int n);
    ColumnVector(const ColumnVector& columnVector);
};

template<typename T>
ColumnVector<T>::ColumnVector() : Matrix<T>(3, 1){}

template<typename T>
ColumnVector<T>::ColumnVector(int n) : Matrix<T>(n, 1){}

template<typename T>
ColumnVector<T>::ColumnVector(const ColumnVector &columnVector) : Matrix<T>(columnVector){}

template<typename D>
ColumnVector<D> &solveSimpleLA(const SquareMatrix<D> &A, const ColumnVector<D> &b) {
    SquareMatrix<D> A_copy(A);
    ColumnVector<D> *out = new ColumnVector<D>(b);
    Matrix<D> copy(A_copy);

    int step = 0;
#ifdef LA_SOLVE_PROCESS_SHOW
    printf("step #%d:\n", step);
    std::cout << A_copy;
    std::cout << *out;
#endif
    step++;

    for (int i = 0; i < A_copy.getRows(); i++) {
        int swap_with = i;
        D max_value = abs_(A_copy[i][i]);

        bool all_rows_zero = (max_value == 0);

        for (int j = i + 1; j < A_copy.getRows(); j++) {
            if (abs_(A_copy[j][i]) > max_value){
                swap_with = j;
                max_value = abs_(A_copy[j][i]);
                all_rows_zero = false;
            }
        }

        if (!all_rows_zero && i != A_copy.getRows() - 1){
            if (swap_with != i){
                PermutationMatrix<D> p(A_copy.getRows(), swap_with, i);
                A_copy = *(SquareMatrix<D>*)&(p * A_copy);
                out = (ColumnVector<D>*)&(p * (*out));
#ifdef LA_SOLVE_PROCESS_SHOW
                printf("step #%d: permutation\n", step);
                std::cout << A_copy;
                std::cout << *out;
#endif
                step++;
            }

            for (int j = i + 1; j < A_copy.getRows(); j++) {
                if (A_copy[j][i] != 0){
                    EliminationMatrix<D> e(A_copy, j, i);
                    A_copy = *(SquareMatrix<D>*)&(e * A_copy);
                    out = (ColumnVector<D>*)&(e * (*out));
#ifdef LA_SOLVE_PROCESS_SHOW
                    printf("step #%d: elimination\n", step);
                    std::cout << A_copy;
                    std::cout << *out;
#endif
                    step++;
                }
            }
        } else if (all_rows_zero){
            return (ColumnVector<D>&)b;
        }
    }

#ifdef LA_SOLVE_PROCESS_SHOW
    printf("Way back:\n");
#endif

    for (int i = A_copy.getRows() - 1; i > 0; i--) {
        for (int j = i - 1; j >= 0; j--) {
            if (A_copy[j][i] != 0){
                EliminationMatrix<D> e(A_copy, j, i);
                A_copy = *(SquareMatrix<D>*)&(e * A_copy);
                out = (ColumnVector<D>*)&(e * (*out));
#ifdef LA_SOLVE_PROCESS_SHOW
                printf("step #%d: elimination\n", step);
                std::cout << A_copy;
                std::cout << *out;
#endif
                step++;
            }
        }
    }

    for (int i = 0; i < A_copy.getRows(); i++) {
        D value = A_copy[i][i];
        for (int j = 0; j < A_copy.getRows(); j++) {
            A_copy[i][j] /= value;
        }
        (*out)[i][0] /= value;
    }

#ifdef LA_SOLVE_PROCESS_SHOW
    printf("Diagonal normalization:\n");
    std::cout << A_copy;
    std::cout << *out;
#endif

    return *out;
}

template<typename T>
ColumnVector<T> &solveSimpleLAJacobiMethod(const SquareMatrix<T>& A, const ColumnVector<T>& b, double e){
    ColumnVector<T> *out = new ColumnVector<T>(b);
    SquareMatrix<T> *copy = new SquareMatrix<T>(A);

    for (int i = 0; i < A.getRows(); i++) {
        T sum = 0;
        for (int j = 0; j < A.getRows(); j++) {
            if (j != i){
                sum += abs_(A[i][j]);
            }
        }

        if (sum >= A[i][i]){
            delete copy;
#ifdef LA_JACOBI_METHOD_SHOW
            printf("The method is not applicable!\n");
#endif
            return (ColumnVector<T>&)b;
        }
    }

    for (int i = 0; i < A.getRows(); i++) {
        for (int j = 0; j < A.getRows(); j++) {
            if (i != j){
                (*copy)[i][j] /= -A[i][i];
            } else{
                (*copy)[i][i] = 0;
            }
        }
        (*out)[i][0] /= A[i][i];
    }

#ifdef LA_JACOBI_METHOD_SHOW
    printf("alpha:\n");
    std::cout << *copy;
    printf("beta:\n");
    std::cout << *out;
#endif

    int step = 0;
#ifdef LA_JACOBI_METHOD_SHOW
    printf("x(%d):\n", step);
    std::cout << *out;
#endif
    step++;

    ColumnVector<T> b_copy(*out);
    double current_e = e;
    while (current_e >= e){
        ColumnVector<T> new_out = *(ColumnVector<T>*)&(*copy * (*out) + b_copy);

        current_e = std::sqrt(((new_out - (*out)).transpose() * (new_out - (*out)))[0][0]);
        (*out) = new_out;

#ifdef LA_JACOBI_METHOD_SHOW
        printf("e: %.4lf\n", current_e);
        printf("x(%d):\n", step);
        std::cout << *out;
#endif

        step++;
    }

    delete copy;
    return *out;
}

template<typename T>
ColumnVector<T> &solveSimpleLASeidelMethod(const SquareMatrix<T>& A, const ColumnVector<T>& b, double e){
    ColumnVector<T> beta(b);
    SquareMatrix<T> alpha(A);
    SquareMatrix<T> B(A.getRows());
    SquareMatrix<T> C(A.getRows());

    for (int i = 0; i < A.getRows(); i++) {
        T sum = 0;
        for (int j = 0; j < A.getRows(); j++) {
            if (j != i){
                sum += abs_(A[i][j]);
            }
        }

        if (sum >= A[i][i]){
#ifdef LA_SEIDEL_METHOD_SHOW
            printf("The method is not applicable!\n");
#endif
            return (ColumnVector<T>&)b;
        }
    }

    for (int i = 0; i < A.getRows(); i++) {
        for (int j = 0; j < A.getRows(); j++) {
            if (i != j){
                alpha[i][j] /= -A[i][i];

                if (i > j){
                    B[i][j] = alpha[i][j];
                } else{
                    C[i][j] = alpha[i][j];
                }
            } else{
                alpha[i][j] = 0;
            }
        }

        beta[i][0] /= A[i][i];
    }

    ColumnVector<T> *x = new ColumnVector<T>(beta);

#ifdef LA_SEIDEL_METHOD_SHOW
    printf("beta:\n");
    std::cout << *x;
    printf("alpha:\n");
    std::cout << alpha;
    printf("B:\n");
    std::cout << B;
    printf("C:\n");
    std::cout << C;
    printf("I-B:\n");
    std::cout << I<T>(B.getRows()) - B;
    printf("(I-B)_-1:\n");
    std::cout << ((SquareMatrix<T>*)&(I<T>(B.getRows()) - B))->inverse();
#endif

    int step = 0;
#ifdef LA_SEIDEL_METHOD_SHOW
    printf("x(%d):\n", step);
    std::cout << *x;
#endif
    step++;

    SquareMatrix<T> L(A.getRows());
    SquareMatrix<T> U(A.getRows());

    for (int i = 0; i < A.getRows(); i++) {
        for (int j = 0; j < A.getRows(); j++) {
            if (i >= j){
                L[i][j] = A[i][j];
            } else{
                U[i][j] = A[i][j];
            }
        }
    }

    SquareMatrix<T> L_inverse = L.inverse();

    double current_e = e;
    while (current_e >= e){
        ColumnVector<T> x_new((ColumnVector<T>&)(L_inverse * (b - U * (*x))));

        current_e = std::sqrt(((x_new - (*x)).transpose() * (x_new - (*x)))[0][0]);
        (*x) = x_new;

#ifdef LA_SEIDEL_METHOD_SHOW
        printf("e: %.4lf\n", current_e);
        printf("x(%d):\n", step);
        std::cout << *x;
#endif
        step++;
    }

    return *x;
}

template<typename T>
ColumnVector<T>& lsRegression(const Matrix<T>& data, int n){
    Matrix<T> A(data.getRows(), n + 1);
    ColumnVector<T> b(data.getRows());

    for (int i = 0; i < A.getRows(); i++) {
        b[i][0] = data[i][1];
        A[i][0] = 1;
        for (int j = 1; j < A.getColumns(); j++) {
            A[i][j] = (T)std::pow(data[i][0], j);
        }
    }

#ifdef LS_REGRESSION_PROCESS_SHOW
    printf("A:\n");
    std::cout << A;
#endif

    b = (ColumnVector<T>&)(A.transpose() * b);
    A = A.transpose() * A;

#ifdef LS_REGRESSION_PROCESS_SHOW
    printf("A_T*A:\n");
    std::cout << A;
#endif

    A = ((SquareMatrix<T>&)A).inverse();

#ifdef LS_REGRESSION_PROCESS_SHOW
    printf("(A_T*A)^-1:\n");
    std::cout << A;
    printf("A_T*b:\n");
    std::cout << b;
#endif

    return (ColumnVector<T>&)(A * b);
}

template<typename T>
Matrix<T>& getMatrix(){
    int rows;
    int columns;

    std::cin >> rows;
    std::cin >> columns;

    Matrix<T> *out = new Matrix<T>(rows, columns);

    std::cin >> *out;

    return *out;
}

template<typename T>
SquareMatrix<T>& getSquareMatrix(){
    int size;

    std::cin >> size;

    SquareMatrix<T> *out = new SquareMatrix<T>(size);

    std::cin >> *out;

    return *out;
}

template<typename T>
ColumnVector<T>& getColumnVector(){
    int size;

    std::cin >> size;

    ColumnVector<T> *out = new ColumnVector<T>(size);

    std::cin >> *out;

    return *out;
}

Matrix<double>& generateData(int n, int x1, int y1, int x2, int y2, double dispersion){
    Matrix<double> *out = new Matrix<double>(n, 2);

    double k = (double)(y2 - y1) / (x2 - x1);
    double b = y1 - k * x1;
    double sixthX = (double)(x2 - x1) / 6;
    double sixthY = (double)(y2 - y1) / 6;

    int factor = 2 * (rand() % 2) - 1;

    double xs1, ys1, xs2, ys2, xs3, ys3;

    xs1 = ((double)(rand() % 100) / 100) * sixthX + x1;
    xs2 = ((double)(rand() % 100) / 100) * sixthX + x1 + 2 * sixthX;
    xs3 = ((double)(rand() % 100) / 100) * sixthX + x1 + 4 * sixthX;

    ys1 = k * xs1 + b + factor * sixthX * (1 + ((double)(rand() % 100) / 100));
    ys2 = k * xs2 + b - factor * sixthX * (1 + ((double)(rand() % 100) / 100));
    ys3 = k * xs3 + b + factor * sixthX * (1 + ((double)(rand() % 100) / 100));

    Matrix<double> data(3, 2);
    data[0][0] = xs1; data[0][1] = ys1;
    data[1][0] = xs2; data[1][1] = ys2;
    data[2][0] = xs3; data[2][1] = ys3;

    ColumnVector<double> parabola = lsRegression(data, 2);

    double minX = x1, maxX = x2;
    double A = parabola[2][0];
    double B = parabola[1][0];
    double C = parabola[0][0];

    if (factor > 0){
        if (B * B - 4 * A * (C - y2) >= 0){
            double x1_ = (-B + std::sqrt(B * B - 4 * A * (C - y2))) / (2 * A);
            double x2_ = (-B - std::sqrt(B * B - 4 * A * (C - y2))) / (2 * A);
            maxX = maxX > x1_ ? x1_ : maxX;
            minX = minX > x2_ ? minX : x2_;
        }
    } else{
        if (B * B - 4 * A * (C - y1) >= 0){
            double x1_ = (-B + std::sqrt(B * B - 4 * A * (C - y1))) / (2 * A);
            double x2_ = (-B - std::sqrt(B * B - 4 * A * (C - y1))) / (2 * A);
            maxX = maxX > x1_ ? x1_ : maxX;
            minX = minX > x2_ ? minX : x2_;
        }
    }

    double areaX = (maxX - minX) / n;

    for (int i = 0; i < n; i++) {
        double x = ((double)(rand() % 100) / 100) * areaX + i * areaX + minX;
        double y = parabola[0][0] + parabola[0][1] * x + parabola[0][2] * x * x;
        y = dispersion - ((double)(rand() % 100) / 100) * dispersion * 2 + y;

        (*out)[i][0] = x;
        (*out)[i][1] = y;
    }

    return *out;
}

int main() {
#ifdef WIN32
    FILE *pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE *pipe = popen(GNUPLOT_NAME, "w");
#endif
    srand(time(nullptr));

    if (pipe != nullptr)
    {
        const int MAX_X = 10;
        const int MIN_X = 0;
        const int MAX_Y = 10;
        const int MIN_Y = 0;
        const int MIN_NUMBER_OF_POINTS = 15;
        const int MAX_NUMBER_OF_POINTS = 20;
        const double dispersion = 0.8;

        std::string input;
        std::cout << "(Y) - Generate data\n";
        std::cout << "(N) - Load data.txt\n";
        std::cin >> input;



        int n = rand() % (MAX_NUMBER_OF_POINTS - MIN_NUMBER_OF_POINTS + 1) + MIN_NUMBER_OF_POINTS;
        Matrix<double> data;

        if (input == "Y"){
            std::ofstream outfile ("data.txt");
            data = generateData(n, MIN_X, MIN_Y, MAX_X, MAX_Y, dispersion);

            outfile << data.getRows() << " " << data.getColumns() << "\n";
            outfile << data << "\n";

            outfile.close();
        } else if (input == "N"){
            std::ifstream input_file("data.txt");

            int rows, columns;
            input_file >> rows;
            input_file >> columns;

            data = Matrix<double>(rows, columns);
            input_file >> data;
        } else{
            std::cout << "Fatal error!\n";
            exit(404);
        }

        ColumnVector<double> polynomial2 = lsRegression(data, 4);
        ColumnVector<double> polynomial1 = lsRegression(data, 2);

        fprintf(pipe, "set xrange [%d: %d]\n", MIN_X, MAX_X);
        fprintf(pipe, "set yrange [%d: %d]\n", MIN_Y, MAX_Y);
        fprintf(pipe, "set grid\n");
        fprintf(pipe, "set xlabel 'X'\n");
        fprintf(pipe, "set ylabel 'Y'\n");
        fprintf(pipe, "set xtics %d,1,%d\n", MIN_X, MAX_X);
        fprintf(pipe, "set ytics %d,1,%d\n", MIN_Y, MAX_Y);
        fprintf(pipe, "plot '-' with points pointtype 5 title 'Data', ");
        fprintf(pipe, "%lf*x*x + %lf*x + %lf title 'Polynomial with degree 2', ", polynomial1[2][0], polynomial1[1][0], polynomial1[0][0]);
        fprintf(pipe, "%lf*x*x*x*x + %lf*x*x*x + %lf*x*x + %lf*x + %lf title 'Polynomial with degree 4'\n", polynomial2[4][0], polynomial2[3][0], polynomial2[2][0], polynomial2[1][0], polynomial2[0][0]);

        for (int i = 0; i < n; i++) {
            fprintf(pipe, "%lf %lf\n", data[i][0], data[i][1]);
        }
        fprintf(pipe, "%s\n", "e");
        fflush(pipe);

        printf("Press any key...\n");
        int useless;
        std::cin >> useless;
#ifdef WIN32
        _pclose(pipe);
#else
        pclose(pipe);
#endif
    }
    else{
        std::cout << "Could not open pipe" << std::endl;
    }

    return 0;
}