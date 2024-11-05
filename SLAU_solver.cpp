#include <algorithm>
#include <iostream> // TEST
#include <fstream>
#include <cmath>
#include <omp.h>
#include <iomanip>

using namespace std;

void SpMv(int N, int*& IA, int*& JA, double*& A, double*& x, double*& ans){
    #pragma omp parallel for
    for(int i = 0; i < N; ++i){
        double cur_sum = 0;
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            cur_sum += A[j] * x[JA[j]];
        }
        ans[i] = cur_sum;
    }
}

double dot(int N, double*& x, double*& y){
    double out = 0;
    #pragma omp parallel
    {
        double cur_res = 0;
        #pragma omp for
        for (int i = 0; i < N; ++i){
            cur_res += x[i] * y[i];
        }
        const int nt = omp_get_num_threads();
        const int tn = omp_get_thread_num();
        for(int i = 0; i < nt; ++i){
            #pragma omp barrier
            if (i == tn){
                out += cur_res;
            }
        }
    }
    return out;
}

void axpy(int N, double alpha, double*& x, double*& y, double*& ans) {
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        ans[i] = alpha * x[i] + y[i];  
    }
}

void vector_cp(int N, double*& dest, double*& src){
    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        dest[i] = src[i];
    }
}

void vector_fill(int N, double*& x, double alpha){
    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        x[i] = alpha;
    }
}

void get_M(int N, int*& IA, int*& JA, double*& A, double*& M){
    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            if (i == JA[j]){
                M[j] = 1.0 / A[j];
            }
        }
    }
}

double l2(int N, double*& x){
    double res = 0;
    for(int i = 0; i < N; ++i){
        res += x[i] * x[i];
    }
    return sqrt(res);
}

double solve(int N, int*& IA, int*& JA, double*& A, double*& b, double eps, int maxit, double*& x, int &n){
    int k = 0;
    double *x_k, *r_k, *M, *z_k, *p_k, *q_k, *x_k_prev, *p_k_prev, *r_k_prev, *tmp, ro_k, ro_k_prev, beta_k, alpha_k, t;
    x_k = new double[N];
    r_k = new double[N];
    M = new double[IA[N]];
    z_k = new double[N];
    p_k = new double[N];
    q_k = new double[N];
    x_k_prev = new double[N];
    p_k_prev = new double[N];
    r_k_prev = new double[N];
    //tmp = new double[N];

    vector_fill(N, x_k_prev, 0); // x_0
    //vector_fill(N, x_k, 0); // for prev
    vector_cp(N, r_k_prev, b); // r_0

    // TEST
    cout << "doubled_E: " << IA[N] << endl;
    t = omp_get_wtime();
    SpMv(N, IA, JA, M, r_k_prev, z_k);
    t = omp_get_wtime() - t;
    cout << "SpMv took: " << setprecision(3) << t << " seconds" << endl;

    /*

    do{
        ++k;
        vector_fill(IA[N], M, 0);
        get_M(N, IA, JA, A, M);
        SpMv(N, IA, JA, M, r_k_prev, z_k);
        ro_k = dot(N, r_k_prev, z_k);
        

        if(k == 1){
            vector_cp(N, p_k, z_k);
            vector_cp(N, p_k_prev, p_k);
        }
        else{
            beta_k = ro_k / ro_k_prev;
            //vector_cp(N, tmp, p_k);
            axpy(N, beta_k, p_k_prev, z_k, p_k);
            vector_cp(N, p_k_prev, p_k);
        }

        SpMv(N, IA, JA, A, p_k, q_k);
        alpha_k = ro_k / dot(N, p_k, q_k);
        //vector_cp(N, tmp, x_k);

        axpy(N, alpha_k, p_k, x_k_prev, x_k);
        vector_cp(N, x_k_prev, x_k);
        axpy(N, -alpha_k, q_k, r_k_prev, r_k);
        vector_cp(N, r_k_prev, r_k);

        ro_k_prev = ro_k;
    }
    while(ro_k > eps * eps && k < maxit);

    x = new double[N];
    vector_cp(N, x, x_k);
    n = k;

    return l2(N, r_k);
    */

    delete [] x_k;
    delete [] r_k;
    delete [] M;
    delete [] z_k;
    delete [] p_k;
    delete [] q_k;
    delete [] x_k_prev;
    delete [] p_k_prev;
    delete [] r_k_prev;
    //delete [] tmp;
}


int oldInd2New(int Nx, int Ny, int K1, int K2, int i, int j){
    int I, n, r, trian_cnt, new_I;
    int K = K1 + K2;
    I = i * Nx + j;
    n = I / K;
    r = I % K;
    trian_cnt = n * K2 + max(0, r - K1);
    new_I = I + trian_cnt;
    return new_I;
}

void generate(int Nx, int Ny, int K1, int K2, int& N, int*& IA, int*& JA){
    int K = K1 + K2;
    int trian_cnt = ((Ny * Nx) / K) * K2 + max(0, (Ny * Nx) % K - K1);
    int doubled_E = 0;
    N = Ny * Nx + trian_cnt;

    //int* v_types = new int[N];
    int* cnt_neigh = new int[N];

    #pragma omp parallel
    {
        int doubled_E_local = 0;
        #pragma omp for
        for(int i = 0; i < Ny; ++i){
            for(int j = 0; j < Nx; ++j){
                int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
                if ((i * Nx + j) % K >= K1){
                    //v_types[new_I] = 1; // верхний треугольник
                    //v_types[new_I + 1] = 2; // нижний треугольник

                    cnt_neigh[new_I] = 2;
                    if (i - 1 >= 0){ // верхний сосед
                        ++cnt_neigh[new_I];
                    }
                    if (j - 1 >= 0){ // левый сосед
                        ++cnt_neigh[new_I];
                    }
                    doubled_E_local += cnt_neigh[new_I];

                    cnt_neigh[new_I + 1] = 2;
                    if (i + 1 < Ny){ // нижний сосед
                        ++cnt_neigh[new_I + 1]; 
                    }
                    if (j + 1 < Nx){
                        ++cnt_neigh[new_I + 1]; // правый сосед
                    }
                    doubled_E_local += cnt_neigh[new_I + 1];
                }
                else{
                    //v_types[new_I] = 0; // обычная клетка
                    cnt_neigh[new_I] = 1;
                    if (i - 1 >= 0){ // верхний сосед
                        ++cnt_neigh[new_I];
                    }
                    if (i + 1 < Ny){ // нижний сосед
                        ++cnt_neigh[new_I]; 
                    }
                    if (j - 1 >= 0){ // левый сосед
                        ++cnt_neigh[new_I];
                    }
                    if (j + 1 < Nx){
                        ++cnt_neigh[new_I]; // правый сосед
                    }
                    doubled_E_local += cnt_neigh[new_I];
                }
            }
        }

        const int nt = omp_get_num_threads();
        const int tn = omp_get_thread_num();
        for(int i = 0; i < nt; ++i){
            #pragma omp barrier
            if (i == tn){
                doubled_E += doubled_E_local;
            }
        }
    }

    IA = new int[N + 1];
    JA = new int[doubled_E];

    for(int i = 0; i < N; ++i){
        if(i == 0){
            IA[0] = 0;
        }
        IA[i] = IA[i - 1] + cnt_neigh[i - 1];
    }
    IA[N] = doubled_E;

    #pragma omp parallel for
    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            if ((i * Nx + j) % K >= K1){
                int cur_JA = IA[new_I];
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1;
                        ++cur_JA;
                    }else{
                        JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = new_I - 1;
                    ++cur_JA;
                }
                JA[cur_JA] = new_I; // главная диагональ
                ++cur_JA;
                JA[cur_JA] = new_I + 1; // нижний треугольник = правый сосед
                ++cur_JA;
                
                JA[cur_JA] = new_I; // верхний треугольник = левый сосед
                ++cur_JA;
                JA[cur_JA] = new_I + 1; // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = new_I + 2;
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i + 1, j);
                    ++cur_JA;
                }            
            }
            else{
                int cur_JA = IA[new_I];
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1;
                        ++cur_JA;
                    }else{
                        JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = new_I - 1;
                    ++cur_JA;
                }
                JA[cur_JA] = new_I; // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = new_I + 1;
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = oldInd2New(Nx, Ny, K1, K2, i + 1, j);
                    ++cur_JA;
                }
            }
        }
    }
}

void fill(int N, int*& IA, int*& JA, double*& A, double*& b){
    double DIAG_COEFF = 1.234;

    int* diag = new int[N];
    A = new double[IA[N]];
    b = new double[N];

    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            if (i == JA[j]){
                diag[i] = j;
                A[j] = 0;
                continue;
            }
            A[j] = cos(i * JA[j] + i + JA[j]);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            if (i == JA[j]){
                continue;
            }
            A[diag[i]] += abs(A[j]);
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        A[diag[i]] *= DIAG_COEFF;
        b[i] = sin(i);
    }

    delete [] diag;
}



int main(int argc, char** argv){
    double EPS = 1e-5;
    int MAXIT = 100;

    int Nx, Ny, K1, K2, N, doubled_E, T, n;
    int *IA, *JA;
    double *A, *b, *x, t, res;
    ofstream logFile("output_log.txt", ios::out);
    if (!logFile.is_open()){
        return 1;
    }
    if (argc < 2){
        // TODO print help
    }
    ifstream inFile(argv[1], ios::in);
    inFile >> Nx >> Ny >> K1 >> K2 >> T;

    omp_set_num_threads(T);

    //logFile << "N: " << N <<endl;

    t = omp_get_wtime();
    generate(Nx, Ny, K1, K2, N, IA, JA);
    t = omp_get_wtime() - t;
    logFile << setprecision(3) << "Generate took: " << t << " seconds" << endl;
    // for (int i = 0; i < N+1; ++i){
    //     logFile << IA[i] << " ";
    // }
    logFile << endl;
    doubled_E = IA[N];
    // for (int i = 0; i < doubled_E; ++i){
    //     logFile << JA[i] << " ";
    // }  
    logFile << endl;
    t = omp_get_wtime();
    fill(N, IA, JA, A, b);
    t = omp_get_wtime() - t;
    logFile << setprecision(3) << "Fill took: " << t << " seconds" << endl;

    logFile << "N: " << N << endl;
    // for (int i = 0; i < doubled_E; ++i){
    //     logFile << A[i] << " ";
    // }  

    res = solve(N, IA, JA, A, b, EPS, MAXIT, x, n);

    // for(int i = 0; i < N; ++i){
    //     logFile << x[i] << " ";
    // }

    logFile.close();
    inFile.close();

    delete [] IA;
    delete [] JA;
    delete [] A;
    delete [] b;
    //delete [] x;
    return 0;
}