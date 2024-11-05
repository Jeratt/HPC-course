#include <algorithm>
#include <iostream> // TEST
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;

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
}



int main(int argc, char** argv){
    int Nx, Ny, K1, K2, N, doubled_E, T;
    int *IA, *JA;
    double *A, *b;
    ofstream logFile("error_log.txt", ios::out);
    if (!logFile.is_open()){
        return 1;
    }
    if (argc < 2){
        // TODO print help
    }
    ifstream inFile(argv[1], ios::in);
    inFile >> Nx >> Ny >> K1 >> K2 >> T;

    omp_set_num_threads(T);

    logFile << "N: " << N <<endl;
    generate(Nx, Ny, K1, K2, N, IA, JA);
    for (int i = 0; i < N+1; ++i){
        logFile << IA[i] << " ";
    }
    logFile << endl;
    doubled_E = IA[N];
    for (int i = 0; i < doubled_E; ++i){
        logFile << JA[i] << " ";
    }  
    logFile << endl;
    fill(N, IA, JA, A, b);
    for (int i = 0; i < doubled_E; ++i){
        logFile << A[i] << " ";
    }  

    logFile.close();
    inFile.close();

    delete [] IA;
    delete [] JA;
    return 0;
}