#include <algorithm>
#include <iostream> // TEST
#include <fstream>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <cstdarg>
#include <cstdlib> 
#include <iomanip>

#define crash(...) exit(Crash(__VA_ARGS__))

using namespace std;

int Crash(const char *fmt, ...) {
    va_list ap;
    cerr << "\nEpic fail: \n";

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);

    cerr << "\n";
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, -1);
    return 0;
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

int setHalo(int& N, int new_I, int*& Part, int*& L2G, int*& G2L){
    if (Part[new_I] != 1){ // halo
        L2G[N] = new_I;
        G2L[new_I] = N;
        ++N;
    }
    return G2L[new_I];
}

void generate(int p_id, int Nx, int Ny, int K1, int K2, int Px, int Py, int& N, int& N0, int*& Part, int*& L2G, int*& G2L, int*& IA, int*& JA){
    int P = Px * Py;
    int K = K1 + K2;
    int trian_cnt = ((Ny * Nx) / K) * K2 + max(0, (Ny * Nx) % K - K1);
    int doubled_E = 0;
    int N_all = Ny * Nx + trian_cnt;

    //int* v_types = new int[N];
    int* cnt_neigh = new int[N_all];

    // MPI SPECIFIC - FILL PART
    N = 0;
    N0 = 0;
    Part = new int[N_all];
    int process_size = N_all / P;

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int old_I = i * Nx + j;
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            int cur_p = old_I / process_size;
            if (cur_p >= P){
                if (old_I % P == p_id){
                    Part[new_I] = 1;
                    ++N0;
                }
                else{
                    Part[new_I] = 0;
                }
            }
            else if (cur_p == p_id){
                Part[new_I] = 1;
                ++N0;
            }
            else{
                Part[new_I] = 0;
            }

            if (old_I % K >= K1){
                Part[new_I + 1] = Part[new_I];
            }
        }
    }

    N = N0;

    // #pragma omp parallel for
    // for (int i = 0; i < N_all; ++i){
    //     int cur_p = i / process_size;
    //     if (cur_p >= P){
    //         if (i % P == p_id){
    //             Part[i] = 1;
    //             ++N;
    //             ++N0;
    //         }
    //         else{
    //             Part[i] = 0;
    //         }
    //     }
    //     else if (cur_p == p_id){
    //         Part[i] = 1;
    //         ++N;
    //         ++N0;
    //     }
    //     else{
    //         Part[i] = 0;
    //     }
    // }

    // MPI SPECIFIC - FILL PART

    #pragma omp parallel
    {
        #pragma omp for reduction(+:doubled_E)
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
                    doubled_E += cnt_neigh[new_I] * Part[new_I];

                    cnt_neigh[new_I + 1] = 2;
                    if (i + 1 < Ny){ // нижний сосед
                        ++cnt_neigh[new_I + 1]; 
                    }
                    if (j + 1 < Nx){
                        ++cnt_neigh[new_I + 1]; // правый сосед
                    }
                    doubled_E += cnt_neigh[new_I + 1] * Part[new_I + 1];
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
                    doubled_E += cnt_neigh[new_I] * Part[new_I];
                }
            }
        }
    }

    IA = new int[N0];
    L2G = new int[N0];
    G2L = new int[N_all];
    JA = new int[doubled_E];

    #pragma omp parallel for
    for(int i = 0; i < N_all; ++i){
        G2L[i] = -1;
    }

    int prev_global, local_ind = 0;
    for(int i = 0; i < N_all; ++i){
        if (Part[i] == 1){
            if (local_ind == 0){
                IA[0] = 0;
            }
            else{
                IA[local_ind] = IA[local_ind - 1] + cnt_neigh[prev_global];
            }
            L2G[local_ind] = i;
            G2L[i] = local_ind;
            prev_global = i;
            ++local_ind;
        }
    }
    IA[local_ind] = doubled_E;


    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            int old_I = i * Nx + j;
            int cur_JA = IA[G2L[new_I]];
            if (old_I % K >= K1){
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1, Part, L2G, G2L);
                        //JA[cur_JA] = G2L[oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1];
                        ++cur_JA;
                    }else{
                        //JA[cur_JA] = G2L[oldInd2New(Nx, Ny, K1, K2, i - 1, j)];
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part, L2G, G2L);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = setHalo(N, new_I - 1, Part, L2G, G2L);
                    ++cur_JA;
                }
                JA[cur_JA] = setHalo(N, new_I, Part, L2G, G2L); // главная диагональ
                ++cur_JA;
                JA[cur_JA] = setHalo(N, new_I + 1, Part, L2G, G2L); // нижний треугольник = правый сосед
                ++cur_JA;
                
                JA[cur_JA] = setHalo(N, new_I, Part, L2G, G2L); // верхний треугольник = левый сосед
                ++cur_JA;
                JA[cur_JA] = setHalo(N, new_I + 1, Part, L2G, G2L); // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = setHalo(N, new_I + 2, Part, L2G, G2L);
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part, L2G, G2L);
                    ++cur_JA;
                }            
            }
            else{
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1, Part, L2G, G2L);
                        ++cur_JA;
                    }else{
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part, L2G, G2L);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = setHalo(N, new_I - 1, Part, L2G, G2L);
                    ++cur_JA;
                }
                JA[cur_JA] = setHalo(N, new_I, Part, L2G, G2L); // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = setHalo(N, new_I + 1, Part, L2G, G2L);
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part, L2G, G2L);
                    ++cur_JA;
                }
            }
        }
    }
}


int main(int argc, char** argv){
    int mpi_res;
    mpi_res = MPI_Init(&argc, &argv);
    if(mpi_res!= MPI_SUCCESS)
    {
        crash("MPI_Init failed (code %d)\n", mpi_res);
    }

    int NumProc;
    mpi_res = MPI_Comm_size(MPI_COMM_WORLD,&NumProc); // узнаем число процессов
    if(mpi_res!= MPI_SUCCESS)
    {
        crash("MPI_Comm_size failed (code %d)\n", mpi_res);
    }

    int MyID;
    mpi_res = MPI_Comm_rank(MPI_COMM_WORLD,&MyID); // узнаем номер данного процесса
    if(mpi_res!= MPI_SUCCESS)
    crash("MPI_Comm_rank failed (code %d)\n", mpi_res);

    double EPS = 1e-5;
    int MAXIT = 100;

    int Nx, Ny, K1, K2, Px, Py, N, N0, doubled_E, T, n;
    int *IA, *JA, *Part, *L2G, *G2L;
    double *A, *b, *x, t, res;
    
    ofstream logFile("output_log.txt", ios::out);
    if (!logFile.is_open()){
        return 1;
    }
    if (argc < 6){
        cout << "Запуск программы должен осуществляться в виде <program> Nx Ny K1 K2 T, где Nx, Ny - параметры расчётной сетки, K1, K2 - параметры разбиения на треугольники, T - число потоков"<<endl;
        cout << "На выходе создаётся файл output_log.txt с отладочной информацией" << endl;
        return 2;
    }
    Nx = atoi(argv[1]);
    Ny = atoi(argv[2]);
    K1 = atoi(argv[3]);
    K2 = atoi(argv[4]);
    Px = atoi(argv[5]);
    Py = atoi(argv[6]);
    T = atoi(argv[7]);
    if (Nx == 0 || Ny == 0 || K1 == 0 || K2 == 0 || T == 0){
        cout <<"Введены некорректные параметры. Запустите программу без параметров, чтобы увидеть help";
        return 3;
    }

    omp_set_num_threads(T);

    t = omp_get_wtime();
    generate(MyID, Nx, Ny, K1, K2, Px, Py, N, N0, Part, L2G, G2L, IA, JA);
    t = omp_get_wtime() - t;
    logFile << setprecision(5) << "Generate took: " << t << " seconds" << endl;

    logFile.close();

    delete [] IA;
    cerr << "good IA" << endl;
    delete [] JA;
    cerr << "good JA" << endl;
    delete [] Part;
    cerr << "good PART" << endl;
    delete [] L2G; 
    cerr << "good L2G" << endl;
    delete [] G2L;
    cerr << "good G2L" << endl;
    // delete [] A;
    // delete [] b;
    // delete [] x;
    return 0;
}