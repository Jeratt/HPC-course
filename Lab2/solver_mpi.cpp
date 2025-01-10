#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include <cstdarg>
#include <cstdlib> 
#include <iomanip>
#include <vector>

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

class tCommScheme {
private:
    vector<int> Send;       // List of cells to send to all neighbors
    vector<int> Recv;       // List of cells to receive from all neighbors
    vector<int> SendOffset;     // Offsets for send lists for each neighbor
    vector<int> RecvOffset;     // Offsets for receive lists for each neighbor
    vector<int> Neighbours;     // List of neighbor process IDs
    MPI_Comm MyComm;                 // MPI communicator for the group

public:
    // Constructor
    tCommScheme(MPI_Comm comm = MPI_COMM_WORLD) : MyComm(comm) {}

    int* GetSendList() {
        return Send.data();
    }

    int* GetRecvList() {
        return Recv.data();
    }

    int* GetSendOffset() {
        return SendOffset.data();
    }

    int* GetRecvOffset() {
        return RecvOffset.data();
    }

    int* GetListOfNeigbours() {
        return Neighbours.data();
    }

    int GetNumOfNeighbours() {
        return Neighbours.size();
    }

    MPI_Comm GetMyComm() const {
        return MyComm;
    }

    void Initialize(int N,
                    int N0,
                    int P,
                    int*& IA,
                    int *& JA,
                    int *& Part,
                    int *& L2G,
                    int *& G2L)
    {
        SendOffset.clear();
        RecvOffset.clear();
        Neighbours.clear();

        vector<vector<int>> SendToProcess;
        vector<vector<int>> RecvFromProcess;
        for(int i=0;i<P;++i){
            SendToProcess.push_back(vector<int>());
            RecvFromProcess.push_back(vector<int>());
        }

        for (int i = 0; i < N0; ++i){
            for (int j = IA[i]; j < IA[i + 1]; ++j){
                int p = Part[JA[j]];
                if (p < 0){
                    SendToProcess[-(p + 1)].push_back(i);
                    RecvFromProcess[-(p+1)].push_back(JA[j]);
                }
            }
        }

        for (int p = 0; p < P; ++p){
            if (!SendToProcess[p].empty()){
                Neighbours.push_back(p);

                vector<int> globalSend, globalRecv;
                for (int i = 0; i < SendToProcess[p].size(); ++i) {
                    int idx = SendToProcess[p][i];
                    globalSend.push_back(L2G[idx]);
                }
                for (int i = 0; i < RecvFromProcess[p].size(); ++i) {
                    int idx = RecvFromProcess[p][i];
                    globalRecv.push_back(L2G[idx]);
                }

                sort(globalSend.begin(), globalSend.end());
                sort(globalRecv.begin(), globalRecv.end());
                globalSend.erase(unique(globalSend.begin(), globalSend.end()), globalSend.end());
                globalRecv.erase(unique(globalRecv.begin(), globalRecv.end()), globalRecv.end());

                vector<int> localSend, localRecv;
                for (int i = 0; i < globalSend.size(); ++i) {
                    int globalIdx = globalSend[i];
                    localSend.push_back(G2L[globalIdx]);
                }
                for (int i = 0; i < globalRecv.size(); ++i) {
                    int globalIdx = globalRecv[i];
                    localRecv.push_back(G2L[globalIdx]);
                }

                SendToProcess[p] = move(localSend);
                RecvFromProcess[p] = move(localRecv);
            }
        }

        int send_ind = 0, recv_ind = 0;

        for (int p = 0; p < Neighbours.size(); ++p){
            SendOffset.push_back(send_ind);
            RecvOffset.push_back(recv_ind);
            for (int i = 0; i < SendToProcess[p].size(); ++i){
                Send.push_back(SendToProcess[p][i]);
                ++send_ind;
            }
            for(int i = 0; i < RecvFromProcess[p].size(); ++i){
                Recv.push_back(RecvFromProcess[p][i]);
                ++recv_ind;
            }
        }   

    }
};

template <typename VarType /* тип значений */> void Update(VarType *V, // Входной массив значений в вершинах/ячейках, который надо обновить
tCommScheme &CS/*какая-то структура, описывающая схему обменов*/){
const int B = CS.GetNumOfNeighbours(); // число соседей
    if(B==0) return; // нет соседей - нет проблем
    // tCommScheme - какая-то структура, замените ее на ваш вариант
    // приведем все к POD типам и неймингу, как было в тексте выше
    int *Send = CS.GetSendList(); // список ячеек на отправку по всем соседям
    int *Recv = CS.GetRecvList(); // список ячеек на прием по всем соседям
    int *SendOffset = CS.GetSendOffset(); // смещения списков по каждому соседу на отправку
    int *RecvOffset = CS.GetRecvOffset(); // смещения списков по каждому соседу на прием
    int *Neighbours = CS.GetListOfNeigbours(); // номера процессов соседей
    MPI_Comm MCW = CS.GetMyComm(); // коммуникатор для данной группы (MPI_COMM_WORLD)
    int sendCount=SendOffset[B]; // размер общего списка на отправку по всем соседям
    int recvCount=RecvOffset[B]; // размер общего списка на прием по всем соседям
    // MPI данные - сделаем статиками, поскольку это высокочастотная функция,
    // чтобы каждый раз не реаллокать (так делать вовсе не обязательно).
    static vector<VarType> SENDBUF, RECVBUF; // буферы на отправку и прием по всем соседям
    static vector<MPI_Request> REQ; // реквесты для неблокирующих обменов
    static vector<MPI_Status> STS; // статусы для неблокирующих обменов
    // ресайзим, если надо
    if(2*B > (int)REQ.size()){ REQ.resize(2*B); STS.resize(2*B); }
    if(sendCount>(int)SENDBUF.size()) SENDBUF.resize(sendCount);
    if(recvCount>(int)RECVBUF.size()) RECVBUF.resize(recvCount);
    int nreq=0; // сквозной счетчик реквестов сообщений
    // инициируем получение сообщений
    for(int p=0; p<B; p++){
        int SZ = (RecvOffset[p+1]-RecvOffset[p])*sizeof(VarType); // размер сообщения
        if(SZ<=0) continue; // если нечего слать - пропускаем соседа
        int NB_ID = Neighbours[p]; // узнаем номер процесса данного соседа
        int mpires = MPI_Irecv(&RECVBUF[RecvOffset[p]*sizeof(VarType)], SZ, MPI_CHAR,
        NB_ID, 0, MCW, &(REQ[nreq]));
        ASSERT(mpires==MPI_SUCCESS, "MPI_Irecv failed");
        //ASSERT - какой-то макрос проверки-авоста, замените на ваш способ проверки
        nreq++;
    }
    // пакуем данные с интерфейса по единому списку сразу по всем соседям
    #pragma omp parallel for // в параллельном режиме с целью ускорения (К.О.)
    for(int i=0; i<sendCount; ++i){
        SENDBUF[i] = V[Send[i]/*номер ячейки на отправку*/];
    }
    // инициируем отправку сообщений
    for(int p=0; p<B; p++){
        int SZ =(SendOffset[p+1]-SendOffset[p])*sizeof(VarType); // размер сообщения
        if(SZ<=0) continue; // если нечего принимать - пропускаем соседа
        int NB_ID = Neighbours[p]; // узнаем номер процесса данного соседа
        int mpires = MPI_Isend(&SENDBUF[SendOffset[p]*sizeof(VarType)], SZ, MPI_CHAR,
        NB_ID, 0, MCW, &(REQ[nreq]));
        ASSERT(mpires==MPI_SUCCESS, "MPI_Isend failed");
        nreq++;
    }
    if(nreq>0){ // ждем завершения всех обменов
    int mpires = MPI_Waitall(nreq, &REQ[0], &STS[0]);
    ASSERT(mpires==MPI_SUCCESS, "MPI_Waitall failed");
    }
    // разбираем данные с гало ячеек по единому списку сразу по всем соседям
    #pragma omp parallel for
    for(int i=0; i<recvCount; ++i) 
    {
        V[Recv[i]/*номер ячейки на прием*/] = RECVBUF[i];
    }
}


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
    #pragma omp parallel reduction(+:out)
    {
        #pragma omp for
        for (int i = 0; i < N; ++i){
            out += x[i] * y[i];
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

double solve(int N, int N0, int NumProc, tCommScheme Com, int*& IA, int*& JA, double*& A, double*& b, double eps, int maxit, double*& x, int &n){
    int k = 0;
    double *x_k, *r_k, *M, *z_k, *p_k, *q_k, *x_k_prev, *p_k_prev, *r_k_prev, *tmp, ro_k, ro_k_local, ro_k_prev, beta_k, alpha_k, alpha_k_local, t, t_collect, t_collect_local;
    x_k = new double[N];
    r_k = new double[N];
    M = new double[IA[N0]];
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
    
    t_collect_local = 0;
    for(int i = 0; i < 100; ++i){
        t = omp_get_wtime();
        SpMv(N, IA, JA, M, r_k_prev, z_k);
        t = omp_get_wtime() - t;
        t_collect_local += t;
    }
    t_collect_local /= 100.0;
    MPI_Allreduce(&t_collect_local, &t_collect, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cout << "SpMv took: " << setprecision(5) << t_collect / NumProc << " seconds" << endl;

    t_collect_local = 0;
    for(int i = 0; i < 1000; ++i){
        t = omp_get_wtime();
        axpy(N, 1.23, r_k_prev, x_k_prev, p_k);
        t = omp_get_wtime() - t;
        t_collect_local += t;
    }
    t_collect_local /= 1000.0;
    MPI_Allreduce(&t_collect_local, &t_collect, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cout << "axpy took: " << setprecision(5) << t_collect / NumProc << " seconds" << endl;

    t_collect_local = 0;
    for(int i = 0; i < 100; ++i){
        t = omp_get_wtime();
        dot(N, r_k_prev, x_k_prev);
        t = omp_get_wtime() - t;
        t_collect_local += t;
    }
    t_collect_local /= 100.0;
    MPI_Allreduce(&t_collect_local, &t_collect, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    cout << "dot took: " << setprecision(5) << t_collect / NumProc << " seconds" << endl;


    do{
        cout << "TEST SOLVE 0" << endl;
        ++k;
        vector_fill(IA[N0], M, 0);
        get_M(N0, IA, JA, A, M);
        cout << "TEST SOLVE 00" << endl;
        Update(&r_k_prev, Com);
        SpMv(N0, IA, JA, M, r_k_prev, z_k);
        ro_k_local = dot(N, r_k_prev, z_k);
        cout << "TEST SOLVE 0000" << endl;
        MPI_Allreduce(&ro_k_local, &ro_k, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        cout << "TEST SOLVE 1" << endl;

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

        cout << "TEST SOLVE 2" << endl;

        Update(&p_k, Com);
        SpMv(N0, IA, JA, A, p_k, q_k);
        alpha_k_local = ro_k / dot(N, p_k, q_k);
        MPI_Allreduce(&alpha_k_local, &alpha_k, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //vector_cp(N, tmp, x_k);

        axpy(N, alpha_k, p_k, x_k_prev, x_k);
        vector_cp(N, x_k_prev, x_k);
        axpy(N, -alpha_k, q_k, r_k_prev, r_k);
        vector_cp(N, r_k_prev, r_k);

        ro_k_prev = ro_k;

        cout << "TEST SOLVE 3" << endl;
    }
    while(ro_k > eps * eps && k < maxit);

    // TEST
    // cout << "Number of iterations: " << k <<endl;

    x = new double[N];
    vector_cp(N, x, x_k);
    n = k;

    return l2(N, r_k);
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

/*
i = element of i-th process, i >= 0
i = halo of (-(i+1)-th process, i < 0)
*/

int setHalo(int& N, int new_I, int*& Part, bool*& checkedHalo, int*& L2G, int*& G2L){
    if (Part[new_I] < 0 && checkedHalo[new_I] == false){ // halo
        L2G[N] = new_I;
        G2L[new_I] = N;
        ++N;
        checkedHalo[new_I] = true;    
    }
    return G2L[new_I];
}

void countHalo(int p_id, int& N_halo, int base_ind, int ind, int*& Part){
    if (Part[base_ind] != p_id){
        return;
    }
    if (Part[ind] != p_id){
        if (Part[ind] >= 0){
            Part[ind] = -Part[ind] - 1;
            ++N_halo;
        }
    }
}

void generate(int p_id, int Nx, int Ny, int K1, int K2, int Px, int Py, int& N, int& N0, int*& Part, int*& L2G, int*& G2L, int*& IA, int*& JA){
    int P = Px * Py;
    int K = K1 + K2;
    int trian_cnt = ((Ny * Nx) / K) * K2 + max(0, (Ny * Nx) % K - K1);
    int doubled_E = 0;
    int N_all = Ny * Nx + trian_cnt;
    int N_halo = 0;

    //int* v_types = new int[N];
    int* cnt_neigh = new int[N_all];

    // MPI SPECIFIC - FILL PART
    N = 0;
    N0 = 0;
    Part = new int[N_all];
    bool* checkedHalo = new bool[N_all];
    for(int i = 0; i < N_all; ++i){
        checkedHalo[i] = false;
    }

    int process_size = Ny * Nx / P;

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int old_I = i * Nx + j;
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            int cur_p = old_I / process_size;
            if (cur_p >= P){
                cur_p = old_I % P;
            }
            Part[new_I] = cur_p;
            if (cur_p == p_id){
                ++N0;
            }

            if (old_I % K >= K1){
                Part[new_I + 1] = Part[new_I];
                if (Part[new_I] == p_id){
                    ++N0;
                }
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

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            if ((i * Nx + j) % K >= K1){
                //v_types[new_I] = 1; // верхний треугольник
                //v_types[new_I + 1] = 2; // нижний треугольник

                cnt_neigh[new_I] = 2;
                if (i - 1 >= 0){ // верхний сосед
                    ++cnt_neigh[new_I];
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part);
                }
                if (j - 1 >= 0){ // левый сосед
                    ++cnt_neigh[new_I];
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i, j - 1), Part);
                }
                doubled_E += cnt_neigh[new_I] * (Part[new_I] == p_id ? 1 : 0);

                cnt_neigh[new_I + 1] = 2;
                if (i + 1 < Ny){ // нижний сосед
                    ++cnt_neigh[new_I + 1]; 
                    countHalo(p_id, N_halo, new_I + 1, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part);
                }
                if (j + 1 < Nx){
                    ++cnt_neigh[new_I + 1]; // правый сосед
                    countHalo(p_id, N_halo, new_I + 1, oldInd2New(Nx, Ny, K1, K2, i, j + 1), Part);
                }
                doubled_E += cnt_neigh[new_I + 1] * (Part[new_I + 1] == p_id ? 1 : 0);
            }
            else{
                //v_types[new_I] = 0; // обычная клетка
                cnt_neigh[new_I] = 1;
                if (i - 1 >= 0){ // верхний сосед
                    ++cnt_neigh[new_I];
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part);
                }
                if (i + 1 < Ny){ // нижний сосед
                    ++cnt_neigh[new_I]; 
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part);
                }
                if (j - 1 >= 0){ // левый сосед
                    ++cnt_neigh[new_I];
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i, j - 1), Part);
                }
                if (j + 1 < Nx){
                    ++cnt_neigh[new_I]; // правый сосед
                    countHalo(p_id, N_halo, new_I, oldInd2New(Nx, Ny, K1, K2, i, j + 1), Part);
                }
                doubled_E += cnt_neigh[new_I] * (Part[new_I] == p_id ? 1 : 0);
            }
        }
    }

    cout << endl << endl << "p_id: " << p_id << endl;
    cout << "N_halo: " << N_halo << endl;
    cout << N0 << endl;

    IA = new int[N0 + 1];
    L2G = new int[N0 + N_halo];
    G2L = new int[N_all];
    JA = new int[doubled_E];

    #pragma omp parallel for
    for(int i = 0; i < N_all; ++i){
        G2L[i] = -1;
    }

    int prev_global, local_ind = 0;
    for(int i = 0; i < N_all; ++i){
        if (Part[i] == p_id){
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
    cout << "local_ind: " << local_ind << endl;
    cout << "doubled_E: " << doubled_E << endl;
    IA[local_ind] = doubled_E;

    cout << "TEST SET 1" << endl;

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            int new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            int old_I = i * Nx + j;
            if (Part[new_I] != p_id){
                continue;
            }
            int cur_JA = IA[G2L[new_I]];
            if (old_I % K >= K1){
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1, Part, checkedHalo, L2G, G2L);
                        //JA[cur_JA] = G2L[oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1];
                        ++cur_JA;
                    }else{
                        //JA[cur_JA] = G2L[oldInd2New(Nx, Ny, K1, K2, i - 1, j)];
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part, checkedHalo, L2G, G2L);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = setHalo(N, new_I - 1, Part, checkedHalo, L2G, G2L);
                    ++cur_JA;
                }
                JA[cur_JA] = setHalo(N, new_I, Part, checkedHalo, L2G, G2L); // главная диагональ
                ++cur_JA;
                JA[cur_JA] = setHalo(N, new_I + 1, Part, checkedHalo, L2G, G2L); // нижний треугольник = правый сосед
                ++cur_JA;
                
                JA[cur_JA] = setHalo(N, new_I, Part, checkedHalo, L2G, G2L); // верхний треугольник = левый сосед
                ++cur_JA;
                JA[cur_JA] = setHalo(N, new_I + 1, Part, checkedHalo, L2G, G2L); // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = setHalo(N, new_I + 2, Part, checkedHalo, L2G, G2L);
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part, checkedHalo, L2G, G2L);
                    ++cur_JA;
                }            
            }
            else{
                if (i - 1 >= 0){ // верхний сосед
                    if (((i - 1) * Nx + j) % K >= K1){ // если сверху треугольник
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1, Part, checkedHalo, L2G, G2L);
                        ++cur_JA;
                    }else{
                        JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i - 1, j), Part, checkedHalo, L2G, G2L);
                        ++cur_JA;
                    }
                }
                if (j - 1 >= 0){ // левый сосед
                    JA[cur_JA] = setHalo(N, new_I - 1, Part, checkedHalo, L2G, G2L);
                    ++cur_JA;
                }
                JA[cur_JA] = setHalo(N, new_I, Part, checkedHalo, L2G, G2L); // главная диагональ
                ++cur_JA;
                if (j + 1 < Nx){ // правый сосед
                    JA[cur_JA] = setHalo(N, new_I + 1, Part, checkedHalo, L2G, G2L);
                    ++cur_JA; 
                }
                if (i + 1 < Ny){ // нижний сосед
                    JA[cur_JA] = setHalo(N, oldInd2New(Nx, Ny, K1, K2, i + 1, j), Part, checkedHalo, L2G, G2L);
                    ++cur_JA;
                }
            }
        }
    }

    // TEST
    int cnt = 0;
    for(int i = 0; i < N_all; ++i){
        if (Part[i] < 0 && checkedHalo[i] == false){
            ++cnt;
        }
    }
    cout << "HALO TEST: " << cnt;

    cout << "N: " << N << endl;
    cout << "N0: " << N0 << endl;
    cout << "N_halo: " << N_halo << endl;

    N = N0 + N_halo;

    // if (N0 + N_halo != N){
    //     cout << "N: " << N << endl;
    //     cout << "N0: " << N0 << endl;
    //     cout << "N_halo: " << N_halo << endl;
    // }

}

void fill(int N, int N0, int*& IA, int*& JA, int*& L2G, double*& A, double*& b){
    double DIAG_COEFF = 1.234;

    int* diag = new int[N0];
    A = new double[IA[N0]];
    b = new double[N];

    // cout << "TEST 1" << endl;
    // cout << "L2G: " << L2G[N - 1] << endl;
    // cout << "JA: " << JA[IA[N0] - 1] << endl;
    // cout << "IA[N0]: " << IA[N0] << endl;
    
    #pragma omp parallel for
    for (int i = 0; i < N0; ++i){
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            if (i == JA[j]){
                diag[i] = j;
                A[j] = 0;
                continue;
            }
            A[j] = cos(i * L2G[JA[j]] + i + L2G[JA[j]]);
        }
    }

    // cout << "TEST 2" << endl;

    #pragma omp parallel for
    for (int i = 0; i < N0; ++i){
        for (int j = IA[i]; j < IA[i + 1]; ++j){
            if (i == JA[j]){
                continue;
            }
            A[diag[i]] += abs(A[j]);
        }
    }

    // cout << "TEST 3" << endl;

    #pragma omp parallel for
    for (int i = 0; i < N0; ++i){
        A[diag[i]] *= DIAG_COEFF;
    }
    
    // cout << "TEST 4" << endl;

    #pragma omp parallel for
    for (int i = 0; i < N; ++i){
        b[i] = sin(L2G[i]);
    }

    // cout << "TEST 5" << endl;

    delete [] diag;
}

tCommScheme com(int N, int N0, int P, int*& IA, int*& JA, int*& Part, int*& L2G, int*& G2L){
    tCommScheme out = tCommScheme(MPI_COMM_WORLD);
    out.Initialize(N, N0, P, IA, JA, Part, L2G, G2L);
    return out;
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

    cout << "NUM PROC: " << NumProc << endl;

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
    if (Px * Py != NumProc){
        cout << "Invalid process amount!" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        return 0;
    }

    omp_set_num_threads(T);

    t = omp_get_wtime();
    generate(MyID, Nx, Ny, K1, K2, Px, Py, N, N0, Part, L2G, G2L, IA, JA);
    t = omp_get_wtime() - t;
    cout << setprecision(5) << "Generate took: " << t << " seconds" << endl;

    logFile << "N0:" << N0 << endl;

    t = omp_get_wtime();
    fill(N, N0, IA, JA, L2G, A, b);
    t = omp_get_wtime() - t;
    cout << setprecision(5) << "Fill took: " << t << " seconds" << endl;

    t = omp_get_wtime();
    tCommScheme Com = com(N, N0, NumProc, IA, JA, Part, L2G, G2L);
    t = omp_get_wtime() - t;
    cout << setprecision(5) << "com took: " << t << " seconds" << endl;

    t = omp_get_wtime();
    res = solve(N, N0, NumProc, Com, IA, JA, A, b, EPS, MAXIT, x, n);
    t = omp_get_wtime() - t;

    cout << "Solve took: " << setprecision(5) << t << " seconds" << endl;
    // logFile << "start free" << endl;
    // delete [] IA;
    // logFile << "good IA" << endl;
    // delete [] JA;
    // logFile << "good JA" << endl;
    // delete [] Part;
    // logFile << "good PART" << endl;
    // delete [] L2G; 
    // logFile << "good L2G" << endl;
    // delete [] G2L;
    // logFile << "good G2L" << endl;
    // delete [] A;
    // delete [] b;
    // delete [] x;
    logFile.close();
    return 0;
}