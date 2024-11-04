#include <algorithm>
#include <iostream> // TEST
#include <fstream>

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
    N = oldInd2New(Nx, Ny, K1, K2, Ny - 1, Nx - 1) + 1;
    if ((Ny - 1) * Nx + (Nx - 1) % K >= K1){
        ++N;
    }
    int doubled_E = N; // сразу плюсуем главную диагональ
    cout << "doubled_E: " << doubled_E << endl;
    int cnt_IA = 0, cnt_JA = 0;
    int new_I;
    bool split;

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            if (i - 1 >= 0) {
                ++doubled_E;
            }
            if (i + 1 < Ny) {
                ++doubled_E;
            }
            if (j - 1 >= 0) {
                ++doubled_E;
            }
            if (j + 1 < Nx) {
                ++doubled_E;
            }
            if ((i * Nx + j) % K >= K1){
                doubled_E += 2; // если клетка разбивается - +2 ненулевых элемента
            }
        }
    }

    IA = new int[N + 1]; // + 2 т.к. ещё резервируем под последний элемент, равный размеру JA
    JA = new int[doubled_E];

    // cout << N + 2 << " " << doubled_E << endl;

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            split = (i * Nx + j) % K >= K1;
            if (split){ // если клетка разбивается на треугольники - обработать их вместе
                new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
                IA[cnt_IA] = cnt_JA;
                if (i - 1 >= 0) {
                    if (((i - 1) * Nx + j) % K >= K1){ // если верхний сосед разбивается на треугольник - соеднинение с правым треугольником
                        JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1;
                    }
                    else{
                        JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j);
                    }
                    ++cnt_JA;
                }
                if (j - 1 >= 0) { // левый сосед
                    JA[cnt_JA] = new_I - 1;
                    ++cnt_JA;
                }
                // добавление связи "с собой" - главная диагональ
                JA[cnt_JA] = new_I;
                ++cnt_JA;
                JA[cnt_JA] = new_I + 1; // правый сосед точно есть
                ++cnt_JA;
                ++cnt_IA;
                ++new_I;

                IA[cnt_IA] = cnt_JA;
                JA[cnt_JA] = new_I - 1; // левый сосед точно есть
                ++cnt_JA;
                // добавление связи "с собой" - главная диагональ
                JA[cnt_JA] = new_I;
                ++cnt_JA;
                if (j + 1 < Nx) { // правый сосед
                    JA[cnt_JA] = new_I + 1;
                    ++cnt_JA;
                }
                if (i + 1 < Ny) { // нижний сосед
                    JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i + 1, j);
                    ++cnt_JA;
                }
                ++cnt_IA;
            }
            else{
                new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
                IA[cnt_IA] = cnt_JA;
                if (i - 1 >= 0) {
                    if (((i - 1) * Nx + j) % K >= K1){ // если верхний сосед разбивается на треугольник - соеднинение с правым треугольником
                        JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1;
                    }
                    else{
                        JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j);
                    }
                    ++cnt_JA;
                }

                if (j - 1 >= 0) { // левый сосед
                    JA[cnt_JA] = new_I - 1;
                    ++cnt_JA;
                }

                // добавление связи "с собой" - главная диагональ
                JA[cnt_JA] = new_I;
                ++cnt_JA;

                if (j + 1 < Nx) { // правый сосед
                    JA[cnt_JA] = new_I + 1;
                    ++cnt_JA;
                }

                if (i + 1 < Ny) { // нижний сосед
                    JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i + 1, j);
                    ++cnt_JA;
                }
                ++cnt_IA;
            }
        }
    }
    cout <<"doubled_E: " << doubled_E << "  cnt_JA: " <<  cnt_JA << endl;

    // cout << cnt_IA << "   KEK    "<< cnt_JA << endl;
    IA[cnt_IA] = doubled_E; // последний элемент IA - число ненулевых в портрете
    // cout << "TEST 4" << endl;
}



int main(int argc, char** argv){
    int Nx, Ny, K1, K2, N, doubled_E;
    int *IA, *JA;
    ofstream logFile("error_log.txt", ios::out);
    if (!logFile.is_open()){
        return 1;
    }
    if (argc < 2){
        // TODO print help
    }
    ifstream inFile(argv[1], ios::in);
    inFile >> Nx >> Ny >> K1 >> K2;

    generate(Nx, Ny, K1, K2, N, IA, JA);
    for (int i = 0; i < N+1; ++i){
        logFile << IA[i] << " ";
    }
    logFile << endl;
    doubled_E = IA[N];
    cout << doubled_E << endl;
    for (int i = 0; i < doubled_E; ++i){
        logFile << JA[i] << " ";
    }  

    logFile.close();
    inFile.close();

    delete [] IA;
    cout << "TEST Y" << endl;
    delete [] JA;
    cout << "TEST X" << endl;
    return 0;
}