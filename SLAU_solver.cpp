#include <algorithm>
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

void generate(int Nx, int Ny, int K1, int K2, int& N, int** IA, int** JA){
    N = oldInd2New(Nx, Ny, K1, K2, Ny - 1, Nx - 1);
    int doubled_E = N; // сразу плюсуем главную диагональ
    int K = K1 + K2;
    int cnt_IA = 0, cnt_JA = 0;
    int new_I;

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

    *IA = new int[N + 2]; // + 2 т.к. ещё резервируем под последний элемент, равный размеру JA
    *JA = new int[doubled_E];

    for(int i = 0; i < Ny; ++i){
        for(int j = 0; j < Nx; ++j){
            new_I = oldInd2New(Nx, Ny, K1, K2, i, j);
            *IA[cnt_IA] = cnt_JA;
            if (i - 1 >= 0) {
                if (((i - 1) * Nx + j) % K >= K1){ // если верхний сосед разбивается на треугольник - соеднинение с правым треугольником
                    *JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j) + 1;
                }
                else{
                    *JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i - 1, j);
                }
                ++cnt_JA;
            }

            if (j - 1 >= 0) { // левый сосед
                *JA[cnt_JA] = new_I - 1;
                ++cnt_JA;
            }

            // добавление связи "с собой" - главная диагональ
            *JA[cnt_JA] = new_I;
            ++cnt_JA;

            if (j + 1 < Nx) { // правый сосед
                *JA[cnt_JA] = new_I + 1;
                ++cnt_JA;
            }

            if (i + 1 < Ny) { // нижний сосед
                *JA[cnt_JA] = oldInd2New(Nx, Ny, K1, K2, i + 1, j);
                ++cnt_JA;
            }
            ++cnt_IA;
        }
    }
    *IA[cnt_IA] = doubled_E; // последний элемент IA - число ненулевых в портрете
}



int main(int argc, char** argv){
    
    return 0;
}