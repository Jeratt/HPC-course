Запуск программы должен осуществляться в виде \<program\> Nx Ny K1 K2 T, где Nx, Ny - параметры расчётной сетки, K1, K2 - параметры разбиения на треугольники, T - число потоков

Задача помещается в очередь с помощью .lsf скриптов. Для их создания для удобства прикладывается скрипт script_maker.py

При запуске без/с неверными параметрами будет выведен help с инфомрмацией, идентичной первой строке данного файла

Компилляция - g++ SLAU_solver.cpp -fopenmp -o a.out

Отладочная информация при наличии помещается в output_log.txt
При желании раскомментируйте необходимые секции для их печати в отладку, например (410-412 строки SLAU_solver.cpp) для печати вектора-ответа:
```
    // for(int i = 0; i < N; ++i){
    //     logFile << x[i] << " ";
    // }
```