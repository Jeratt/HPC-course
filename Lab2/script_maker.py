N_grid = [100, 300, 1000, 3000]
T_grid = [1, 2, 4, 8, 16, 32]
p_grid = [1, 2, 4, 8, 16, 32]

for n in N_grid:
    for p in p_grid:
        with open(f"solver_job_{n}_{p}.lsf", 'w') as fin:
            if p > 1:
                fin.write(f"#BSUB -n {p//2} -R \"span[ptile=2]\"\n")
            else:
                fin.write(f"#BSUB -n {p}\n")
            fin.write("#BSUB -J \"OpenMP_job\"\n")
            fin.write("#BSUB -o \"OpenMP_job%J.out\"\n")
            fin.write("#BSUB -e \"OpenMP_job%J.err\"\n")
            fin.write(f"#BSUB -R \"affinity[core({1})]\"\n")
            if p == 1:
                fin.write(f"mpiexec -np 1 ../a.out {n} {n} 100 1 1 1 {1}\n")
            elif p == 2:
                fin.write(f"mpiexec ../a.out {n} {n} 100 1 1 2 {1}\n");
            elif p == 4:
                fin.write(f"mpiexec ../a.out {n} {n} 100 1 2 2 {1}\n");
            elif p == 8:
                fin.write(f"mpiexec ../a.out {n} {n} 100 1 4 2 {1}\n");
            elif p == 16:
                fin.write(f"mpiexec ../a.out {n} {n} 100 1 4 4 {1}\n");