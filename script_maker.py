N_grid = [100, 300, 1000, 3000]
T_grid = [1, 2, 4, 8, 16]

for n in N_grid:
    for t in T_grid:
        with open(f"solver_job_{n}_{t}.lsf", 'w') as fin:
            fin.write("#BSUB -J \"OpenMP_job\"\n")
            fin.write("#BSUB -o \"OpenMP_job%J.out\"\n")
            fin.write("#BSUB -e \"OpenMP_job%J.err\"\n")
            if t < 16:
                fin.write(f"#BSUB -R \"affinity[core({t})]\"\n")
            else:
                fin.write(f"#BSUB -R \"affinity[core(8)]\"\n")
                fin.write(f"OMP_NUM_THREADS={t}\n")
            fin.write(f"/polusfs/lsf/openmp/launchOpenMP.py ../a.out {n} {n} 3 3 {t}")