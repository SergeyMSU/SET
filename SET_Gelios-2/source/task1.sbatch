#! /bin/bash
#SBATCH --job-name="Hello"
#SBATCH --gpus=0
#SBATCH --nodes=1                      # Количество используемых узлов
#SBATCH --ntasks=1                     # Количество MPI процессов
#SBATCH --cpus-per-task=20             # Количество ядер CPU, выделенных для одного процесса
#SBATCH --time=0-10:00
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_b"

# Executable
srun ./a.out
