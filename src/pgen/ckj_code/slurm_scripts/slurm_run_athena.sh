#!/bin/bash
#SBATCH -J slurm_run_athena          # 任务名称
#SBATCH -o log/Athena++.%j.out       # 标准输出文件
#SBATCH --qos=low                    # 质量服务等级
#SBATCH --nodes=1                    # 分配的节点数
#SBATCH --ntasks-per-node=64          # 每个节点的任务数
#SBATCH --ntasks=64                   # 总任务数


# 导入MPI运行环境
module purge
export MPICH_CXX=icpx
module load compiler hdf5/1.12.1-p-oneapi_2023.0 mpi fftw
module list

# 输出一些信息
#TODO 输出当前的日期、时间、主机名、用户、SLURM相关的信息
echo "job id: $SLURM_JOB_ID"
echo "job name: $SLURM_JOB_NAME"
echo "number of tasks: $SLURM_NTASKS"
echo "number of nodes: $SLURM_JOB_NUM_NODES"
echo "working path = $(pwd)"
echo "PATH = $PATH"

# 检查是否只有一个 *.athinput 文件

athinput_files=( *.athinput )   # 使用 globbing 找到所有匹配 *.athinput 的文件，并将它们放入一个数组中
if (( ${#athinput_files[@]} != 1 )); then  # 检查数组的长度
    echo "Error: Expected exactly one *.athinput file, but found ${#athinput_files[@]}"
    exit 1
fi
# 如果只有一个匹配的文件，脚本将继续执行
athinput_file="${athinput_files[0]}"
echo "Found exactly one *.athinput file: ${athinput_file}"


export UCX_TLS=ud,sm,self # https://github.com/openucx/ucx/issues/4742，解决 OFI 大规模时报错
export I_MPI_PMI_LIBRARY=/lib64/libpmi2.so  # 设置 PMI 的路径。也可以 /lib64/libpmi.so，不过 SLURM 官网推荐 pmi2
export SLURM_MPI_TYPE=pmi2 # 等价于 srun 的时候加入参数 --mpi=pmi2。不加会报错 PMI2_Job_GetId returned 14

# 生成 machinefile
srun hostname -s | sort -n > slurm.hosts


# 保存一些信息
# athena -d output -i test_cooling.athinput time/tlim=0
mkdir -p info                                              # 创建 info 路径用于存放各种信息
mpirun -n 1 athena -c > info/configure.txt                 # 打印 Athena++ 的配置信息。MPI 并行时，需要用 mpirun，否则单独 athena -c 的话 MPI 会报错。
mpirun -n 1 athena -d output -i "$athinput_file" -m "$SLURM_NTASKS" > info/MeshBlocks.txt 
mv mesh_structure.dat info/                                # 保存网格结构信息

module load vtune

MPI_CMD="mpirun -n \"$SLURM_NTASKS\" -machinefile slurm.hosts"

# VTune 的命令。
# 不过注意：对于长时间的模拟，加上 vtune 可能会在结束时 MPI_Finalize 时报错 Segmentation fault，暂时不知道怎么解决。
VTUNE_CMD=""
# VTUNE_CMD="vtune -collect hotspots -r vtune_result/athena -trace-mpi -data-limit=200 -target-duration-type=short"   # 注释掉这一行，即可关闭 vtune。# 其他可选的选项： -duration 60

ATHENA_CMD="athena -d output -i \"$athinput_file\""
# ATHENA_CMD="athena -d output -i \"$athinput_file\" -r \"Bondi.00001.rst\""

# 执行MPI并行计算程序
echo "[Running Athena++]: "
# mpirun -n "$SLURM_NTASKS" -machinefile slurm.hosts athena -d output -i "$athinput_file"
# mpirun -n "$SLURM_NTASKS" -machinefile slurm.hosts athena -d output -i "$athinput_file" -r "Bondi.00001.rst"
# srun athena -d output -i test_cooling.athinput #FIXME 不对，这个不能跨节点？

# 组装命令并执行
RUN_COMMAND="$MPI_CMD $VTUNE_CMD $ATHENA_CMD"
echo "[Running command]: $RUN_COMMAND"
eval "$RUN_COMMAND"


# 运行完之后，直接做后处理

VTUNE_CMD=""
# VTUNE_CMD="vtune -collect hotspots -knob sampling-mode=sw -r vtune_result/post_processing -trace-mpi"   # 注释掉这一行，即可关闭 vtune。# 其他可选的选项： -duration 60

# PYTHON_SCRIPT="$HOME/Codes/athena_post_processing/post_processing/scripts/single_SN.py"
PYTHON_SCRIPT="$HOME/Codes/athena_post_processing/post_processing/scripts/Bondi+random_SN.py"
PYTHON_CMD="python \"$PYTHON_SCRIPT\""

# 组装命令并执行
echo "[Running post-processing script]: $PYTHON_SCRIPT"
Post_Processing_COMMAND="$MPI_CMD $VTUNE_CMD $PYTHON_CMD"
echo "[Running command]: $Post_Processing_COMMAND"
eval "$Post_Processing_COMMAND"

# # srun --jobid=1948647 -N 24 -n 256  athena -d output -i Bondi+SAMR.athinput -r output/Bondi.00000.rst