#!/bin/bash
#SBATCH -J slurm_run_athena          # 任务名称
#SBATCH -o log/Athena++.%j.out       # 标准输出文件
#SBATCH --qos=low                    # 质量服务等级
#SBATCH --nodes=1                    # 分配的节点数
#SBATCH --ntasks-per-node=64          # 每个节点的任务数
#SBATCH --ntasks=64                  # 总任务数

# 临时启用的选项
##SBATCH --partition=C064M1024G        # 设置分区


# 以下不启用，留作参考，用多个 # 号注释掉
##SBATCH --partition=C064M1024G        # 设置分区
##SBATCH --nodelist=node1,node2,node3  # 指定节点
##SBATCH --exclude=node4,node5         # 排除节点


# 输出一些信息
echo "====================== Info ======================="
echo "当前时间: $(date)"
echo "主机名: $(hostname)"
echo "用户名: $(whoami)"
echo "运行脚本: $0"
echo "工作目录: $(pwd)"
echo "PATH: $PATH"
echo
echo "Slurm 相关信息:"
echo "作业 ID: $SLURM_JOB_ID"
echo "作业名称: $SLURM_JOB_NAME"
echo "提交目录: $SLURM_SUBMIT_DIR"
echo "分配的节点数: $SLURM_JOB_NUM_NODES"
echo "每个节点的任务数: $SLURM_NTASKS_PER_NODE"
echo "总任务数: $SLURM_NTASKS"
echo "节点列表: $SLURM_JOB_NODELIST"
echo "分区: $SLURM_JOB_PARTITION"
echo "每个 CPU 的内存: $SLURM_MEM_PER_CPU MB"
echo "每个节点的内存: $SLURM_MEM_PER_NODE"
echo "作业账户: $SLURM_JOB_ACCOUNT"
echo "作业 QOS: $SLURM_JOB_QOS"
echo "===================================================="
echo  # 输出一个空行


# 重置 SECONDS 变量，用于计时
SECONDS=0

# 获取脚本开始时间
start_time=$(date +%s)


# 检查是否只有一个 *.athinput 文件
athinput_files=( *.athinput )   # 使用 globbing 找到所有匹配 *.athinput 的文件，并将它们放入一个数组中
if (( ${#athinput_files[@]} != 1 )); then  # 检查数组的长度
    echo "Error: Expected exactly one *.athinput file, but found ${#athinput_files[@]}"
    exit 1
fi
# 如果只有一个匹配的文件，脚本将继续执行
athinput_file="${athinput_files[0]}"
echo "Found exactly one *.athinput file: ${athinput_file}"



# 导入MPI运行环境
module purge
export MPICH_CXX=icpx
module load compiler hdf5/1.12.1-p-oneapi_2023.0 mpi fftw
module load vtune
module list

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


MPI_CMD="mpirun -n \"$SLURM_NTASKS\" -machinefile slurm.hosts"  #FUTURE 可以尝试用 srun。之前的尝试，似乎不能跨节点？

# VTune 的命令。
# 不过注意：对于长时间的模拟，加上 vtune 可能会在结束时 MPI_Finalize 时报错 Segmentation fault，暂时不知道怎么解决。
VTUNE_CMD=""
# VTUNE_CMD="-gtool \"vtune -collect hotspots -r vtune_result/Athena++.$SLURM_JOB_ID -trace-mpi -data-limit=200 -target-duration-type=long : 0-1 -- \" "   # 注释掉这一行，即可关闭 vtune。# 其他可选的选项： -duration 60

ATHENA_CMD="athena -d output -i \"$athinput_file\""
# ATHENA_CMD="$ATHENA_CMD -r \"output/Bondi+SN.final.rst\""   # 从某个 rst 文件开始继续模拟


# 组装命令并执行
RUN_COMMAND="$MPI_CMD $VTUNE_CMD $ATHENA_CMD"
echo "[Running Athena++]: "
echo "开始运行 Athena++ 的时间: $SECONDS 秒"
echo "[Running command]: $RUN_COMMAND"
eval "$RUN_COMMAND"

echo "Athena++ 运行完毕"
echo "当前时间: $(date)"
echo "当前脚本耗时: $SECONDS 秒"



# 运行完之后，直接做后处理
#TODO 可以考虑调用 post_processing.sh，但要考虑到 SLURM 环境变量、MPI_CMD 等的传递问题

VTUNE_CMD=""
# VTUNE_CMD="vtune -collect hotspots -knob sampling-mode=sw -r vtune_result/post_processing.$SLURM_JOB_ID -trace-mpi"   # 注释掉这一行，即可关闭 vtune。# 其他可选的选项： -duration 60

# PYTHON_SCRIPT="$HOME/Codes/athena_post_processing/post_processing/scripts/single_SN.py"
PYTHON_SCRIPT="$HOME/Codes/athena_post_processing/post_processing/scripts/Bondi+random_SN.py"
PYTHON_CMD="python \"$PYTHON_SCRIPT\""

# 组装命令并执行
echo "[Running post-processing script]: $PYTHON_SCRIPT"
Post_Processing_COMMAND="$MPI_CMD $VTUNE_CMD $PYTHON_CMD"
echo "[Running command]: $Post_Processing_COMMAND"
eval "$Post_Processing_COMMAND"

# # srun --jobid=1948647 -N 24 -n 256  athena -d output -i Bondi+SAMR.athinput -r output/Bondi.00000.rst


# 输出脚本运行时间
echo
echo "脚本运行完毕"
echo "当前时间: $(date)"
echo "运行时长: $SECONDS 秒"