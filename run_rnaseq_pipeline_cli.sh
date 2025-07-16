#!/bin/bash

# =================================================================
# RNA-seq 并行分析流程脚本 (命令行版)
#
# 功能:
# 1. 使用 Trim Galore 对原始 FASTQ 文件进行质量修剪
# 2. 使用 HISAT2 将修剪后的 reads 比对到参考基因组
# 3. 使用 Samtools 处理比对文件
# 4. 使用 featureCounts 进行基因计数
#
# 特点:
# - 命令行参数: 通过参数灵活配置所有路径和资源。
# - 并行处理: 可同时处理多个样本以加快速度。
# - 自动路径: 关键路径已配置好。
#
# 作者: Cpg
# 日期: 2025-07-16
# =================================================================

# --- 1. 帮助信息和参数默认值 ---

usage() {
    echo "Usage: $0 -g <genome_index_prefix> -a <annotation_gtf> -i <input_dir> -o <output_dir> [-j <max_jobs>] [-n <threads_per_job>] [-s <steps_to_skip>]"
    echo ""
    echo "Options:"
    echo "  -g    [必须] HISAT2 基因组索引的前缀路径 (e.g., /path/to/ref/mm10)"
    echo "  -a    [必须] GTF 注释文件路径 (e.g., /path/to/annotation.gtf)"
    echo "  -i    [必须] 包含原始 FASTQ 文件的输入目录"
    echo "  -o    [必须] 用于存放所有结果的输出目录"
    echo "  -j    [可选] 并行处理的样本作业数 (默认: 1)"
    echo "  -n    [可选] 每个作业内部使用的线程数 (默认: 8)"
    echo "  -s    [可选] 跳过指定的步骤，用逗号分隔 (e.g., '1,4'). Steps: 1=TrimGalore, 2=HISAT2, 3=Samtools, 4=featureCounts"
    echo "  -h    显示此帮助信息"
    echo ""
    exit 1
}

# 默认值
MAX_JOBS=1
THREADS_PER_JOB=8
SKIP_STEPS=""

# --- 2. 解析命令行参数 ---

while getopts ":g:a:i:o:j:n:s:h" opt; do
    case ${opt} in
        g) REF_GENOME=$OPTARG ;;
        a) GTF_FILE=$OPTARG ;;
        i) RAW_DIR=$OPTARG ;;
        o) PROJ_DIR=$OPTARG ;;
        j) MAX_JOBS=$OPTARG ;;
        n) THREADS_PER_JOB=$OPTARG ;;
        s) SKIP_STEPS=$OPTARG ;;
        h) usage ;;
        \?) echo "无效选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数." >&2; usage ;;
    esac
done

# 检查所有必需的参数是否已提供
if [ -z "${REF_GENOME}" ] || [ -z "${GTF_FILE}" ] || [ -z "${RAW_DIR}" ] || [ -z "${PROJ_DIR}" ]; then
    echo "错误: 缺少必需的参数。"
    usage
fi

# --- 3. 设置目录 ---

CLEAN_DIR="${PROJ_DIR}/cleandata"
ALIGN_DIR="${PROJ_DIR}/alignment"
COUNT_DIR="${PROJ_DIR}/counts"
UNMAPPED_DIR="${PROJ_DIR}/unmapped_reads"

# 创建输出目录
echo "--- 创建输出目录 ---"
mkdir -p ${PROJ_DIR}
mkdir -p ${CLEAN_DIR}
mkdir -p ${ALIGN_DIR}
mkdir -p ${COUNT_DIR}
mkdir -p ${UNMAPPED_DIR}
echo "项目根目录: ${PROJ_DIR}"
echo "目录创建完成: cleandata, alignment, counts"
echo ""


# --- 4. 定义单个样本的处理函数 ---

process_sample() {
    local fq1=$1
    local fq2=$2
    local sample=$3
    
    echo "================================================="
    echo ">>> [PID: $$] 开始处理样本: ${sample}"
    echo "================================================="

    # --- 步骤 1: Trim Galore ---
    if [[ ! ",${SKIP_STEPS}," =~ ",1," ]]; then
        echo "[${sample}] 步骤 1: 使用 Trim Galore 进行质量修剪..."
        trim_galore --paired --fastqc --cores ${THREADS_PER_JOB} -o ${CLEAN_DIR} ${fq1} ${fq2} > /dev/null 2>&1
        echo "[${sample}] 修剪完成."
    else
        echo "[${sample}] 跳过步骤 1: Trim Galore"
    fi

    # 修正: Trim Galore 的输出文件名是基于输入文件名加上 _val_1.fq.gz / _val_2.fq.gz
    local base_fq1_name=$(basename ${fq1} .fastq.gz)
    local base_fq2_name=$(basename ${fq2} .fastq.gz)
    local clean_fq1="${CLEAN_DIR}/${base_fq1_name}_val_1.fq.gz"
    local clean_fq2="${CLEAN_DIR}/${base_fq2_name}_val_2.fq.gz"

    # --- 步骤 2: HISAT2 比对 ---
    if [[ ! ",${SKIP_STEPS}," =~ ",2," ]]; then
        echo "[${sample}] 步骤 2: 使用 HISAT2 进行比对..."
        hisat2 -p ${THREADS_PER_JOB} -x ${REF_GENOME} \
               -1 ${clean_fq1} \
               -2 ${clean_fq2} \
               -S ${ALIGN_DIR}/${sample}.sam \
               --un-conc-gz ${UNMAPPED_DIR}/${sample}.unmapped.fastq.gz
        echo "[${sample}] 比对完成."
    else
        echo "[${sample}] 跳过步骤 2: HISAT2"
    fi

    # --- 步骤 3: Samtools 转换和排序 ---
    if [[ ! ",${SKIP_STEPS}," =~ ",3," ]]; then
        echo "[${sample}] 步骤 3: 使用 Samtools 转换和排序 BAM..."
        samtools view -@ ${THREADS_PER_JOB} -bS ${ALIGN_DIR}/${sample}.sam | \
        samtools sort -@ ${THREADS_PER_JOB} -o ${ALIGN_DIR}/${sample}.sorted.bam
        
        samtools index ${ALIGN_DIR}/${sample}.sorted.bam
        
        rm ${ALIGN_DIR}/${sample}.sam
        
        echo "[${sample}] BAM 文件处理完成."
    else
        echo "[${sample}] 跳过步骤 3: Samtools"
    fi
    echo ">>> [PID: $$] 样本 ${sample} 处理完毕 <<<"
}


# --- 5. 主流程循环 (并行化) ---

echo "--- 开始并行处理所有样本 (最多 ${MAX_JOBS} 个作业同时进行) ---"

# 设置一个计数器来跟踪正在运行的作业
job_count=0

for fq1 in ${RAW_DIR}/*R1.raw.fastq.gz; do
    # 如果正在运行的作业达到了最大值，就等待任何一个作业完成
    if (( job_count >= MAX_JOBS )); then
        wait -n
        ((job_count--))
    fi

    # 启动一个新作业
    (
        fq2=${fq1/R1.raw.fastq.gz/R2.raw.fastq.gz}
        sample=$(basename ${fq1} .R1.raw.fastq.gz | cut -d'-' -f2)
        
        # 在后台执行处理函数
        process_sample "$fq1" "$fq2" "$sample"
    ) &
    
    ((job_count++))
done

# 等待所有剩余的后台作业完成
echo "--- 等待所有样本处理作业完成... ---"
wait
echo "--- 所有样本均已处理完毕. ---"
echo ""


# --- 6. featureCounts 基因计数 ---
if [[ ! ",${SKIP_STEPS}," =~ ",4," ]]; then
    echo "================================================="
    echo ">>> 所有样本比对完成，开始进行基因计数..."
    echo "================================================="

    # 获取所有排序后的 BAM 文件列表
    BAM_FILES=$(ls ${ALIGN_DIR}/*.sorted.bam)

    # Determine the number of available cores in a portable way
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        N_CORES=$(nproc)
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        N_CORES=$(sysctl -n hw.ncpu)
    else
        # Default to a safe number if OS is not recognized
        echo "Warning: Could not determine OS type. Defaulting to 4 threads for featureCounts."
        N_CORES=4
    fi

    # featureCounts can utilize all available cores
    featureCounts -T ${N_CORES} -p -t exon -g gene_id \
                  -a ${GTF_FILE} \
                  -o ${COUNT_DIR}/gene_counts.txt \
                  ${BAM_FILES}
else
    echo "--- 跳过步骤 4: featureCounts ---"
fi

echo "--- 分析流程全部完成！---"
echo "最终的基因计数矩阵保存在: ${COUNT_DIR}/gene_counts.txt"
echo "您可以在 gene_counts.txt 文件中查看每个基因在所有样本中的 read 计数值。"}
