#!/bin/bash

# =================================================================
# 目标细菌序列丰度计算脚本 (BWA 对齐版)
#
# 功能:
# 1. 使用 BWA-MEM 将未比对的 reads 与目标 FASTA 序列进行对齐。
# 2. 计算成功比对到目标序列的 reads 数量。
# 3. 输出每个样本中目标序列的相对丰度。
#
# 作者: Cpg
# 日期: 2025-07-16
# =================================================================

# --- 1. 帮助信息和参数 ---

usage() {
    echo "Usage: $0 -i <unmapped_reads_dir> -f <target_fasta>"
    echo ""
    echo "Options:"
    echo "  -i    [必须] 包含未比对上的 FASTQ 文件的目录 (e.g., /path/to/unmapped_reads)"
    echo "  -f    [必须] 包含目标细菌序列的 FASTA 文件"
    echo "  -h    显示此帮助信息"
    echo ""
    exit 1
}

# --- 2. 解析命令行参数 ---

while getopts ":i:f:h" opt; do
    case ${opt} in
        i) UNMAPPED_DIR=$OPTARG ;;
        f) TARGET_FASTA=$OPTARG ;;
        h) usage ;;
        \?) echo "无效选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数." >&2; usage ;;
    esac
done

# 检查所有必需的参数是否已提供
if [ -z "${UNMAPPED_DIR}" ] || [ -z "${TARGET_FASTA}" ]; then
    echo "错误: 缺少必需的参数。"
    usage
fi

# 检查依赖工具
for tool in bwa samtools; do
    if ! command -v $tool &> /dev/null; then
        echo "错误: 依赖工具 '${tool}' 未安装或不在 PATH 中。"
        echo "请先通过 Conda 环境 (environment.yml) 安装。"
        exit 1
    fi
done

# --- 3. BWA 索引 ---

# 创建一个临时目录来存放索引
INDEX_DIR=$(mktemp -d)
# 确保脚本退出时删除临时目录
trap 'rm -rf -- "$INDEX_DIR"' EXIT

BWA_INDEX="${INDEX_DIR}/target_sequence"

echo "--- 正在为目标序列创建 BWA 索引... ---"
bwa index -p "${BWA_INDEX}" "${TARGET_FASTA}"
echo "索引创建完成。"


# --- 4. 主流程 ---

echo "--- 开始计算目标序列在每个样本中的丰度 (使用 BWA 对齐) ---"

REPORT_FILE="bacterial_abundance_report_bwa.txt"
echo -e "Sample\tMapped_Reads\tTotal_Reads\tAbundance_Percentage" > ${REPORT_FILE}

for fq1 in ${UNMAPPED_DIR}/*.1.unmapped.fastq.gz; do
    sample=$(basename ${fq1} .1.unmapped.fastq.gz)
    fq2=${fq1/.1.unmapped.fastq.gz/.2.unmapped.fastq.gz}
    
    echo "--- 正在处理样本: ${sample} ---"

    # 计算总 reads 数
    total_reads=$(($(zcat ${fq1} | wc -l) / 4))

    # 使用 BWA-MEM 对齐并使用 Samtools 计数映射上的 reads
    # -F 4 标志用于排除未映射的 reads
    mapped_reads=$(bwa mem -t 4 "${BWA_INDEX}" "${fq1}" "${fq2}" 2>/dev/null | samtools view -c -F 4 -)

    # 计算丰度 (避免除以零)
    if [ ${total_reads} -gt 0 ]; then
        abundance=$(echo "scale=10; (${mapped_reads} / ${total_reads}) * 100" | bc)
    else
        abundance=0
    fi

    echo "样本 ${sample}:"
    echo "  - 成功比对到目标的 Reads: ${mapped_reads}"
    echo "  - 总 Reads 数: ${total_reads}"
    echo "  - 丰度 (%): ${abundance}"
    echo ""

    # 将结果写入报告
    echo -e "${sample}\t${mapped_reads}\t${total_reads}\t${abundance}" >> ${REPORT_FILE}
done

echo "--- 丰度计算完成 --- "
echo "详细报告已保存到: ${REPORT_FILE}"
