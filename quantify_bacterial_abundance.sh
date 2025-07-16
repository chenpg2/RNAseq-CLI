#!/bin/bash

# =================================================================
# 目标细菌序列丰度计算脚本
#
# 功能:
# 1. 在指定的未比对上的 FASTQ 文件目录中，计算目标细菌序列的丰度。
#
# 作者: Gemini
# 日期: 2025-07-16
# =================================================================

# --- 1. 帮助信息和参数 ---

usage() {
    echo "Usage: $0 -i <unmapped_reads_dir> -f <fasta_file>"
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
        f) FASTA_FILE=$OPTARG ;;
        h) usage ;;
        \?) echo "无效选项: -$OPTARG" >&2; usage ;;
        :) echo "选项 -$OPTARG 需要一个参数." >&2; usage ;;
    esac
done

# 检查所有必需的参数是否已提供
if [ -z "${UNMAPPED_DIR}" ] || [ -z "${FASTA_FILE}" ]; then
    echo "错误: 缺少必需的参数。"
    usage
fi

# 检查 FASTA 文件是否存在
if [ ! -f "${FASTA_FILE}" ]; then
    echo "错误: FASTA 文件不存在: ${FASTA_FILE}"
    exit 1
fi

# 从 FASTA 文件中提取序列 (忽略标题行和换行符)
TARGET_SEQUENCE=$(grep -v "^>" "${FASTA_FILE}" | tr -d '\n')

if [ -z "${TARGET_SEQUENCE}" ]; then
    echo "错误: 未能从 FASTA 文件中读取序列。"
    exit 1
fi

# --- 3. 主流程 ---

echo "--- 开始计算目标序列在每个样本中的丰度 ---"

# 创建一个报告文件
REPORT_FILE="bacterial_abundance_report.txt"
echo -e "Sample\tTarget_Count\tTotal_Reads\tAbundance" > ${REPORT_FILE}

for fq in ${UNMAPPED_DIR}/*.unmapped.fastq.gz; do
    sample=$(basename ${fq} .unmapped.fastq.gz)
    
    echo "--- 正在处理样本: ${sample} ---"

    # 使用 zgrep 高效地在压缩文件中搜索
    target_count=$(zgrep -c "${TARGET_SEQUENCE}" ${fq})

    # FASTQ 文件每 4 行代表一个 read
    total_reads=$(($(zcat ${fq} | wc -l) / 4))

    # 计算丰度 (避免除以零)
    if [ ${total_reads} -gt 0 ]; then
        abundance=$(echo "scale=10; ${target_count} / ${total_reads}" | bc)
    else
        abundance=0
    fi

    echo "样本 ${sample}:"
    echo "  - 找到目标序列 ${target_count} 次"
    echo "  - 总 reads 数: ${total_reads}"
    echo "  - 丰度: ${abundance}"
    echo ""

    # 将结果写入报告
    echo -e "${sample}\t${target_count}\t${total_reads}\t${abundance}" >> ${REPORT_FILE}
done


echo "--- 丰度计算完成 --- "
echo "详细报告已保存到: ${REPORT_FILE}"
