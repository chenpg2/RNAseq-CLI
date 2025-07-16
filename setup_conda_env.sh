#!/bin/bash

# =================================================================
# Conda 环境一键部署脚本
#
# 功能:
# 1. 检查 Conda 是否已安装。
# 2. 使用 environment.yml 文件创建或更新 Conda 环境。
#
# 作者: Cpg
# 日期: 2025-07-16
# =================================================================

# --- 1. 检查 Conda 是否安装 ---
if ! command -v conda &> /dev/null; then
    echo "错误: Conda 未安装或未在您的 PATH 中。"
    echo "请先安装 Miniconda 或 Anaconda。"
    echo "访问: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html"
    exit 1
fi

# --- 2. 从 YML 文件创建环境 ---
ENV_FILE="environment.yml"

if [ ! -f "${ENV_FILE}" ]; then
    echo "错误: 未找到 ${ENV_FILE} 文件。"
    echo "请确保此脚本与 environment.yml 文件在同一目录下。"
    exit 1
fi

# 从 YML 文件中获取环境名称
ENV_NAME=$(grep 'name:' ${ENV_FILE} | cut -d' ' -f2)

echo "--- 准备创建或更新 Conda 环境: ${ENV_NAME} ---"

# 检查环境是否已存在
if conda env list | grep -q "^${ENV_NAME}"; then
    echo "环境 '${ENV_NAME}' 已存在。将尝试更新它。"
    conda env update --file ${ENV_FILE} --prune
else
    echo "正在创建新环境 '${ENV_NAME}'..."
    conda env create --file ${ENV_FILE}
fi

# --- 3. 完成 ---
if [ $? -eq 0 ]; then
    echo ""
    echo "--- Conda 环境 '${ENV_NAME}' 已成功设置！---"
    echo ""
    echo "要激活此环境, 请运行:"
    echo "    conda activate ${ENV_NAME}"
    echo ""
    echo "要停用此环境, 请运行:"
    echo "    conda deactivate"
    echo ""
else
    echo ""
    echo "--- 环境设置失败。请检查上面的错误信息。 ---"
    echo ""
fi
