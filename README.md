# RNA-seq Analysis Pipeline

This repository contains a simple and efficient pipeline for processing paired-end RNA-seq data. The pipeline automates the steps from raw FASTQ files to a gene count matrix, leveraging parallel processing to significantly speed up the analysis of multiple samples.

## Features

- **Command-Line Interface**: Easy to use with command-line arguments to specify paths and resources.
- **Parallel Processing**: Processes multiple samples concurrently to maximize CPU usage and reduce overall runtime.
- **Automated Workflow**: Integrates quality control, read alignment, and gene counting into a single script.
- **Standard Tools**: Built using widely adopted bioinformatics tools:
    - [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for quality and adapter trimming.
    - [HISAT2](http://daehwankimlab.github.io/hisat2/) for fast and sensitive alignment.
    - [Samtools](http://www.htslib.org/) for processing alignment files.
    - [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) for accurate gene-level quantification.

## Prerequisites

Before running the pipeline, ensure you have [Conda](https://docs.conda.io/en/latest/miniconda.html) installed. All required bioinformatics tools will be installed automatically into a dedicated Conda environment.

## Environment Setup

To ensure all dependencies are handled correctly, we strongly recommend using Conda to create a dedicated environment for this pipeline.

1.  **Create the Environment**: The `environment.yml` file lists all the required software. The `setup_conda_env.sh` script automates the environment creation.

    Make the setup script executable and run it:
    ```bash
    chmod +x setup_conda_env.sh
    ./setup_conda_env.sh
    ```
    This will create a new Conda environment named `rnaseq_analysis_env`.

2.  **Activate the Environment**: Before running the main pipeline, activate the newly created environment:
    ```bash
    conda activate rnaseq_analysis_env
    ```

## Directory Structure

The pipeline expects the following directory structure:

```
/path/to/your/project/
├── rawdata/
│   ├── sample1_R1.raw.fastq.gz
│   ├── sample1_R2.raw.fastq.gz
│   ├── sample2_R1.raw.fastq.gz
│   └── sample2_R2.raw.fastq.gz
├── reference/
│   ├── annotation/
│   │   └── gencode.vM25.annotation.gtf
│   └── genome/
│       ├── GRCm38.p6.genome.fa
│       └── (hisat2_index_files)
└── run_rnaseq_pipeline_cli.sh
```

- `rawdata/`: Contains the raw paired-end FASTQ files.
- `reference/`: Contains the reference genome (`.fa`), HISAT2 index files, and gene annotation (`.gtf`).

## Usage

The main script `run_rnaseq_pipeline_cli.sh` is controlled via command-line arguments.

### Command-Line Options

```
Usage: ./run_rnaseq_pipeline_cli.sh -g <genome_index_prefix> -a <annotation_gtf> -i <input_dir> -o <output_dir> [-j <max_jobs>] [-n <threads_per_job>]

Options:
  -g    [必须] HISAT2 基因组索引的前缀路径 (e.g., /path/to/ref/mm10)
  -a    [必须] GTF 注释文件路径 (e.g., /path/to/annotation.gtf)
  -i    [必须] 包含原始 FASTQ 文件的输入目录
  -o    [必须] 用于存放所有结果的输出目录
  -j    [可选] 并行处理的样本作业数 (默认: 8)
  -n    [可选] 每个作业内部使用的线程数 (默认: 4)
  -h    显示此帮助信息
```

### Example

1.  **Make the script executable:**
    ```bash
    chmod +x run_rnaseq_pipeline_cli.sh
    ```

2.  **Run the pipeline:**
    ```bash
    ./run_rnaseq_pipeline_cli.sh \
        -g /path/to/your/project/reference/genome/mm10 \
        -a /path/to/your/project/reference/annotation/gencode.vM25.annotation.gtf \
        -i /path/to/your/project/rawdata \
        -o /path/to/your/project/output_results \
        -j 12 \
        -n 8
    ```

## Output

The pipeline will create the specified output directory (`-o`) with the following subdirectories:

-   `cleandata/`: Contains trimmed FASTQ files and FastQC reports.
-   `alignment/`: Contains sorted and indexed BAM files for each sample.
-   `counts/`: Contains the final gene count matrix named `gene_counts.txt`. This file can be directly used for downstream differential gene expression analysis.


## License

This project is licensed under the MIT License.

```