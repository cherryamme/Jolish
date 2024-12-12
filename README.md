# Jolish
![jolish](assets/Jolish-logo2.png)
Jolish is a powerful tool designed to enhance the quality of third-generation sequencing FASTQ data. By employing advanced algorithms, Jolish accurately corrects errors and optimizes read quality, ensuring more reliable downstream analysis.

## Features

- **Error Correction**: Accurately corrects errors in sequencing data.
- **Quality Optimization**: Enhances read quality for better analysis.
- **High Performance**: Efficient processing with multi-threading support.

## Installation

Clone the repository and build with Cargo:

```bash
git clone https://github.com/yourusername/Jolish.git
cd Jolish
cargo build --release
```
## Example

Correct file

```bash
# correct fastq file
Jolish -i example_reads.fq -o corrected_reads.fq
# or correct bam file
Jolish -i example_reads.fq.bam -o corrected_reads.fq
```


## Usage

After building, run Jolish using the command line:

```bash
./target/release/Jolish -i input.fq -o output.fq
```

- `-i`, `--inputs <INPUTS>...`：输入文件路径
- `-o`, `--outfile <OUTFILE>`：输出文件名（默认：outfile.fq.gz）
- `-M`, `--max-reads <MAX_READS>`：最大处理的读取数
- `-t`, `--threads <THREADS>`：线程数（默认：10）
- `-m`, `--min-length <MIN_LENGTH>`：按最小长度过滤读取（默认：100）
- `-w`, `--chunk-window <CHUNK_WINDOW>`：窗口大小（默认：100）
- `-c`, `--chunk-size <CHUNK_SIZE>`：块大小（默认：100）
- `--cr <CORRECT_RATIO>`：校正比例（默认：0.3）
- `--cfr <CORRECT_FIX_RATIO>`：校正修复比例（默认：0.3）
- `--id_sep <SPLIT_PARAM>`：ID 分隔符（默认：%）
- `-p`, `--type <TYPE_INDEX>`：引物类型索引（默认：1）
- `-s`, `--strand <ORIENT_INDEX>`：链方向索引（默认：3）
- `-r`, `--region <REGION>`：提取的区域（例如：chr1:1000-2000）
- `-h`, `--help`：打印帮助信息
- `-V`, `--version`：打印版本信息

## Dependencies

- Rust and Cargo for building the project.
- Optional: `samtools`, `minimap2` for downstream analysis.

## Contributing

Contributions are welcome. Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
```
