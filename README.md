# HumanBeGone-ONT 🧬

**HumanBeGone-ONT** is a specialized variant of the HumanBeGone bioinformatics pipeline tailored specifically for Oxford Nanopore Technologies (ONT) long-read sequencing data. It rapidly decontaminates long-read FASTQ files by combining initial QC (via **Fastplong**) with robust k-mer mapping (**Kraken2**) and highly sensitive long-read aligners (**Minimap2** or **Winnowmap**) to scrub human sequences from your dataset securely.

## 🚀 Features
- **Long-Read QC**: Quality filtering uniquely designed for noisy Nanopore sequences via Fastplong.
- **Dual Alignment Architecture**: Flexible target modes. Use `--fast` for rapid decontamination (Minimap2) or `--high-accuracy` for rigorous, repetitive-element-aware alignment mapping (Winnowmap).
- **Directory & File Seamless Processing**: SampleSheets can point directly to single `.fastq.gz` pipelines or to raw structured folders outputted by sequencing hardware (natively utilizing `fastcat` under the hood).
- **Dynamic Interactive Reporting**: Automatically generates an interactive Chart.js HTML report mapping Nanopore survival across the workflow.

---

## 🛠️ Installation & Initialization 

The environment requires initialization and massive index databases to operate. We provide a fully automated `init.sh` script to streamline this!

1. Clone this repository locally.
2. Run the initialization script. This operates the integrated `humanbegone-ont.yml` conda environment, dynamically hooks into it, pulls down the massive Kraken2 databases, reference FASTA formats, and Winnowmap K-mers securely from Zenodo, and finally seamlessly generates your `.mmi` indexes automatically:
   ```bash
   ./init.sh
   ```
3. Always activate the physical conda environment when attempting to process data:
   ```bash
   conda activate humanbegone-ont
   ```

---

## 💻 Usage

The pipeline processes native sequences mapped via a `.csv` target metadata matrix.

**Basic Invocation**
```bash
./humanbegone-ont.sh <samplesheet.csv> --fast|--high-accuracy [OPTIONS]
```

### 1. SampleSheet Matrix Format
Your input target must strictly be a comma-separated `.csv` specifying exactly a `Sample_Name` mapping mapping to the location of the reads:
```csv
Sample_Name,ONT_Reads_Path
Sample01,path/to/sample01.fastq.gz
Sample02,path/to/barcode02_folder/
```
*Note: Your `ONT_Reads_Path` can be a naked `.fastq` file, a `.fastq.gz` file, or a full directory. If you provide a directory, the pipeline will systematically deploy `fastcat` behind the scenes to concatenate all valid files dynamically.*

### 2. Help Commands & Options Array

```bash
./humanbegone-ont.sh --help
```

**Mandatory Flags:**
* `--fast` or `--high-accuracy` : Dictates the baseline alignment logic. `--fast` utilizes Minimap2. `--high-accuracy` deploys Winnowmap.
* `--kraken-db <path>` : Targeted path specifying the extracted Kraken2 database directory.
* `--mm2-index <path>` : (*Required if using `--fast`*) Path to the generated Minimap2 index.
* `--wn-kmer <path>` : (*Required if using `--high-accuracy`*) Path to the repetitive k-mer table.
* `--ref-fasta <path>` : (*Required if using `--high-accuracy`*) Path to the reference human FASTA sequence.

**Optional Modifiers:**
* `--skip-fastplong` : Physically bypass the entire Fastplong QC module. Sequences will pipe completely raw natively into Kraken2. 
* `--threads <int>` : Parallel multithreading limits. *(Default: 8)*.
* `--output-dir <path>` : Custom destination directory mapping. All executed scratch data builds into a hidden workspace relative to your SampleSheet location and transports permanently to this explicit output endpoint upon conclusion. *(Default: `Results/` folder alongside your SampleSheet)*.

### 🌟 Example Executions

**Running in FAST Mode (Minimap2):**
```bash
./humanbegone-ont.sh ont_samples.csv --fast \
    --threads 24 \
    --kraken-db ./kraken_T2T_db/ \
    --mm2-index ./T2T_Human.mmi \
    --output-dir /storage/ont_cleaned/
```

**Running in HIGH ACCURACY Mode (Winnowmap):**
```bash
./humanbegone-ont.sh ont_samples.csv --high-accuracy \
    --threads 24 \
    --kraken-db ./kraken_T2T_db/ \
    --wn-kmer ./T2T_repetitive_kmers.txt \
    --ref-fasta ./T2T_Unique.fna \
    --output-dir /storage/ont_cleaned/
```

---

## 📊 Pipeline Outputs

Once seamlessly processed, your target arrays will neatly deploy hierarchically into organized endpoint metrics:
* `Results/Fastp/<SampleName>/` *(Quality filtering JSON logs, HTML metrics, and trimmed read `.gz` arrays)*
* `Results/Kraken/<SampleName>/` *(Terminal sequence tracking logs predicting scrub thresholds)*
* `Results/Aligner/<SampleName>/` *(The finalized, structurally stripped `.fastq.gz` clean payload alongside alignment `.log` mapping info)*
* `Results/Reports/` *(The final master accumulated `Summary.csv` mapping structure and your dynamically togglable `processing_report.html` interactive visual chart!)*
