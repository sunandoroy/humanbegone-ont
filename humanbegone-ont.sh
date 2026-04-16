#!/bin/bash

START_TIME=$(date +%s)

# Usage: 
#   ./humanbegone-ont.sh samplesheet.csv --fast|--high-accuracy [OPTIONS]

INPUT=""
MODE=""
SKIP_FASTP=false
THREADS=12
K2_DB=""
MM2_INDEX=""
WN_KMER_TABLE=""
REF_FASTA=""
OUTPUT_DIR=""

print_help() {
    echo "Usage: $0 <samplesheet_csv> [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help                Show this help menu"
    echo "  --fast | --high-accuracy  (required) Set ONT processing mode (minimap2 vs winnowmap)"
    echo "  --kraken-db <path>        (required) Path to Kraken2 database directory"
    echo "  --mm2-index <path>        (required for --fast) Path to Minimap2 index"
    echo "  --wn-kmer <path>          (required for --high) Path to repetitive k-mer table"
    echo "  --ref-fasta <path>        (required for --high) Path to Reference Fasta"
    echo "  --output-dir <path>       Custom directory to output Results folder"
    echo "  --skip-fastplong          Skip FastpLong processing"
    echo "  --threads <int>           Number of threads (default: 12)"
    exit 0
}

# Parse arguments cleanly allowing for dynamic flags and key-value pairs
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) print_help ;;
        --fast) MODE="--fast"; shift ;;
        --high-accuracy) MODE="--high-accuracy"; shift ;;
        --threads) THREADS="$2"; shift 2 ;;
        --kraken-db) K2_DB="$2"; shift 2 ;;
        --mm2-index) MM2_INDEX="$2"; shift 2 ;;
        --wn-kmer) WN_KMER_TABLE="$2"; shift 2 ;;
        --ref-fasta) REF_FASTA="$2"; shift 2 ;;
        --skip-fastplong) SKIP_FASTP=true; shift ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        *)
            if [ -z "$INPUT" ]; then
                INPUT="$1"
            else
                echo "Unknown argument/input: $1"
                echo "Run with --help for more information."
                exit 1
            fi
            shift
            ;;
    esac
done

if [ -z "$INPUT" ] || [ -z "$MODE" ] || [ -z "$K2_DB" ]; then
    echo "Error: Missing required basic arguments (Input, Mode, Kraken DB)."
    echo "Run with --help for more information."
    exit 1
fi

if [[ "$MODE" == "--fast" && -z "$MM2_INDEX" ]]; then
    echo "Error: Missing --mm2-index for fast mode."
    exit 1
elif [[ "$MODE" == "--high-accuracy" && ( -z "$WN_KMER_TABLE" || -z "$REF_FASTA" ) ]]; then
    echo "Error: Missing --wn-kmer and/or --ref-fasta for high-accuracy mode."
    exit 1
fi

if [[ ! -f "$INPUT" || "$INPUT" != *.csv ]]; then
    echo "Error: Input must be a valid SampleSheet CSV file."
    exit 1
fi

# Preserve original directory
ORIGINAL_DIR="$PWD"

# Convert dynamic relative paths into static absolute paths
if [[ "$INPUT" != /* ]]; then INPUT="$ORIGINAL_DIR/$INPUT"; fi
if [[ -n "$K2_DB" && "$K2_DB" != /* ]]; then K2_DB="$ORIGINAL_DIR/$K2_DB"; fi
if [[ -n "$MM2_INDEX" && "$MM2_INDEX" != /* ]]; then MM2_INDEX="$ORIGINAL_DIR/$MM2_INDEX"; fi
if [[ -n "$WN_KMER_TABLE" && "$WN_KMER_TABLE" != /* ]]; then WN_KMER_TABLE="$ORIGINAL_DIR/$WN_KMER_TABLE"; fi
if [[ -n "$REF_FASTA" && "$REF_FASTA" != /* ]]; then REF_FASTA="$ORIGINAL_DIR/$REF_FASTA"; fi
if [[ -n "$OUTPUT_DIR" && "$OUTPUT_DIR" != /* ]]; then OUTPUT_DIR="$ORIGINAL_DIR/$OUTPUT_DIR"; fi

# Extract working directory dynamically mapped to SampleSheet origin
BASE_INPUT_DIR="$(dirname "$INPUT")"
WORK_DIR="${BASE_INPUT_DIR}/work"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR" || exit 1

# Initialize global Summary CSV array
SUMMARY_FILE="scrubbing_summary_$(date +%Y%m%d_%H%M%S).csv"
echo "Sample_Name,Initial_Reads,After_fastp,Kraken_Human_Removed,Aligner_Human_Removed,Final_Clean_Reads,Percent_Scrubbed" > "$SUMMARY_FILE"

echo "==========================================="
echo "   ONT-SCRUB"
echo "   Mode: $MODE"
echo "==========================================="

# --- PRE-FLIGHT CHECK ---
TOOLS="kraken2 samtools awk"
if [ "$SKIP_FASTP" = false ]; then TOOLS="$TOOLS fastplong gojq"; fi
if [[ "$MODE" == "--fast" ]]; then TOOLS="$TOOLS minimap2"; else TOOLS="$TOOLS winnowmap"; fi

for tool in $TOOLS; do
	if ! command -v "$tool" &> /dev/null; then
		echo "ERROR: $tool is not installed or not in PATH."
		exit 1
	fi
done
# Checking fastcat just in case they drop folders in mapping
if ! command -v "fastcat" &> /dev/null; then
    echo "WARNING: fastcat is not installed. Will not be able to process directories in SampleSheet natively."
fi

# --- CORE SCRUBBING FUNCTION ---
process_sample() {
	local NAME=$1
	local READS_PATH=$2

    # Resolve relative paths relative to the Input file's directory ($BASE_INPUT_DIR)
    if [[ ! -z "$READS_PATH" && "$READS_PATH" != /* ]]; then
        READS_PATH="${BASE_INPUT_DIR}/$READS_PATH"
    fi

    echo "-------------------------------------------"
	echo "Processing Sample: $NAME"
    echo "Reads Path = '$READS_PATH'"
    
    if [ ! -e "$READS_PATH" ]; then
        echo "Error: Target reads path not found for $NAME at $READS_PATH"
        return 1
    fi

    # 1. Input Processing (fastcat if folder, direct bind if simple file)
    echo "[1/4] Assimilating Input..."
    if [ -d "$READS_PATH" ]; then
        fastcat "$READS_PATH" -s "$NAME" > "${NAME}_raw.fastq"
    else
        # Determine if it's compressed or naked
        if [[ "$READS_PATH" == *.gz ]]; then
            zcat -f "$READS_PATH" > "${NAME}_raw.fastq"
        else
            cat "$READS_PATH" > "${NAME}_raw.fastq"
        fi
    fi
    rm -r fastcat-histograms
    TOTAL_LINES=$(wc -l < "${NAME}_raw.fastq" | awk '{print $1}')
    INITIAL_READS=$((TOTAL_LINES / 4))

    if [ "$INITIAL_READS" -eq 0 ]; then
        echo "Error: No reads extracted for $NAME. Skipping."
        rm -f "${NAME}_raw.fastq"
        return 1
    fi

    if [ "$SKIP_FASTP" = false ]; then
        # 2. fastplong (QC)
        echo "[2/4] Running fastplong QC..."
        fastplong -i "${NAME}_raw.fastq" -o "${NAME}_qc.fastq.gz" -q 10 --json "${NAME}_fastp.json" --html "${NAME}_fastp.html" --thread $THREADS 2> "${NAME}_fastp.log"
        QC_READS=$(gojq '.summary.after_filtering.total_reads' "${NAME}_fastp.json")
        QC_READS=${QC_READS:-0}
        K2_INPUT="${NAME}_qc.fastq.gz"
    else
        echo "[2/4] Skipping fastplong QC..."
        QC_READS=$INITIAL_READS
        K2_INPUT="${NAME}_raw.fastq"
    fi

    # 3. Kraken 2 (K-mer Scrub)
    echo "[3/4] Running Kraken 2..."
    kraken2 --db "$K2_DB" --threads $THREADS "$K2_INPUT" --confidence 0.05 --unclassified-out "${NAME}_k2_clean.fastq" --report "${NAME}_k2_report.txt" > "${NAME}_k2.log" 2>&1
    K2_CLEAN_COUNT=$(grep " sequences unclassified" "${NAME}_k2.log" | awk '{print $1}')
    K2_CLEAN_COUNT=${K2_CLEAN_COUNT:-0}

    # 4. Alignment Scrub (Minimap2 vs Winnowmap2)
    echo "[4/4] Running Alignment Scrub ($MODE)..."
    if [[ "$MODE" == "--fast" ]]; then
    	minimap2 -ax map-ont -k13 -w5 --secondary=no -t $THREADS "$MM2_INDEX" "${NAME}_k2_clean.fastq" > "${NAME}_mapped.sam" 2> "${NAME}_aligner.log"
    else
    	winnowmap -ax map-ont -k13 -w5 -W "$WN_KMER_TABLE" --secondary=no --split-prefix=yes -t $THREADS "$REF_FASTA" "${NAME}_k2_clean.fastq" > "${NAME}_mapped.sam" 2> "${NAME}_aligner.log"
    fi

    samtools fastq -f 4 -G 2048 "${NAME}_mapped.sam" | gzip > "${NAME}_final_scrubbed.fastq.gz"
    FINAL_CLEAN_COUNT=$(zcat "${NAME}_final_scrubbed.fastq.gz" | wc -l | awk '{print $1}')
    FINAL_CLEAN_COUNT=$((FINAL_CLEAN_COUNT / 4))
    FINAL_CLEAN_COUNT=${FINAL_CLEAN_COUNT:-0}

    # Calculate Stats
    K2_REMOVED=$((QC_READS - K2_CLEAN_COUNT))
    ALIGN_REMOVED=$((K2_CLEAN_COUNT - FINAL_CLEAN_COUNT))
    TOTAL_REMOVED=$((INITIAL_READS - FINAL_CLEAN_COUNT))
    
    if [ "$INITIAL_READS" -gt 0 ]; then
        PCT_REMOVED=$(awk "BEGIN {printf \"%.2f\", ($TOTAL_REMOVED/$INITIAL_READS)*100}")
    else
        PCT_REMOVED=0
    fi

    # Append to Summary File
    echo "$NAME,$INITIAL_READS,$QC_READS,$K2_REMOVED,$ALIGN_REMOVED,$FINAL_CLEAN_COUNT,$PCT_REMOVED%" >> "$SUMMARY_FILE"

    # Cleanup intermediate sample files in active loop
    rm -f "${NAME}_raw.fastq" "${NAME}_mapped.sam" "${NAME}_k2_clean.fastq"
}

# --- INPUT ROUTING (SAMPLESHEET ONLY) ---
sed -i 's/\r//g' "$INPUT"
while IFS=',' read -r col1 col2 trash; do
    # Skip headers
    if [[ "$col1" =~ ^[Ss]ample ]]; then continue; fi
    # Execute sample route natively
    if [[ -n "$col1" && -n "$col2" ]]; then
        process_sample "$col1" "$col2"
    fi
done < "$INPUT"

# --- HTML REPORT GENERATION MODULE ---
echo -e "\nScrubbing complete. Final summary saved to: $SUMMARY_FILE"
column -s, -t < "$SUMMARY_FILE"

REPORT_FILE="scrubbing_report_$(date +%Y%m%d_%H%M%S).html"
echo -e "\nGenerating HTML report: $REPORT_FILE"

cat << 'EOF' > "$REPORT_FILE"
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HumanBeGone-ONT processing Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body { font-family: 'Inter', sans-serif; margin: 40px; background-color: #f8f9fa; color: #333; }
        .container { background-color: white; padding: 30px; border-radius: 12px; box-shadow: 0 10px 15px rgba(0,0,0,0.05); max-width: 1200px; margin: auto; }
        h1 { text-align: center; color: #2c3e50; font-weight: 600; margin-bottom: 20px; }
        .controls { text-align: center; margin-bottom: 20px; }
        button { background-color: #3498db; color: white; border: none; padding: 10px 20px; border-radius: 6px; cursor: pointer; font-size: 14px; font-weight: bold; transition: background-color 0.2s; }
        button:hover { background-color: #2980b9; }
        .chart-container { position: relative; height: 60vh; width: 100%; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Scrubbing Summary Report</h1>
        <div class="controls">
            <button id="toggleViewBtn">Switch to Percentage View</button>
        </div>
        <div class="chart-container">
            <canvas id="scrubChart"></canvas>
        </div>
    </div>
    <script>
        const data = [
EOF

# Parse the summary file to generate JSON-like data handling edge cases
tail -n +2 "$SUMMARY_FILE" | while IFS=',' read -r name initial after_fastp kraken aligner final pct; do
    fastp_removed=$((initial - after_fastp))
    fastp_removed=${fastp_removed:-0}
    kraken=${kraken:-0}
    aligner=${aligner:-0}
    final=${final:-0}
    initial=${initial:-0}
    
    if [ "$initial" -gt 0 ]; then
        p_fastp=$(awk -v f="$fastp_removed" -v i="$initial" 'BEGIN { printf "%.2f", (f/i)*100 }')
        p_kraken=$(awk -v k="$kraken" -v i="$initial" 'BEGIN { printf "%.2f", (k/i)*100 }')
        p_align=$(awk -v b="$aligner" -v i="$initial" 'BEGIN { printf "%.2f", (b/i)*100 }')
        p_final=$(awk -v f="$final" -v i="$initial" 'BEGIN { printf "%.2f", (f/i)*100 }')
    else
        p_fastp=0 p_kraken=0 p_align=0 p_final=0
    fi
    
    echo "            { label: '$name', count: { fastp: $fastp_removed, kraken: $kraken, align: $aligner, final: $final }, pct: { fastp: $p_fastp, kraken: $p_kraken, align: $p_align, final: $p_final } }," >> "$REPORT_FILE"
done

cat << 'EOF' >> "$REPORT_FILE"
        ];

        const labels = data.map(d => d.label);
        let isCountView = true;
        
        function getDatasets(viewType) {
            return [
                { 
                    label: 'Final Clean Reads', 
                    data: data.map(d => d[viewType].final), 
                    backgroundColor: 'rgba(75, 192, 192, 0.8)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 1
                },
                { 
                    label: 'Lost in Aligner', 
                    data: data.map(d => d[viewType].align), 
                    backgroundColor: 'rgba(255, 205, 86, 0.8)',
                    borderColor: 'rgba(255, 205, 86, 1)',
                    borderWidth: 1
                },
                { 
                    label: 'Lost in Kraken2', 
                    data: data.map(d => d[viewType].kraken), 
                    backgroundColor: 'rgba(255, 159, 64, 0.8)',
                    borderColor: 'rgba(255, 159, 64, 1)',
                    borderWidth: 1
                },
                { 
                    label: 'Lost in FastpLong', 
                    data: data.map(d => d[viewType].fastp), 
                    backgroundColor: 'rgba(255, 99, 132, 0.8)',
                    borderColor: 'rgba(255, 99, 132, 1)',
                    borderWidth: 1
                }
            ];
        }

        const ctx = document.getElementById('scrubChart').getContext('2d');
        const chart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: labels,
                datasets: getDatasets('count')
            },
            options: {
                indexAxis: 'y',
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    x: {
                        stacked: true,
                        title: { display: true, text: 'Number of Reads', font: { weight: 'bold' } },
                        grid: { color: 'rgba(0,0,0,0.05)' }
                    },
                    y: {
                        stacked: true,
                        title: { display: true, text: 'Samples', font: { weight: 'bold' } },
                        grid: { display: false }
                    }
                },
                plugins: {
                    title: { display: false },
                    tooltip: { mode: 'index', intersect: false },
                    legend: { position: 'top' }
                }
            }
        });
        
        document.getElementById('toggleViewBtn').addEventListener('click', function() {
            isCountView = !isCountView;
            this.textContent = isCountView ? "Switch to Percentage View" : "Switch to Count View";
            
            chart.data.datasets = getDatasets(isCountView ? 'count' : 'pct');
            chart.options.scales.x.title.text = isCountView ? 'Number of Reads' : 'Percentage of Reads (%)';
            if (!isCountView) {
                chart.options.scales.x.max = 100;
            } else {
                delete chart.options.scales.x.max;
            }
            chart.update();
        });
    </script>
</body>
</html>
EOF

echo "Done! You can view the report in your browser."

# --- OUTPUT ARCHIVING ---
echo "==========================================="
echo "   Organizing outputs into Results directory"
echo "==========================================="
mkdir -p Results/Reports
mv "$SUMMARY_FILE" "$REPORT_FILE" Results/Reports/ 2>/dev/null || true

# Parse summary file to organize by sample
if [ -f "Results/Reports/$SUMMARY_FILE" ]; then
    tail -n +2 "Results/Reports/$SUMMARY_FILE" | while IFS=',' read -r name trash; do
        mkdir -p "Results/Fastp/$name" "Results/Kraken/$name" "Results/Aligner/$name"
        mv "${name}_fastp"* "${name}_qc"*.fq.gz "Results/Fastp/$name/" 2>/dev/null || true
        mv "${name}_k2_report.txt" "${name}_k2.log" "Results/Kraken/$name/" 2>/dev/null || true
        mv "${name}_aligner.log" "${name}_final"*.fastq* "Results/Aligner/$name/" 2>/dev/null || true
    done
fi

if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${BASE_INPUT_DIR}/Results"
fi

ABS_LOCAL=$(readlink -f ./Results)
mkdir -p "$OUTPUT_DIR"
ABS_FINAL=$(readlink -f "$OUTPUT_DIR")

if [ "$ABS_LOCAL" != "$ABS_FINAL" ]; then
    echo "Relocating Results to output directory: $OUTPUT_DIR"
    mv Results/* "$OUTPUT_DIR/" 2>/dev/null || true
    rmdir Results 2>/dev/null || true
fi

# Clean Temp Env
cd "$ORIGINAL_DIR" || exit 1
rm -rf "$WORK_DIR"

echo "Cleanup complete! All files logically sorted in $OUTPUT_DIR/"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
RUNTIME_H=$((RUNTIME / 3600))
RUNTIME_M=$(((RUNTIME % 3600) / 60))
RUNTIME_S=$((RUNTIME % 60))
echo "Runtime : ${RUNTIME_H}h ${RUNTIME_M}m ${RUNTIME_S}s"
