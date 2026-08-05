# LongAllele pipeline configuration
# Usage: bash longallele.sh config_template.sh
# Copy this file, fill in your paths and settings, then run.

# ──── Required ────────────────────────────────────────────────────────────────
SCOTCH_TARGET="/path/to/scotch_output"   # SCOTCH output directory
BAM_PATH="/data/bam/nexus-rna-002-tumor-long-read_minimap2_mdtagged_sorted.bam"
REF_FASTA="/data/fasta/GRCh38.p14.genome.chr17.fa.gz"
OUTPUT_DIR="/sample/"

# ──── Parallelization ─────────────────────────────────────────────────────────
N_GENE_JOBS=50   # SLURM array size for steps 1–3 (number of gene-parallel tasks)
N_SAMPLES=1      # number of BAM files; set >1 and use space-separated lists above
                 # for multi-sample analysis

# ──── Optional inputs ─────────────────────────────────────────────────────────
CELL_TYPE_DF=""  # path to CSV with Cell/CellType columns (leave empty if not needed)
PREFIX=""        # output filename prefix (leave empty for none)

# ──── Variant calling (steps 1–2) ────────────────────────────────────────────
DEPTH=3         # minimum read depth at SNV position
N_ALT_COUNT=3   # minimum alt-allele read count
MIN_MAPQ=20      # minimum mapping quality
MIN_BASEQ=5      # minimum base quality

# ──── EM haplotyping (step 3) ─────────────────────────────────────────────────
SEED=42                    # random seed
MAX_ITER=50                # maximum EM iterations per gene
TOL=1e-3                   # convergence tolerance
HET_FILTER=0.99            # heterozygosity probability threshold
CLF_INIT=true              # use SNV classifier scores to initialize EM (recommended)
# RNA_EDITING_DB=""        # override bundled hg38 editing DB (leave commented for hg38)
# SNV_CLASSIFIER=""        # path to trained SNV classifier .joblib (optional)
# HIGH_ARTIFACT_MODE=false # enable nascent-RNA leak filters for snRNA-seq (default off)

# ──── Downstream analysis (step 5) ───────────────────────────────────────────
EVENT_MODE="all_events"    # all_events | switching_events | fdr_events
SNV_EVENT_DISTANCE=50      # ±bp exonic distance for SNV–event linking
EVENT_MIN_READS=10         # minimum weighted reads per event test
N_WORKERS=4                # parallel workers for step 5
# ASTU_SIG_ONLY=false      # restrict step 5 to ASTU-significant genes (default off)

# ──── SLURM resources ─────────────────────────────────────────────────────────
PARTITION="cpu"   # shared across all steps

# Steps 1, 2, 1.5: lightweight per-gene / per-sample pileup
MEM_12="16G"  ; TIME_12="4:00:00"  ; CPUS_12=2
# Step 3: EM haplotyping — more memory for large genes
MEM_3="32G"   ; TIME_3="8:00:00"   ; CPUS_3=4
# Step 4: aggregation across all genes — single job, high memory
MEM_4="64G"   ; TIME_4="4:00:00"   ; CPUS_4=8
# Step 5: downstream analysis — parallel workers, high memory
MEM_5="64G"   ; TIME_5="6:00:00"   ; CPUS_5=8