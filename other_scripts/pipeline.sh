#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

ORTHOSIM_ROOT="/path/to/OrthoTrain_docs"
OUTPUT_BASE="/path/to/pipeline_output"

REAL_PROTEOMES="/path/to/real/proteomes"
REAL_SPECIES_TREE="/path/to/real/species_tree.nwk"
REAL_OUTPUT_OF="/path/to/real/orthofinder_output"
REAL_OUTPUT_FASTOMA="/path/to/real/fastoma_output"
REAL_OUTPUT_BROCCOLI="/path/to/real/broccoli_output"
REAL_OUTPUT_SONIC="/path/to/real/sonicparanoid_output"

PFAM_PROTEOME="${ORTHOSIM_ROOT}/PFAM/pfam_genomes/Saccharomyces_cerevisiae.fa"

N_ORTHOGROUPS=5000
THREADS=64

GA_EXTRA_FLAGS=""
# GA_EXTRA_FLAGS="--generations 100 --sample-size 50 --num-parents-mating 10"

# =============================================================================

ORTHOSIM="python ${ORTHOSIM_ROOT}/OrthoSim.py"
BENCHMARK="python ${ORTHOSIM_ROOT}/benchmark.py"
BENCHMARK_ORTHOLOGS="python ${ORTHOSIM_ROOT}/benchmark_orthologs.py"
RUN_TOOLS="python ${ORTHOSIM_ROOT}/run_tools.py"

TOOLS=(OF3 FASTOMA BROCOLI SONICPARANOID2)

declare -A BENCH_TOOL_NAME=(
    [OF3]=orthofinder [FASTOMA]=fastoma
    [BROCOLI]=broccoli [SONICPARANOID2]=sonicparanoid2
)
declare -A REAL_TOOL_OUTPUT=(
    [OF3]="$REAL_OUTPUT_OF"       [FASTOMA]="$REAL_OUTPUT_FASTOMA"
    [BROCOLI]="$REAL_OUTPUT_BROCCOLI" [SONICPARANOID2]="$REAL_OUTPUT_SONIC"
)

inferred_tool_output() {
    case "$2" in
        OF3)            echo "${1}/proteome_files/OrthoFinder/Results_I1.5" ;;
        FASTOMA)        echo "${1}/fastoma_output" ;;
        BROCOLI)        echo "${1}" ;;
        SONICPARANOID2) echo "${1}/sp_default" ;;
    esac
}

log() { echo; echo "=== $* ==="; echo; }

mkdir -p "$OUTPUT_BASE"

# --- Step 3: Get-Parameters and GA ---
for TOOL in "${TOOLS[@]}"; do
    log "Get-Parameters: ${TOOL}"
    GP_CMD=(
        $ORTHOSIM --Get-Parameters
        --output      "${OUTPUT_BASE}/training_${TOOL}"
        --tool        "$TOOL"
        --tool-output "${REAL_TOOL_OUTPUT[$TOOL]}"
        --threads     "$THREADS"
    )
    [[ "$TOOL" != "OF3" ]] && GP_CMD+=(--tools-proteomes "$REAL_PROTEOMES")
    [[ "$TOOL" == "BROCOLI" || "$TOOL" == "SONICPARANOID2" ]] && \
        GP_CMD+=(--tool-tree "$REAL_SPECIES_TREE")
    "${GP_CMD[@]}"

    log "GA: ${TOOL}"
    # shellcheck disable=SC2086
    $ORTHOSIM --GA \
        --output                   "${OUTPUT_BASE}/ga_${TOOL}" \
        --Orthogroup-train-results "${OUTPUT_BASE}/training_${TOOL}" \
        --Orthogroups              "$N_ORTHOGROUPS" \
        --threads                  "$THREADS" \
        $GA_EXTRA_FLAGS
done

# --- Step 4: Simulate ---
for TOOL in "${TOOLS[@]}"; do
    log "Simulate: ${TOOL}"
    $ORTHOSIM --Orthogroup-simulation \
        --Config      "${OUTPUT_BASE}/ga_${TOOL}/config.txt" \
        --output      "${OUTPUT_BASE}/sim_${TOOL}" \
        --Orthogroups "$N_ORTHOGROUPS" \
        --threads     "$THREADS"
done

# --- Step 5: Run all tools on each simulation ---
for TRAIN_TOOL in "${TOOLS[@]}"; do
    log "Run tools: sim_${TRAIN_TOOL}"
    $RUN_TOOLS \
        --sim-folder "${OUTPUT_BASE}/sim_${TRAIN_TOOL}" \
        --threads    "$THREADS"
done

# --- Step 6: Benchmark (4x4 matrix) ---
BENCH_DIR="${OUTPUT_BASE}/benchmarks"
mkdir -p "$BENCH_DIR"

for TRAIN_TOOL in "${TOOLS[@]}"; do
    SIM_DIR="${OUTPUT_BASE}/sim_${TRAIN_TOOL}"
    for INFER_TOOL in "${TOOLS[@]}"; do
        LABEL="${TRAIN_TOOL}__${INFER_TOOL}"
        TOOL_OUT=$(inferred_tool_output "$SIM_DIR" "$INFER_TOOL")
        BENCH="${BENCH_TOOL_NAME[$INFER_TOOL]}"

        $BENCHMARK \
            --output-folder "$SIM_DIR" --tool-output "$TOOL_OUT" \
            --tool "$BENCH" --csv "${BENCH_DIR}/og_${LABEL}.csv"

        $BENCHMARK_ORTHOLOGS \
            --output-folder "$SIM_DIR" --tool-output "$TOOL_OUT" \
            --tool "$BENCH" --csv "${BENCH_DIR}/orthologs_${LABEL}.csv"
    done
done

log "Done. Results in: ${BENCH_DIR}"
