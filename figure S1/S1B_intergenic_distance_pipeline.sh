#!/usr/bin/env bash
#
# run_intergenic_distance_pipeline.sh
#
# Compute intergenic read distance histograms stratified by gene class.
#
# Pipeline Overview:
#   1. Parse GTF to extract gene coordinates with biotype and 3'UTR annotations
#   2. Filter BAM files (remove unmapped, secondary, duplicates; MAPQ >= 30)
#   3. Extract 3' endpoints of aligned reads (strand-aware)
#   4. Identify strictly intergenic reads (no overlap with any gene body)
#   5. Compute signed distance to nearest gene 5' or 3' end
#   6. Bin distances into histograms stratified by gene class
#
# Gene Classes:
#   - protein_coding_3UTR:   protein-coding genes with annotated 3'UTR
#   - protein_coding_no3UTR: protein-coding genes without annotated 3'UTR
#   - lncRNA_3UTR:           lncRNA genes with annotated 3'UTR
#   - lncRNA_no3UTR:         lncRNA genes without annotated 3'UTR
#
# Output:
#   Per-sample histogram TSV with columns:
#     sample, gene_group, region, bin_start_bp, bin_end_bp, bin_mid_bp, count, total_good_reads
#
# Usage:
#   ./run_intergenic_distance_pipeline.sh GTF BAM_DIR OUT_DIR [BIN] [MIN] [MAX]
#
# Arguments:
#   GTF      - Ensembl GTF annotation file
#   BAM_DIR  - Directory containing input BAM files (*.bam)
#   OUT_DIR  - Output directory for results
#   BIN      - Histogram bin width in bp (default: 500)
#   MIN      - Minimum distance in bp (default: -10000)
#   MAX      - Maximum distance in bp (default: 10000)
#
# Requirements:
#   - samtools (>= 1.10)
#   - bedtools (>= 2.29)
#   - GNU awk, sort, paste

set -euo pipefail

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 GTF BAM_DIR OUT_DIR [BIN] [MIN] [MAX]" >&2
  exit 1
fi

GTF="$1"
BAM_DIR="$2"
OUT_DIR="$3"
BIN="${4:-500}"
MIN="${5:--10000}"
MAX="${6:-10000}"

mkdir -p "${OUT_DIR}"
mkdir -p "${OUT_DIR}/filtered_bam" "${OUT_DIR}/beds" "${OUT_DIR}/hist"

echo "[INFO] GTF     = ${GTF}"
echo "[INFO] BAM_DIR = ${BAM_DIR}"
echo "[INFO] OUT_DIR = ${OUT_DIR}"
echo "[INFO] BIN     = ${BIN} bp, MIN = ${MIN} bp, MAX = ${MAX} bp"

# -----------------------------------------------------------------------------
# Step 1: Build gene BED files from GTF
# -----------------------------------------------------------------------------
# Creates three BED files with gene coordinates tagged by biotype and 3'UTR status:
#   - genes.bed:    Full gene bodies
#   - genes_5p.bed: 5' endpoints (TSS)
#   - genes_3p.bed: 3' endpoints (TES)
#
# Name field format: gene_id|biotype|3UTR_status
#   biotype: protein_coding, lncRNA, or other
#   3UTR_status: 3UTR or no3UTR

GENE_BED="${OUT_DIR}/genes.bed"
GENE_5P_BED="${OUT_DIR}/genes_5p.bed"
GENE_3P_BED="${OUT_DIR}/genes_3p.bed"

echo "[INFO] Building gene BEDs with biotype and 3'UTR annotations"

# AWK script for GTF parsing (used for all three BED outputs)
read -r -d '' AWK_GTF_PARSER << 'AWKSCRIPT' || true
BEGIN { FS = OFS = "\t" }

function trim(s) {
  gsub(/^[ \t]+|[ \t]+$/, "", s)
  return s
}

function get_attr(attr_str, key,    n, i, f) {
  n = split(attr_str, f, ";")
  for (i = 1; i <= n; i++) {
    f[i] = trim(f[i])
    if (f[i] ~ ("^" key "[ \t]")) {
      sub(("^" key "[ \t]+"), "", f[i])
      gsub(/"/, "", f[i])
      return f[i]
    }
  }
  return ""
}

function norm_biotype(bt) {
  if (bt == "protein_coding") return "protein_coding"
  if (bt == "lncRNA" || bt == "lincRNA") return "lncRNA"
  return "other"
}

function label(gid,    bt, utr) {
  bt = (gid in biotype ? biotype[gid] : "other")
  utr = (gid in has3utr ? "3UTR" : "no3UTR")
  return gid "|" bt "|" utr
}

# Pass 1: Collect biotype and 3'UTR information
FNR == NR {
  if ($3 == "gene") {
    gid = get_attr($9, "gene_id")
    bt = get_attr($9, "gene_biotype")
    if (bt == "") bt = get_attr($9, "gene_type")
    bt = norm_biotype(bt)
    if (gid != "") biotype[gid] = bt
  } else if ($3 == "three_prime_utr") {
    gid = get_attr($9, "gene_id")
    if (gid != "") has3utr[gid] = 1
  }
  next
}

# Pass 2: Emit BED records
$3 == "gene" {
  gid = get_attr($9, "gene_id")
  if (gid == "") gid = "NA"
  nm = label(gid)
  
  if (mode == "GENE") {
    print $1, $4 - 1, $5, nm, 0, $7
  } else if (mode == "5P") {
    if ($7 == "+") { s = $4 - 1; e = $4 } else { s = $5 - 1; e = $5 }
    print $1, s, e, nm, 0, $7
  } else if (mode == "3P") {
    if ($7 == "+") { s = $5 - 1; e = $5 } else { s = $4 - 1; e = $4 }
    print $1, s, e, nm, 0, $7
  }
}
AWKSCRIPT

# Generate BED files
awk -v mode="GENE" "${AWK_GTF_PARSER}" "${GTF}" "${GTF}" > "${GENE_BED}"
awk -v mode="5P" "${AWK_GTF_PARSER}" "${GTF}" "${GTF}" > "${GENE_5P_BED}"
awk -v mode="3P" "${AWK_GTF_PARSER}" "${GTF}" "${GTF}" > "${GENE_3P_BED}"

# Sort BED files for bedtools -sorted operation
echo "[INFO] Sorting gene BEDs"
for f in "${GENE_BED}" "${GENE_5P_BED}" "${GENE_3P_BED}"; do
  LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n "${f}" -o "${f}"
done

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

filter_bam() {
  # Filter BAM: remove unmapped, secondary, supplementary, duplicates; require MAPQ >= 30
  local in_bam="$1"
  local out_bam="$2"
  local tmp_bam="${out_bam}.tmp"

  echo "    Filtering BAM -> ${out_bam##*/}"
  samtools view -F 0x904 -q 30 -b "${in_bam}" > "${tmp_bam}"

  if ! samtools quickcheck -q "${tmp_bam}"; then
    echo "    [ERROR] BAM validation failed for ${tmp_bam}" >&2
    rm -f "${tmp_bam}"
    exit 1
  fi

  mv "${tmp_bam}" "${out_bam}"
}

bam_to_3prime_bed_sorted() {
  # Convert BAM to BED of read 3' endpoints (strand-aware)
  local bam="$1"
  local out_bed="$2"

  bedtools bamtobed -i "${bam}" \
    | awk 'BEGIN { FS = OFS = "\t" }
      {
        chr = $1; s = $2; e = $3; name = $4; score = $5; strand = $6
        if (strand == "+") {
          print chr, e - 1, e, name, score, strand
        } else if (strand == "-") {
          print chr, s, s + 1, name, score, strand
        } else {
          mid = int((s + e) / 2)
          print chr, mid, mid + 1, name, score, "."
        }
      }' \
    > "${out_bed}"

  LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n "${out_bed}" -o "${out_bed}"
}

make_intergenic_bed_sorted() {
  # Extract reads that do not overlap any gene body
  local reads_bed="$1"
  local genes_bed="$2"
  local out_bed="$3"

  bedtools intersect -sorted -a "${reads_bed}" -b "${genes_bed}" -v > "${out_bed}"
  LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n "${out_bed}" -o "${out_bed}"
}

compute_distance_histogram() {
  # Compute distance histograms stratified by gene class
  local sample="$1"
  local reads_bed="$2"
  local genes_5p_bed="$3"
  local genes_3p_bed="$4"
  local out_hist="$5"
  local bin="$6"
  local minv="$7"
  local maxv="$8"
  local total_good_reads="$9"

  local tmp5 tmp3
  tmp5=$(mktemp)
  tmp3=$(mktemp)

  # Find closest 5' and 3' gene endpoints for each read
  bedtools closest -sorted -a "${reads_bed}" -b "${genes_5p_bed}" -d > "${tmp5}"
  bedtools closest -sorted -a "${reads_bed}" -b "${genes_3p_bed}" -d > "${tmp3}"

  # Combine and bin distances
  paste "${tmp5}" "${tmp3}" \
    | awk -v SAMPLE="${sample}" -v TOTAL="${total_good_reads}" \
          -v BIN="${bin}" -v MIN="${minv}" -v MAX="${maxv}" '
      BEGIN {
        FS = OFS = "\t"
        minBin = int(MIN / BIN)
        maxBin = int((MAX - 1) / BIN)

        # Gene classes to track
        keep["protein_coding_3UTR"] = 1
        keep["protein_coding_no3UTR"] = 1
        keep["lncRNA_3UTR"] = 1
        keep["lncRNA_no3UTR"] = 1

        regions["5p"] = 1
        regions["3p"] = 1
      }

      function parse_label(lbl,    a, n, bt, utr, grp) {
        n = split(lbl, a, "\\|")
        bt = (n >= 2 ? a[2] : "other")
        utr = (n >= 3 ? a[3] : "no3UTR")
        grp = bt "_" utr
        return grp
      }

      {
        # Columns from paste: 1-6 read, 7-12 gene5p, 13 dist5; 14-19 read, 20-25 gene3p, 26 dist3
        d5 = $(13) + 0
        d3 = $(26) + 0

        # Choose nearer endpoint; ties go to 5p
        region = "3p"; dist = d3; glabel = $(23)
        if (d5 <= d3) { region = "5p"; dist = d5; glabel = $(10) }

        grp = parse_label(glabel)
        if (!(grp in keep)) next

        signed = (region == "5p" ? -dist : dist)
        if (signed < MIN || signed > MAX) next

        b = int(signed / BIN)
        key = grp SUBSEP region SUBSEP b
        counts[key]++
      }

      END {
        print "sample", "gene_group", "region", "bin_start_bp", "bin_end_bp", "bin_mid_bp", "count", "total_good_reads"

        for (g in keep) {
          for (r in regions) {
            for (b = minBin; b <= maxBin; b++) {
              key = g SUBSEP r SUBSEP b
              c = (key in counts ? counts[key] : 0)
              start = b * BIN
              end = start + BIN
              mid = start + BIN / 2.0
              print SAMPLE, g, r, start, end, mid, c, TOTAL
            }
          }
        }
      }
    ' > "${out_hist}"

  rm -f "${tmp5}" "${tmp3}"
}

# -----------------------------------------------------------------------------
# Step 2: Process BAM files
# -----------------------------------------------------------------------------

shopt -s nullglob

for BAM in "${BAM_DIR}"/*.bam; do
  base=$(basename "${BAM}")
  sample="${base%.bam}"

  echo "[INFO] Processing sample: ${sample}"

  FILTERED_BAM="${OUT_DIR}/filtered_bam/${sample}.filtered.bam"
  READS_3P_BED="${OUT_DIR}/beds/${sample}.3p.bed"
  INTERGENIC_BED="${OUT_DIR}/beds/${sample}.intergenic.bed"
  HIST_OUTPUT="${OUT_DIR}/hist/${sample}.intergenic.updown.by_geneclass.hist.tsv"

  # Filter BAM (reuse if valid)
  if [ -f "${FILTERED_BAM}" ] && samtools quickcheck -q "${FILTERED_BAM}"; then
    echo "    Reusing existing filtered BAM"
  else
    [ -f "${FILTERED_BAM}" ] && rm -f "${FILTERED_BAM}"
    filter_bam "${BAM}" "${FILTERED_BAM}"
  fi

  total_good_reads=$(samtools view -c "${FILTERED_BAM}")

  echo "    Extracting read 3' endpoints"
  bam_to_3prime_bed_sorted "${FILTERED_BAM}" "${READS_3P_BED}"

  echo "    Identifying intergenic reads"
  make_intergenic_bed_sorted "${READS_3P_BED}" "${GENE_BED}" "${INTERGENIC_BED}"

  echo "    Computing distance histograms"
  compute_distance_histogram "${sample}" "${INTERGENIC_BED}" "${GENE_5P_BED}" "${GENE_3P_BED}" \
    "${HIST_OUTPUT}" "${BIN}" "${MIN}" "${MAX}" "${total_good_reads}"
done

echo "[INFO] Pipeline complete. Outputs in: ${OUT_DIR}"
