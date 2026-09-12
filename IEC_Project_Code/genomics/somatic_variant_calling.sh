#!/usr/bin/env bash
# Somatic variant calling
# Align exome reads, recalibrate bases, and perform tumor-only Mutect2 filtering and annotation.

set -euo pipefail

SAMPLE="${1:?SAMPLE required}"
R1="${2:?R1 fastq.gz required}"
R2="${3:?R2 fastq.gz required}"

THREADS="${THREADS:-8}"

BWA_MEM2="${BWA_MEM2:-bwa-mem2}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
GATK="${GATK:-gatk}"
FUNCOTATOR_DS="${FUNCOTATOR_DS:-references/funcotator}"
DO_FUNCOTATOR="${DO_FUNCOTATOR:-1}"

export BWA_USE_SIMD="${BWA_USE_SIMD:-avx2}"

REF="${REF:-references/genome.fasta}"
DICT="${DICT:-${REF%.fasta}.dict}"

DBSNP="${DBSNP:-references/dbsnp.vcf}"
MILLS_INDEL="${MILLS_INDEL:-references/known_indels.vcf.gz}"
ONEKG_INDEL="${ONEKG_INDEL:-references/additional_indels.vcf}"

GNOMAD_AFONLY="${GNOMAD_AFONLY:-references/germline_allele_frequencies.vcf}"
PON_VCF="${PON_VCF:-references/panel_of_normals.vcf}"

DESIGN_DIR="${DESIGN_DIR:-references/capture_design}"
REGIONS_BED="${REGIONS_BED:-${DESIGN_DIR}/Regions.sorted.bed}"
PROBES_BED="${PROBES_BED:-${DESIGN_DIR}/MergedProbes.sorted.bed}"
PADDED_BED="${PADDED_BED:-${DESIGN_DIR}/Padded.sorted.bed}"

OUT_ROOT="${OUT_ROOT:-results/genomics/somatic_variants}"
SOUT="${OUT_ROOT}/${SAMPLE}"
ALIGN_DIR="${SOUT}/align"
QC_DIR="${SOUT}/qc"
VCF_DIR="${SOUT}/mutect2"
mkdir -p "$ALIGN_DIR" "$QC_DIR" "$VCF_DIR"

LOG="${SOUT}/${SAMPLE}.pipeline.log"
exec > >(tee -i "$LOG") 2>&1

echo "=== $(date) START ==="
echo "Sample:   $SAMPLE"
echo "R1:       $R1"
echo "R2:       $R2"
echo "Threads:  $THREADS"
echo "Out:      $SOUT"
echo "REF:      $REF"
echo "Design:   $DESIGN_DIR"

for f in "$R1" "$R2" "$REF" "$DICT" "$REGIONS_BED" "$PROBES_BED" "$PADDED_BED" "$GNOMAD_AFONLY"; do
  [[ -s "$f" ]] || {
    echo "[FATAL] missing: $f" >&2
    exit 2
  }
done
for x in "$BWA_MEM2" "$SAMTOOLS" "$BCFTOOLS"; do
  command -v "$x" >/dev/null 2>&1 || {
    echo "[FATAL] missing executable: $x" >&2
    exit 2
  }
done

SORT_BAM="${ALIGN_DIR}/${SAMPLE}.sorted.bam"

if [[ ! -f "${REF}.bwt.2bit.64" && ! -f "${REF}.bwt.2bit.32" ]]; then
  echo "[FATAL] bwa-mem2 index not found next to reference. Build once with:"
  echo "  ${BWA_MEM2} index ${REF}"
  exit 2
fi

RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}\tPU:${SAMPLE}"

echo "[Align] bwa-mem2 mem -> sort"
"$BWA_MEM2" mem -M -t "$THREADS" -R "$RG" "$REF" "$R1" "$R2" |
  "$SAMTOOLS" sort -@ "$THREADS" -o "$SORT_BAM" -
"$SAMTOOLS" index -@ "$THREADS" "$SORT_BAM"

command -v "$GATK" >/dev/null 2>&1 || {
  echo "[FATAL] GATK not found" >&2
  exit 2
}

ensure_gatk_index() {
  local vcf="$1"
  if [[ ! -s "$vcf" ]]; then return 0; fi
  if [[ -s "${vcf}.tbi" || -s "${vcf}.csi" || -s "${vcf}.idx" ]]; then return 0; fi
  echo "[Index] IndexFeatureFile $vcf"
  "$GATK" IndexFeatureFile -I "$vcf" || true
  if [[ "$vcf" == *.vcf.gz && ! -s "${vcf}.tbi" && ! -s "${vcf}.csi" ]]; then
    "$BCFTOOLS" index -f -t "$vcf" || true
  fi
}

ensure_gatk_index "$GNOMAD_AFONLY"
ensure_gatk_index "$DBSNP"
ensure_gatk_index "$MILLS_INDEL"
ensure_gatk_index "$ONEKG_INDEL"

if [[ -s "$PON_VCF" ]]; then
  if [[ ! -s "${PON_VCF}.tbi" && ! -s "${PON_VCF}.csi" ]]; then
    "$BCFTOOLS" index -f -t "$PON_VCF" || true
  fi
fi

ref_has_chr=0
if grep -q $'^@SQ\tSN:chr' "$DICT"; then ref_has_chr=1; fi

normalize_bed() {
  local in_bed="$1"
  local out_bed="$2"
  awk -v REF_HAS_CHR="$ref_has_chr" 'BEGIN{OFS="\t"}
    $0 ~ /^#/ || $1 ~ /^(track|browser)$/ { next }
    NF<3 { next }
    {
      c=$1; s=$2; e=$3;
      if (REF_HAS_CHR==1) { if (c !~ /^chr/) c="chr" c }
      else { sub(/^chr/,"",c) }
      if (s ~ /^[0-9]+$/ && e ~ /^[0-9]+$/ && e>s) print c,s,e
    }' "$in_bed" | sort -k1,1V -k2,2n -k3,3n >"$out_bed"
}

NREG="${QC_DIR}/Regions.norm.bed"
NPRO="${QC_DIR}/MergedProbes.norm.bed"
NPAD="${QC_DIR}/Padded.norm.bed"
normalize_bed "$REGIONS_BED" "$NREG"
normalize_bed "$PROBES_BED" "$NPRO"
normalize_bed "$PADDED_BED" "$NPAD"

REG_IL="${QC_DIR}/Regions.interval_list"
PRO_IL="${QC_DIR}/MergedProbes.interval_list"
"$GATK" BedToIntervalList -I "$NREG" -O "$REG_IL" -SD "$DICT"
"$GATK" BedToIntervalList -I "$NPRO" -O "$PRO_IL" -SD "$DICT"

DEDUP_BAM="${ALIGN_DIR}/${SAMPLE}.dedup.bam"
DUP_METRICS="${QC_DIR}/${SAMPLE}.markdup.metrics.txt"

echo "[GATK] MarkDuplicates"
"$GATK" MarkDuplicates \
  -I "$SORT_BAM" \
  -O "$DEDUP_BAM" \
  -M "$DUP_METRICS" \
  --CREATE_INDEX true \
  --VALIDATION_STRINGENCY SILENT

BQSR_TABLE="${ALIGN_DIR}/${SAMPLE}.bqsr.table"
BQSR_BAM="${ALIGN_DIR}/${SAMPLE}.bqsr.bam"

DO_BQSR=1
for v in "$DBSNP" "$MILLS_INDEL" "$ONEKG_INDEL"; do
  [[ -s "$v" ]] || DO_BQSR=0
done

if [[ "$DO_BQSR" == "1" ]]; then
  echo "[GATK] BaseRecalibrator"
  "$GATK" BaseRecalibrator \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --known-sites "$DBSNP" \
    --known-sites "$MILLS_INDEL" \
    --known-sites "$ONEKG_INDEL" \
    -O "$BQSR_TABLE"

  echo "[GATK] ApplyBQSR"
  "$GATK" ApplyBQSR \
    -R "$REF" \
    -I "$DEDUP_BAM" \
    --bqsr-recal-file "$BQSR_TABLE" \
    -O "$BQSR_BAM"
  "$SAMTOOLS" index -@ "$THREADS" "$BQSR_BAM"
else
  echo "[WARN] Known-sites missing; skipping BQSR. Using MarkDuplicates BAM as final."
  BQSR_BAM="$DEDUP_BAM"
fi

HS_METRICS="${QC_DIR}/${SAMPLE}.hs_metrics.txt"
echo "[GATK] CollectHsMetrics"
"$GATK" CollectHsMetrics \
  -I "$BQSR_BAM" \
  -O "$HS_METRICS" \
  -R "$REF" \
  --BAIT_INTERVALS "$PRO_IL" \
  --TARGET_INTERVALS "$REG_IL"

UNFILT="${VCF_DIR}/${SAMPLE}.unfiltered.vcf.gz"
F1R2="${VCF_DIR}/${SAMPLE}.f1r2.tar.gz"
OB="${VCF_DIR}/${SAMPLE}.read-orientation-model.tar.gz"
PILEUPS="${VCF_DIR}/${SAMPLE}.pileups.table"
CONTAM="${VCF_DIR}/${SAMPLE}.contamination.table"
FILT="${VCF_DIR}/${SAMPLE}.filtered.vcf.gz"
PASS="${VCF_DIR}/${SAMPLE}.PASS.vcf.gz"

PON_ARGS=()
if [[ -s "$PON_VCF" && (-s "${PON_VCF}.tbi" || -s "${PON_VCF}.idx") ]]; then
  PON_ARGS=(--panel-of-normals "$PON_VCF")
fi

echo "[GATK] Mutect2"
"$GATK" Mutect2 \
  -R "$REF" \
  -I "$BQSR_BAM" -tumor "$SAMPLE" \
  --germline-resource "$GNOMAD_AFONLY" \
  "${PON_ARGS[@]}" \
  -L "$NPAD" \
  --f1r2-tar-gz "$F1R2" \
  --max-mnp-distance 0 \
  -O "$UNFILT"

"$GATK" LearnReadOrientationModel -I "$F1R2" -O "$OB"

"$GATK" GetPileupSummaries \
  -I "$BQSR_BAM" \
  -V "$GNOMAD_AFONLY" \
  -L "$NPAD" \
  -O "$PILEUPS"

"$GATK" CalculateContamination -I "$PILEUPS" -O "$CONTAM"

"$GATK" FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILT" \
  --contamination-table "$CONTAM" \
  --ob-priors "$OB" \
  -O "$FILT"

echo "[VCF] PASS only"
"$BCFTOOLS" view -f PASS -Oz -o "$PASS" "$FILT"
"$BCFTOOLS" index -f -t "$PASS"

if [[ "${DO_FUNCOTATOR:-1}" == "1" ]]; then
  [[ -d "$FUNCOTATOR_DS" ]] || {
    echo "[FATAL] FUNCOTATOR_DS not found: $FUNCOTATOR_DS"
    exit 2
  }
  [[ -s "$PASS" ]] || {
    echo "[FATAL] PASS VCF missing: $PASS"
    exit 2
  }

  VCF_TO_ANNOTATE="${VCF_TO_ANNOTATE:-$PASS}"
  [[ -s "$VCF_TO_ANNOTATE" ]] || {
    echo "[FATAL] VCF_TO_ANNOTATE missing: $VCF_TO_ANNOTATE"
    exit 2
  }

  MAF_OUT="${VCF_DIR}/${SAMPLE}.PASS.funcotator.maf"
  echo "[Funcotator] Annotating $VCF_TO_ANNOTATE -> $MAF_OUT"

  "$GATK" Funcotator \
    --variant "$VCF_TO_ANNOTATE" \
    --reference "$REF" \
    --ref-version hg19 \
    --data-sources-path "$FUNCOTATOR_DS" \
    --output "$MAF_OUT" \
    --output-file-format MAF

  echo "[Funcotator] Wrote: $MAF_OUT"
fi

echo "=== $(date) DONE ==="
echo "Final BAM:  $BQSR_BAM"
echo "HS metrics: $HS_METRICS"
echo "PASS VCF:   $PASS"
echo "MAF:        ${VCF_DIR}/${SAMPLE}.PASS.funcotator.maf"
