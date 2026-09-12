#!/usr/bin/env bash
# Vector junction extraction
# Align paired reads to a combined genome/vector reference and extract junction evidence.

set -euo pipefail

SAMPLE="${1:?SAMPLE required}"
R1="${2:-}"
R2="${3:-}"

THREADS="${THREADS:-16}"
OUT_ROOT="${OUT_ROOT:-results/vector_integration}"

REF="${REF:-references/genome_plus_vector.fa}"
VECTOR_REGEX="${VECTOR_REGEX:-^(ciltacel|axicel|CAR_VECTOR)$}"

MIN_MAPQ="${MIN_MAPQ:-20}"
BIN_SIZE="${BIN_SIZE:-10}"
DO_TRIM="${DO_TRIM:-1}"

SAMTOOLS="${SAMTOOLS:-samtools}"
FASTP="${FASTP:-fastp}"
GATK_JAR="${GATK_JAR:?Set GATK_JAR to a GATK 4 jar}"
JAVA="${JAVA:-java}"

BWA="${BWA:-bwa}"

OUTDIR="$OUT_ROOT/$SAMPLE"
QC="$OUTDIR/qc"
ALIGN="$OUTDIR/align"
JUNC="$OUTDIR/junctions"
LOGS="$OUTDIR/logs"
mkdir -p "$QC" "$ALIGN" "$JUNC" "$LOGS"

LOG="$LOGS/${SAMPLE}.pipeline.log"
exec > >(tee -i "$LOG") 2>&1

echo "=== $(date) START WGS CAR insertion pipeline with checkpoints ==="
echo "Sample: $SAMPLE"
echo "REF: $REF"
echo "OUTDIR: $OUTDIR"
echo "THREADS: $THREADS"
echo "VECTOR_REGEX: $VECTOR_REGEX"
echo "MIN_MAPQ: $MIN_MAPQ"
echo "BIN_SIZE: $BIN_SIZE"
echo "DO_TRIM: $DO_TRIM"

command -v "$SAMTOOLS" >/dev/null 2>&1 || {
  echo "[FATAL] samtools missing: $SAMTOOLS" >&2
  exit 2
}
[[ -s "$REF" && -s "${REF}.fai" ]] || {
  echo "[FATAL] REF or REF.fai missing: $REF" >&2
  exit 2
}
DICT="${REF%.fa}.dict"
[[ -s "$DICT" ]] || DICT="${REF%.fasta}.dict"
[[ -s "$DICT" ]] || {
  echo "[FATAL] REF dict missing (.dict) near $REF" >&2
  exit 2
}
[[ -n "$GATK_JAR" && -s "$GATK_JAR" ]] || {
  echo "[FATAL] GATK jar not found" >&2
  exit 2
}

command -v "$BWA" >/dev/null 2>&1 || {
  echo "[FATAL] bwa not found" >&2
  exit 2
}

RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1"

SORTBAM="$ALIGN/${SAMPLE}.sorted.bam"
CAR_BAM="$ALIGN/${SAMPLE}.car_related.bam"
DEDUP="$ALIGN/${SAMPLE}.dedup.bam"

if [[ -s "$DEDUP" && -s "${DEDUP}.bai" ]]; then
  echo "[CHECKPOINT] Found dedup bam+index. Skipping trim/alignment/subset/MarkDuplicates."
else
  if [[ -s "$CAR_BAM" && -s "${CAR_BAM}.bai" ]]; then
    echo "[CHECKPOINT] Found car_related bam+index. Skipping trim/alignment/subset."
  else
    if [[ -s "$SORTBAM" && -s "${SORTBAM}.bai" ]]; then
      echo "[CHECKPOINT] Found sorted bam+index. Skipping trim/alignment."
    else
      [[ -n "$R1" && -n "$R2" ]] || {
        echo "[FATAL] Need R1/R2 unless sorted/car_related/dedup exists." >&2
        exit 2
      }
      [[ -s "$R1" && -s "$R2" ]] || {
        echo "[FATAL] FASTQs missing" >&2
        exit 2
      }

      IN1="$R1"
      IN2="$R2"
      if [[ "$DO_TRIM" == "1" ]]; then
        TR1="$ALIGN/${SAMPLE}.R1.trim.fastq.gz"
        TR2="$ALIGN/${SAMPLE}.R2.trim.fastq.gz"
        if [[ -s "$TR1" && -s "$TR2" ]]; then
          echo "[CHECKPOINT] Found trimmed FASTQs. Skipping fastp."
        else
          echo "[1/4] fastp trimming"
          "$FASTP" -w "$THREADS" \
            -i "$R1" -I "$R2" \
            -o "$TR1" -O "$TR2" \
            --detect_adapter_for_pe \
            -h "$QC/${SAMPLE}.fastp.html" \
            -j "$QC/${SAMPLE}.fastp.json"
        fi
        IN1="$TR1"
        IN2="$TR2"
      else
        echo "[1/4] trimming disabled (DO_TRIM=0)"
      fi

      echo "[2/4] Align (bwa mem) → sort/index"
      "$BWA" mem -t "$THREADS" -M -R "$RG" "$REF" "$IN1" "$IN2" |
        "$SAMTOOLS" sort -@ "$THREADS" -o "$SORTBAM" -
      "$SAMTOOLS" index -@ "$THREADS" "$SORTBAM"
    fi

    echo "[2.5/4] Subset to CAR-related readnames"
    VECLIST="$LOGS/${SAMPLE}.vector_contigs.txt"
    QNLIST="$LOGS/${SAMPLE}.vector_readnames.txt"

    "$SAMTOOLS" idxstats "$SORTBAM" | awk -v re="$VECTOR_REGEX" '$1 ~ re {print $1}' >"$VECLIST"
    [[ -s "$VECLIST" ]] || {
      echo "[FATAL] No contigs match VECTOR_REGEX=$VECTOR_REGEX" >&2
      exit 2
    }

    "$SAMTOOLS" view -@ "$THREADS" -F 3852 "$SORTBAM" $(tr '\n' ' ' <"$VECLIST") |
      awk -v mq="$MIN_MAPQ" '$5>=mq {print $1}' | sort -u >"$QNLIST"
    [[ -s "$QNLIST" ]] || {
      echo "[FATAL] No reads mapped to vector contigs (after MQ filter)." >&2
      exit 2
    }

    "$SAMTOOLS" view -@ "$THREADS" -N "$QNLIST" -bh "$SORTBAM" >"$CAR_BAM"
    "$SAMTOOLS" index -@ "$THREADS" "$CAR_BAM"
  fi

  echo "[3/4] MarkDuplicates (GATK)"
  METRICS="$ALIGN/${SAMPLE}.markdups.metrics.txt"
  "$JAVA" -jar "$GATK_JAR" MarkDuplicates \
    -I "$CAR_BAM" \
    -O "$DEDUP" \
    -M "$METRICS" \
    --CREATE_INDEX true
  "$SAMTOOLS" index -@ "$THREADS" "$DEDUP"
fi

echo "[4/4] Extract vector-genome junction evidence"
DISC="$JUNC/${SAMPLE}.discordant_pairs.tsv"
SPLT="$JUNC/${SAMPLE}.split_reads.tsv"
MERG="$JUNC/${SAMPLE}.junctions.merged.tsv"
CLST="$JUNC/${SAMPLE}.clusters.tsv"

EXCL_FLAG=3852

"$SAMTOOLS" idxstats "$DEDUP" | awk -v re="$VECTOR_REGEX" '$1 ~ re {found=1} END{exit(found?0:1)}' ||
  {
    echo "[FATAL] No contigs match VECTOR_REGEX=$VECTOR_REGEX in BAM header" >&2
    exit 2
  }

"$SAMTOOLS" view -@ "$THREADS" -f 1 -F "$EXCL_FLAG" "$DEDUP" |
  awk -v re="$VECTOR_REGEX" -v mq="$MIN_MAPQ" 'BEGIN{OFS="\t"}
  {
    q=$1; flag=$2; r=$3; pos=$4; mapq=$5; cig=$6; rn=$7; pn=$8; tlen=$9;

    if(rn=="=") rn=r;

    if(mapq < mq) next;
    if(rn=="*" || pn==0) next;

    r_is_vec  = (r  ~ re);
    rn_is_vec = (rn ~ re);

    if(r_is_vec == rn_is_vec) next;

    if(r_is_vec && !rn_is_vec){
      print q, r, pos, mapq, rn, pn, flag, cig, tlen;
    } else if(!r_is_vec && rn_is_vec){
      print q, rn, pn, mapq, r, pos, flag, cig, tlen;
    }
  }' >"$DISC"

"$SAMTOOLS" view -@ "$THREADS" -F "$EXCL_FLAG" "$DEDUP" |
  awk -v re="$VECTOR_REGEX" -v mq="$MIN_MAPQ" 'BEGIN{OFS="\t"}
  {
    q=$1; r=$3; pos=$4; mapq=$5; cig=$6;
    if(mapq < mq) next;

    sa="";
    for(i=12;i<=NF;i++){
      if($i ~ /^SA:Z:/){ sa=$i; sub(/^SA:Z:/,"",sa); break; }
    }
    if(sa=="") next;

    n=split(sa, arr, ";");
    for(k=1;k<=n;k++){
      if(arr[k]=="") continue;
      m=split(arr[k], f, ",");
      if(m < 5) continue;

      sa_r=f[1]; sa_pos=f[2]; sa_str=f[3]; sa_cig=f[4]; sa_mq=f[5];
      if(sa_mq < mq) continue;

      r_is_vec  = (r    ~ re);
      sa_is_vec = (sa_r ~ re);

      if( (r_is_vec && !sa_is_vec) || (!r_is_vec && sa_is_vec) ){
        print q, r, pos, mapq, cig, sa_r, sa_pos, sa_str, sa_cig, sa_mq;
      }
    }
  }' >"$SPLT"

{
  echo -e "sample\tevidence\tqname\tvec_contig\tvec_pos\tgen_contig\tgen_pos\tvec_mapq\tgen_mapq\tdetails"
  awk -v s="$SAMPLE" 'BEGIN{OFS="\t"}
    {print s,"DISCORDANT",$1,$2,$3,$5,$6,$4,$4,("flag="$7";cigar="$8";tlen="$9)}' "$DISC"
  awk -v s="$SAMPLE" -v re="$VECTOR_REGEX" 'BEGIN{OFS="\t"}
    {
      q=$1; pr=$2; ppos=$3; pmq=$4; pcig=$5;
      sr=$6; spos=$7; sstr=$8; scig=$9; smq=$10;

      pr_is_vec=(pr ~ re);
      sr_is_vec=(sr ~ re);

      if(pr_is_vec && !sr_is_vec){
        print s,"SPLIT",q,pr,ppos,sr,spos,pmq,smq,("pCIGAR="pcig";saStr="sstr";saCIGAR="scig);
      } else if(!pr_is_vec && sr_is_vec){
        print s,"SPLIT",q,sr,spos,pr,ppos,smq,pmq,("pCIGAR="pcig";saStr="sstr";saCIGAR="scig);
      }
    }' "$SPLT"
} >"$MERG"

TMPCL="$JUNC/${SAMPLE}.clusters.tmp"
awk -v bin="$BIN_SIZE" 'BEGIN{FS=OFS="\t"} NR==1{next}
  {
    chr=$6; pos=$7+0;
    b=int(pos/bin)*bin;
    key=chr OFS b;
    total[key]++;

    if($2=="DISCORDANT") dp[key]++; else if($2=="SPLIT") sp[key]++;
  }
  END{
    print "gen_contig","bin_start","total_support","n_discordant","n_split";
    for(k in total){
      d=(k in dp)?dp[k]:0;
      s=(k in sp)?sp[k]:0;
      print k,total[k],d,s;
    }
  }' "$MERG" >"$TMPCL"

{
  head -n 1 "$TMPCL"
  tail -n +2 "$TMPCL" | sort -k1,1V -k2,2n -k3,3nr
} >"$CLST"
rm -f "$TMPCL"

echo "=== $(date) DONE ==="
echo "DEDUP:    $DEDUP"
echo "DISC:     $DISC ($(wc -l <"$DISC") lines)"
echo "SPLIT:    $SPLT ($(wc -l <"$SPLT") lines)"
echo "CLUSTERS: $CLST ($(($(wc -l <"$CLST") - 1)) bins)"
