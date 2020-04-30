if [ $# -lt 3 ]
then
  echo Usage: $0 R1.fq R2.fq out_prefix threads
  exit 1
fi

threads=${5:-1}
R1=$1
R2=$2
OUT_PREFIX=$3
RESUME=${6:-overwrite}
OUTDIR=$(dirname $OUT_PREFIX)
TRIMMED_R1=$OUTDIR/trimmed_$(basename $R1)
TRIMMED_R2=$OUTDIR/trimmed_$(basename $R2)
ref=${4:-/home/bl177/lee/boram/ref/hisat2/grch37_snp_tran/genome_snp_tran}

mkdir -p $OUTDIR

nl="\\\\\n"
fastpCMD="fastp --correction $nl
  -i $R1 $nl
  -I $R2 $nl
  -o $TRIMMED_R1 $nl
  -O $TRIMMED_R2 $nl
  --json ${OUT_PREFIX}.fastp.json $nl
  --thread $threads"

cmd="hisat2 $nl
  --sp 1,0 $nl
  --score-min L,0,-0.5 $nl
  --pen-canintronlen S,9,0.1 $nl
  --pen-noncanintronlen S,9,0.1 $nl
  --max-intronlen 90000 $nl
  --dta $nl
  -k 10 $nl
  --secondary $nl
  --threads $threads $nl
  -x $ref $nl
  -1 $TRIMMED_R1 $nl
  -2 $TRIMMED_R2 $nl
  --time $nl
  --met-file ${OUT_PREFIX}.hisat2.metrics.txt $nl
  -S /dev/stdout | samtools view -Sb -o ${OUT_PREFIX}.hisat2.presort.bam -"

RunCMD() {
  OUTFILE=$1
  if [[ $RESUME = resume && -s $OUTFILE && -e ${OUTFILE}.done ]]
  then
    echo
    echo $OUTFILE exists. Skip.
    return 0
  fi
  shift
  rm -f ${OUTFILE}.fail ${OUTFILE}.done
  CMD=$@
  EXE=${CMD//\\\\\\n/}
  echo
  echo [`date`]
  echo -e $CMD
  if echo $EXE | grep "|" -q
  then
    ${EXE%|*} | ${EXE#*|}
  else
    $EXE
  fi
  ES=$?
  echo
  echo "Exit status: $ES"
  if [ $ES -ne 0 -o ! -s $OUTFILE ]
  then
    echo -e $CMD > ${OUTFILE}.fail
    exit $ES
  fi
  echo -e $CMD > ${OUTFILE}.done
  return $ES
}

RunCMD $TRIMMED_R1 $fastpCMD
RunCMD ${OUT_PREFIX}.hisat2.presort.bam $cmd
RunCMD ${OUT_PREFIX}.hisat2.bam samtools sort -@ $threads -o ${OUT_PREFIX}.hisat2.bam ${OUT_PREFIX}.hisat2.presort.bam
RunCMD ${OUT_PREFIX}.hisat2.bam.bai samtools index -@ $threads ${OUT_PREFIX}.hisat2.bam
rm -f $TRIMMED_R1 $TRIMMED_R2
rm -f ${OUT_PREFIX}.hisat2.presort.bam
