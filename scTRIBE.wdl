workflow main {
  String root
  String fastq1
  String fastq2
  String outdir
  String PISA
  String minimap2
  String juncbed
  String bcftools
  String sambamba
  String ID
  String gtf
  String ref
  String config
  String ?lib
  call makedir {
    input:
    Dir=outdir
  }
  call parseFastq {
    input:
    lib=lib,
    PISA=PISA,
    config=config,
    fastq1=fastq1,
    fastq2=fastq2,
    outdir=makedir.Outdir,
    root=root,
  }
  call fastq2bam {
    input:
    lib=lib,
    PISA=PISA,
    fastq=parseFastq.fastq,
    outdir=outdir,
    minimap2=minimap2,
    juncbed=juncbed,
    gtf=gtf,
    ref=ref,
    root=root
  }
  call sortBam {
    input:
    lib=lib,
    PISA=PISA,
    bam1=fastq2bam.bam,
    sambamba=sambamba,
    gtf=gtf,
    root=root,
    outdir=outdir
  }
  call callVar {
    input:
    bam=sortBam.bam,
    PISA=PISA,
    sambamba=sambamba,
    bcftools=bcftools,
    root=root,
    outdir=outdir,
    ref=ref,
    lib=lib
  }
  call cellCount {
    input:
    bam=sortBam.bam,
    PISA=PISA,
    outdir=outdir,
    root=root,
    lib=lib
  }
  call countMatrix {
    input:
    lib=lib,
    root=root,
    PISA=PISA,
    outdir=outdir,
    anno=sortBam.bam
  }
}
task callVar {
  String bam
  String bcftools
  String sambamba
  String root
  String lib
  String ref
  String PISA
  String outdir
  command {

    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${PISA} pick -tag GN -o ${outdir}/temp/pick_GN.bam -@ 4 ${bam}
    ${sambamba} index -t 4 ${outdir}/temp/pick_GN.bam
    ${bcftools} mpileup -Ou -f ${ref} ${outdir}/temp/pick_GN.bam | ${bcftools} call -Ob -mv -o ${outdir}/outs/var.bcf
    ${bcftools} filter -s LowQual -e '%QUAL<20 || DP>10' -O b -o ${outdir}/outs/filter.bcf ${outdir}/outs/var.bcf
  }
}
task countMatrix {
  String root
  String outdir
  String anno
  String PISA
  String ?lib
  command {

    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi
    mkdir -p ${outdir}/outs/raw_gene_expression
    ${PISA} count -@ 10 -tag CB -anno-tag GN -umi UB -outdir ${outdir}/outs/raw_gene_expression ${anno}
    
    echo "[`date +%F` `date +%T`] workflow end" >> ${outdir}/workflowtime.log
  }
  output {
    String matrix = "${outdir}/outs/raw_gene_expression/matrix.mtx.gz"
    File matrix0 = "${outdir}/outs/raw_gene_expression/matrix.mtx.gz"
  }
}

task cellCount {
  String bam
  String outdir
  String root
  String PISA
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${PISA} attrcnt -cb CB -tags UB,GN -dedup -@ 4 -o ${outdir}/temp/cell_stat.txt -all-tags ${outdir}/outs/final.bam
  }
  output {
    String count="${outdir}/temp/cell_stat.txt"
  }
}
task makedir {
  String Dir
  command {
    echo "[`date +%F` `date +%T`] workflow start" > ${Dir}/workflowtime.log
    mkdir -p ${Dir}
    mkdir -p ${Dir}/outs
    mkdir -p ${Dir}/temp
    mkdir -p ${Dir}/report
  }
  output {
    String Outdir="${Dir}"
  }
}
task parseFastq {
  String config
  String fastq1
  String fastq2
  String outdir
  String root
  String PISA
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    # disable multi-threads support
    ${PISA} parse -q 20 -config ${config} -cbdis ${outdir}/temp/barcode_counts_raw.txt -report ${outdir}/report/sequencing_report.csv ${fastq1} ${fastq2} -1 ${outdir}/temp/reads.fq

  }
  output {
    String count="${outdir}/temp/barcode_counts_raw.txt"
    String fastq="${outdir}/temp/reads.fq"
    String sequencingReport="${outdir}/report/sequencing_report.csv"
      }
}

task fastq2bam {
  String fastq  
  String outdir
  String minimap2
  String juncbed
  String ref
  String PISA
  String root
  String gtf
  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${minimap2} -x splice:hq -uf -t 5 -a --junc-bed ${juncbed} ${ref} ${fastq} 1>${outdir}/temp/aln.sam && \
    ${PISA} sam2bam -adjust-mapq -gtf ${gtf} -o ${outdir}/temp/aln.bam -report ${outdir}/report/alignment_report.csv ${outdir}/temp/aln.sam && \
    rm -f ${outdir}/temp/aln.sam
  }
  output {
    String bam="${outdir}/temp/aln.bam"
    File bam0="${outdir}/temp/aln.bam"
    String alnReport="{outdir}/report/alignment_report.csv"
  }
}
task sortBam {
  String bam1
  String sambamba
  String root
  String outdir
  String PISA
  String gtf

  String ?lib
  command {
    if [ -f ${default=abjdbashj lib} ]; then
    source ${lib}
    fi

    ${sambamba} sort -t 10 -o ${outdir}/temp/sorted.bam ${bam1}
    ${PISA} anno -gtf ${gtf} -o ${outdir}/temp/anno.bam -report ${outdir}/report/anno_report.csv ${outdir}/temp/sorted.bam
    ${PISA} corr -tag UR -new-tag UB -cr -@ 10 -tags-block CB,GN -o ${outdir}/temp/corr.bam ${outdir}/temp/anno.bam
    ${PISA} rmdup -tag CB,UR -@ 10 -k -o ${outdir}/outs/final.bam ${outdir}/temp/corr.bam 
  }
  output {
    String bam="${outdir}/outs/final.bam"
  }
}
