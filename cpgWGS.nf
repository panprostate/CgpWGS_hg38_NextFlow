nextflow.enable.dsl=2



process getNames{
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal)

	output:
	tuple val(tumor), val(normal), emit: files
	tuple env(NAME_MT), env(NAME_WT), emit: sample_names
	
    script:
    """
    NAME_MT=`samtools view -H $tumor | grep "^@RG" | cut -f3 | uniq | sed s/SM://`
	NAME_WT=`samtools view -H $normal | grep "^@RG" | cut -f3 | uniq | sed s/SM://`
    """
}

process CN_init{
	errorStrategy 'finish'


	output:
	val true
	
	script:
	"""
	mkdir -p $params.TMP
	awk -v FS='\t' -v OFS='\t' '{print \$1,0,\$2, 2}' ${params.REF_BASE}/genome.fa.fai > ${params.TMP}/norm.cn.bed 
	awk -v FS='\t' -v OFS='\t' '{print \$1,0,\$2, 2}' ${params.REF_BASE}/genome.fa.fai > ${params.TMP}/tum.cn.bed
	"""
}
	    
process cavemanSetup {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)
    val ready

    output:
	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names

    script:
    """
    caveman.pl \
     -r $params.REF_BASE/genome.fa.fai \
     -ig $params.REF_BASE/caveman/HiDepth.tsv \
     -b $params.REF_BASE/caveman/flagging \
     -ab $params.REF_BASE/vagrent \
     -u $params.REF_BASE/caveman \
     -s $params.SPECIES \
     -sa $params.ASSEMBLY \
     -t $task.cpus \
     -st $params.PROTOCOL \
     -tc $params.TMP/tum.cn.bed \
     -nc $params.TMP/norm.cn.bed \
     -td 5 -nd 2 \
     -tb $tumor \
     -nb $normal \
     -c $params.REF_BASE/caveman/flag.vcf.config.WGS.ini \
     -f $params.REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
     -e 350000 \
     -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
     -x $params.CONTIG_EXCLUDE \
     -p setup
     """ 
}

process compareBamGenotypes{
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)


    script:
    "compareBamGenotypes.pl \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped \
    -nb $normal \
    -j $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/genotyped/result.json \
    -tb $tumor \
    -s $params.REF_BASE/general.tsv \
    -g $params.REF_BASE/gender.tsv"
}


process verifyBamHomChkNormal{
	errorStrategy 'finish'

    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    script:
    "verifyBamHomChk.pl -d 25 \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_WT}/contamination \
    -b $normal \
    -t $task.cpus \
    -j $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_WT}/contamination/result.json \
    -s $params.REF_BASE/verifyBamID_snps.vcf.gz"
}

process cavemanSplit{
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    output:
	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names

    script:
    """
    echo "CWD=\${PWD}" 1<> $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
     -r $params.REF_BASE/genome.fa.fai \
     -ig $params.REF_BASE/caveman/HiDepth.tsv \
     -b $params.REF_BASE/caveman/flagging \
     -ab $params.REF_BASE/vagrent \
     -u $params.REF_BASE/caveman \
     -s $params.SPECIES \
     -sa $params.ASSEMBLY \
     -t $task.cpus \
     -st $params.PROTOCOL \
     -tc $params.TMP/tum.cn.bed \
     -nc $params.TMP/norm.cn.bed \
     -td 5 -nd 2 \
     -tb $tumor \
     -nb $normal \
     -c $params.REF_BASE/caveman/flag.vcf.config.WGS.ini \
     -f $params.REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
     -e 350000 \
     -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
     -x $params.CONTIG_EXCLUDE \
     -p split
     """
}

process ascat {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    output:
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/${NAME_MT}.samplestatistics.txt", emit: sample_stats
    val "${params.OUTPUT_DIR}/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat/${NAME_MT}.copynumber.caveman.csv", emit: CN

    script:
    """
    ASCAT_ADD_ARGS=''

    if [ ! -z ${params.ASCAT_PURITY} ]; then
    ASCAT_ADD_ARGS="\$ASCAT_ADD_ARGS -pu $params.ASCAT_PURITY"
    fi
    if [ ! -z ${params.ASCAT_PLOIDY} ]; then
    ASCAT_ADD_ARGS="\$ASCAT_ADD_ARGS -pi $params.ASCAT_PLOIDY"
    fi

    ascat.pl \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/ascat \
    -t $tumor \
    -n $normal \
    -sg $params.REF_BASE/ascat/SnpGcCorrections.tsv \
    -r $params.REF_BASE/genome.fa \
    -q 20 \
    -g L \
    -l $params.REF_BASE/gender.tsv \
    -rs $params.SPECIES \
    -ra $params.ASSEMBLY \
    -pr $params.PROTOCOL \
    -pl ILLUMINA \
    -c $task.cpus \
    -force \
    \$ASCAT_ADD_ARGS
    """
}

process brass_input {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    output:
	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names

    script:
    "brass.pl -j 4 -k 4 -c $task.cpus \
    -d $params.REF_BASE/brass/HiDepth.bed.gz \
    -f $params.REF_BASE/brass/brass_np.groups.gz \
    -g $params.REF_BASE/genome.fa \
    -s $params.SPECIES -as $params.ASSEMBLY -pr $params.PROTOCOL -pl ILLUMINA \
    -g_cache $params.REF_BASE/vagrent/vagrent.cache.gz \
    -vi $params.REF_BASE/brass/viral.genomic.fa.2bit \
    -mi $params.REF_BASE/brass/all_ncbi_bacteria \
    -b $params.REF_BASE/brass/500bp_windows.gc.bed.gz \
    -ct $params.REF_BASE/brass/CentTelo.tsv \
    -cb $params.REF_BASE/brass/cytoband.txt \
    -t $tumor \
    -n $normal \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass \
    -p input"
}

process brass_cover {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    output:
	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names

    script:
    "brass.pl -j 4 -k 4 -c $task.cpus \
    -d $params.REF_BASE/brass/HiDepth.bed.gz \
    -f $params.REF_BASE/brass/brass_np.groups.gz \
    -g $params.REF_BASE/genome.fa \
    -s $params.SPECIES -as $params.ASSEMBLY -pr $params.PROTOCOL -pl ILLUMINA \
    -g_cache $params.REF_BASE/vagrent/vagrent.cache.gz \
    -vi $params.REF_BASE/brass/viral.genomic.fa.2bit \
    -mi $params.REF_BASE/brass/all_ncbi_bacteria \
    -b $params.REF_BASE/brass/500bp_windows.gc.bed.gz \
    -ct $params.REF_BASE/brass/CentTelo.tsv \
    -cb $params.REF_BASE/brass/cytoband.txt \
    -t $tumor \
    -n $normal \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass \
    -p cover"
}

process prepAscat {
	errorStrategy 'finish'

    input:
    tuple val(NAME_MT), val(NAME_WT)
    val ASCAT_CN
    val MT_sampleStats

    output:
    tuple val("${params.TMP}/${NAME_MT}_vs_${NAME_WT}/norm.cn.bed"), val("${params.TMP}/${NAME_MT}_vs_${NAME_WT}/tum.cn.bed"), emit: CN_files
    env NORM_CONTAM, emit: contam

    script:
    """
    mkdir -p ${params.TMP}/${NAME_MT}_vs_${NAME_WT}/
    perl -ne '@F=(split q{,}, \$_)[1,2,3,4]; \$F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > ${params.TMP}/${NAME_MT}_vs_${NAME_WT}/norm.cn.bed
    perl -ne '@F=(split q{,}, \$_)[1,2,3,6]; \$F[1]-1; print join("\t",@F)."\n";' < $ASCAT_CN > ${params.TMP}/${NAME_MT}_vs_${NAME_WT}/tum.cn.bed
    NORM_CONTAM=`perl -ne 'if(m/^rho\s(.+)\n/){print 1-\$1;}' $MT_sampleStats`
    """
}

process pindel {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)

    output:
    tuple val(NAME_MT), val(NAME_WT), emit: sample_names
    val "$params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.germline.bed", emit: germline

    script:
    "pindel.pl \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel \
    -r $params.REF_BASE/genome.fa \
    -t $tumor \
    -n $normal \
    -s $params.REF_BASE/pindel/simpleRepeats.bed.gz \
    -u $params.REF_BASE/pindel/pindel_np.gff3.gz \
    -f $params.REF_BASE/pindel/${params.PROTOCOL}_Rules.lst \
    -g $params.REF_BASE/vagrent/codingexon_regions.indel.bed.gz \
    -st $params.PROTOCOL \
    -as $params.ASSEMBLY \
    -sp $params.SPECIES \
    -e $params.CONTIG_EXCLUDE \
    -b $params.REF_BASE/pindel/HiDepth.bed.gz \
    -c $task.cpus \
    -sf $params.REF_BASE/pindel/softRules.lst"
}

process caveman_contam {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(normalBed), val(tumorBed)
    val NORM_CONTAM
    

    output:
 	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names


    script:
    """
	echo "CWD=\${PWD}" 1<> $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
    -r $params.REF_BASE/genome.fa.fai \
    -ig $params.REF_BASE/caveman/HiDepth.tsv \
    -b $params.REF_BASE/caveman/flagging \
    -ab $params.REF_BASE/vagrent \
    -u $params.REF_BASE/caveman \
    -s $params.SPECIES \
    -sa $params.ASSEMBLY \
    -t $task.cpus \
    -st $params.PROTOCOL \
    -tc $tumorBed \
    -nc $normalBed \
    -td 5 -nd 2 \
    -tb $tumor \
    -nb $normal \
    -c $params.SNVFLAG \
    -f $params.REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
    -e $params.CAVESPLIT \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
    -x $params.CONTIG_EXCLUDE \
    -k $NORM_CONTAM \
    -no-flagging -noclean
    """
}

process index_pindel_germline {
	errorStrategy 'finish'
	
    input:
    val GERMLINE_BED 

    output:
    val "${GERMLINE_BED}.gz"
    
    script:
    """
    if [ -f $GERMLINE_BED ]; then
    # need to sort and index pindel germline
    sort -k1,1 -k2,2n -k3,3n $GERMLINE_BED | bgzip -c > ${GERMLINE_BED}.gz
    tabix -p bed ${GERMLINE_BED}.gz
    rm -f $GERMLINE_BED
    fi
    """
}

process add_BRASS {
	errorStrategy 'finish'
	
    debug true
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(MT_sampleStats), val(CN)

    output:
	tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names


    script:
    "brass.pl -j 4 -k 4 -c $task.cpus \
    -d $params.REF_BASE/brass/HiDepth.bed.gz \
    -f $params.REF_BASE/brass/brass_np.groups.gz \
    -g $params.REF_BASE/genome.fa \
    -s $params.SPECIES -as $params.ASSEMBLY -pr $params.PROTOCOL -pl ILLUMINA \
    -g_cache $params.REF_BASE/vagrent/vagrent.cache.gz \
    -vi $params.REF_BASE/brass/viral.genomic.fa.2bit \
    -mi $params.REF_BASE/brass/all_ncbi_bacteria \
    -b $params.REF_BASE/brass/500bp_windows.gc.bed.gz \
    -ct $params.REF_BASE/brass/CentTelo.tsv \
    -cb $params.REF_BASE/brass/cytoband.txt \
    -t $tumor \
    -n $normal \
    -ss $MT_sampleStats \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/brass"
}

process pindel_annotation {
	errorStrategy 'finish'

    input:
    tuple val(NAME_MT), val(NAME_WT)


    script:
    """
    rm -f $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.annot.vcf.gz*

    AnnotateVcf.pl -t -c $params.REF_BASE/vagrent/vagrent.cache.gz \
    -i $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.flagged.vcf.gz \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/pindel/${NAME_MT}_vs_${NAME_WT}.annot.vcf
    """
}

process caveman_flag {
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)
    tuple val(normalBed), val(tumorBed)
    val NORM_CONTAM
    val GERMLINE_BED

    output:
    tuple val(tumor), val(normal), emit: files
	tuple val(NAME_MT), val(NAME_WT), emit: sample_names
    val "$params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.flagged.muts.vcf.gz", emit: vcf

    script:
    """
	echo "CWD=\${PWD}" 1<> $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/tmpCaveman/caveman.cfg.ini
    caveman.pl \
    -r $params.REF_BASE/genome.fa.fai \
    -ig $params.REF_BASE/caveman/HiDepth.tsv \
    -b $params.REF_BASE/caveman/flagging \
    -ab $params.REF_BASE/vagrent \
    -u $params.REF_BASE/caveman \
    -s $params.SPECIES \
    -sa $params.ASSEMBLY \
    -t $task.cpus \
    -st $params.PROTOCOL \
    -tc $tumorBed \
    -nc $normalBed \
    -td 5 -nd 2 \
    -tb $tumor \
    -nb $normal \
    -c $params.SNVFLAG \
    -f $params.REF_BASE/caveman/flagging/flag.to.vcf.convert.ini \
    -e $params.CAVESPLIT \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman \
    -x $params.CONTIG_EXCLUDE \
    -k $NORM_CONTAM \
    -in $GERMLINE_BED \
    -p flag
    """
}
process caveman_annotation {
	errorStrategy 'finish'
	
    input:

    tuple val(NAME_MT), val(NAME_WT)
    val caveman_flag

    script:
    """
    rm -f $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf.gz*

    AnnotateVcf.pl -t -c $params.REF_BASE/vagrent/vagrent.cache.gz \
    -i $caveman_flag \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}_vs_${NAME_WT}/caveman/${NAME_MT}_vs_${NAME_WT}.annot.muts.vcf
    """
}
process verifyBamHomChkTumor{
	errorStrategy 'finish'
	
    input:
    tuple val(tumor), val(normal) 
    tuple val(NAME_MT), val(NAME_WT)
    val ASCAT_CN

    script:
    "verifyBamHomChk.pl -d 25 \
    -o $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}/contamination \
    -b $tumor \
    -t $task.cpus \
    -a $ASCAT_CN \
    -j $params.OUTPUT_DIR/${params.PROTOCOL}_${NAME_MT}/contamination/result.json \
    -s $params.REF_BASE/verifyBamID_snps.vcf.gz"
}



workflow {
	input = Channel.fromPath(params.sampleFile).splitCsv()
	getNames(input)
	CN = CN_init()
	cavemanSetup(getNames.out.files, getNames.out.sample_names, CN)
	cavemanSplit(cavemanSetup.out.files, cavemanSetup.out.sample_names)
	ascat(getNames.out.files, getNames.out.sample_names)
	compareBamGenotypes(getNames.out.files, getNames.out.sample_names)
	prepAscat(getNames.out.sample_names, ascat.out.CN, ascat.out.sample_stats)
	caveman_contam(cavemanSplit.out.files, cavemanSplit.out.sample_names, prepAscat.out.CN_files, prepAscat.out.contam)
	brass_input(getNames.out.files, getNames.out.sample_names)
	brass_cover(brass_input.out.files, brass_input.out.sample_names)
	add_BRASS(brass_cover.out.files, brass_cover.out.sample_names, ascat.out.sample_stats)
	pindel(getNames.out.files, getNames.out.sample_names)
	index_pindel_germline(pindel.out.germline)
	pindel_annotation(pindel.out.sample_names)
	caveman_flag(caveman_contam.out.files, caveman_contam.out.sample_names,  prepAscat.out.CN_files, prepAscat.out.contam, index_pindel_germline)
	caveman_annotation(caveman_flag.out.sample_names, caveman_flag.out.vcf)
	verifyBamHomChkTumor(getNames.out.files, getNames.out.sample_names, ascat.out.CN)
	verifyBamHomChkNormal(getNames.out.files, getNames.out.sample_names)
	

	
}
