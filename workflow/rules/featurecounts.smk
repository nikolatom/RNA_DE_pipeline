rule feature_counts:
    input:
        sam=sorted(set(expand("../../results/star/pe/{sample}-{unit}/{sample}_Aligned.out.sorted.bam",sample=units.sample_name, unit=units.unit_name))),
        annotation="resources/genome.gtf",
        # optional input
        # chr_names="",          # implicitly sets the -A flag
        fasta="resources/genome.fasta" # implicitly sets the -G flag,

    output:
        "../../results/featurecounts/all.featureCounts",
        "../../results/featurecounts/all.featureCounts.summary",
        "../../results/featurecounts/all.featureCounts.jcounts", 

    threads:
        38
    params:
        tmp_dir="",   # implicitly sets the --tmpDir flag
        r_path="",    # implicitly sets the --Rpath flag
        extra="-p -t exon -g gene_id -s 2 -B -C --minOverlap 60"
    log:
        "logs/featurecounts/all.log"

    conda:
        "../wrappers/executive_wrappers/subread/featurecounts/environment.yaml"
    script:
        "../wrappers/executive_wrappers/subread/featurecounts/wrapper.py"

