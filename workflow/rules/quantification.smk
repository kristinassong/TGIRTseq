rule coco_cc:
    input:
        gtf = rules.coco_ca.output.corrected_gtf,
        bam = rules.star_align.output.bam
    output:
        counts = "results/coco_cc/{sample}.tsv"
    threads:
        32
    params:
        coco_path = "resources/coco/bin"
    log:
        "results/logs/coco/coco_cc_{sample}.log"
    conda:
        "../envs/coco.yaml"
    message:
        "Quantify the number of counts, counts per million (CPM) and transcript per million (TPM) for each gene using CoCo correct_count (cc)."
    shell:
        "python {params.coco_path}/coco.py cc "
        "--countType both "
        "--thread {threads} "
        "--strand 1 "
        "--paired "
        "{input.gtf} "
        "{input.bam} "
        "{output.counts} "
        "&> {log}"


rule merge_coco_cc:
    input:
        counts = expand(rules.coco_cc.output.counts, sample=SAMPLES) 
    output:
        merged_counts = "results/coco_cc/merged_counts_{comp}.tsv",
        merged_cpm = "results/coco_cc/merged_cpm_{comp}.tsv",
        merged_tpm = "results/coco_cc/merged_tpm_{comp}.tsv"
    conda:
        "../envs/coco.yaml"
    message:
        "Merge CoCo correct count outputs into one file."
    script:
        "../scripts/merge_coco_cc_output.py"


rule coco_cb:
    input:
        bam = rules.star_align.output.bam,
        chrNameLength = rules.star_index.output
    output:
        unsorted_bedgraph = "results/coco_cb/{sample}_unsorted.bedgraph"
    params:
        pb = rules.install_pairedBamToBed12.params,
        coco = "resources/coco/bin"
    conda:
        "../envs/coco.yaml"
    threads:
        28
    message:
        "Create bedgraphs from BAM."
    shell:
        "export PATH=$PWD/{params.pb}:$PATH && "
        "python {params.coco}/coco.py cb "
        "-u " # UCSC compatible (adds a chr and track info)
        "-t {threads} "
        "-c 2500000 " # Chunk size, default value
        "{input.bam} "
        "{output.unsorted_bedgraph} "
        "{input.chrNameLength}"


rule sort_bg:
    input:
        unsorted_bedgraph = rules.coco_cb.output.unsorted_bedgraph
    output:
        sorted_bedgraph = "results/coco_cb/{sample}.bedgraph"
    message:
        "Sort bedgraphs and reformat chrM to chrMT."
    shell:
        "sort -k1,1 -k2,2n {input.unsorted_bedgraph} "
        "| sed 's/chrM/chrMT/g' > {output.sorted_bedgraph}"