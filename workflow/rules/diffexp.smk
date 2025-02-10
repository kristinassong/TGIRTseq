rule deseq2:
    input:
        quant = rules.merge_coco_cc.output.merged_counts,
        samples = "resources/design.tsv",
        comparisons = "resources/comparisons.tsv"
    output:
        out_files = "results/deseq2/{comp}.csv"
    params:
        out_dir = "results/deseq2"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/DESeq2.R"


rule volcano_plot:
    input:
        DE_output = rules.deseq2.output.out_files,
        tpm = rules.merge_coco_cc.output.merged_tpm
    output:
        volcano = "results/deseq2/{comp}.svg",
        up_genes = "results/deseq2/{comp}_sig_DE_up.tsv",
        down_genes = "results/deseq2/{comp}_sig_DE_down.tsv"
    params:
        pval_threshold = 0.05
    conda:
        "../envs/plots.yaml"
    message:
        "Create a volcano plot using deseq2 output for {wildcards.comp}."
    script:
        "../scripts/volcano_plot.py"