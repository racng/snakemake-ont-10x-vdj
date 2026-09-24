# ----------------------------------------------------- #
# COMMON PRE-ASSEMBLY RULES                             #
# ----------------------------------------------------- #

# copy preset to workdir 
# -----------------------------------------------------
rule copy_preset:
    input:
        yaml="workflow/presets/{preset}.yaml"
    output:
        yaml=temp("{preset}.yaml")
    shell:
        "cp {input.yaml} {output.yaml}"

# align ont reads
# -----------------------------------------------------
rule align_ont:
    input:
        fastq="results/cutadapt/{sample}.fastq.gz",
        yaml="{preset}.yaml"
    output:
        vdjca="results/{preset}/{sample}/alignments.vdjca",
        txt="results/{preset}/{sample}/align.report.txt",
        json="results/{preset}/{sample}/align.report.json"
    wildcard_constraints:
        preset="|".join(ONT_PRESETS)
    threads:
        config["threads"]["align"]
    resources:
        mem_gb=config["mem_gb"]["align"]
    log:
        "logs/{preset}.{sample}.align.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr align -p local:{wildcards.preset} -f "
        "--species hsa "
        "-t {threads} -Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input.fastq} "
        "{output.vdjca} &> {log}"

# correct barcodes/umi
# -----------------------------------------------------
# no threads parameter
rule refine:
    input:
        "results/{preset}/{sample}/alignments.vdjca"
    output:
        vdjca="results/{preset}/{sample}/alignments.refined.vdjca",
        txt="results/{preset}/{sample}/refine.report.txt",
        json="results/{preset}/{sample}/refine.report.json"
    resources:
        mem_gb=config["mem_gb"]["refine"]
    log:
        "logs/{preset}.{sample}.refine.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr refineTagsAndSort -f "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.vdjca} &> {log}"

# ----------------------------------------------------- #
# AMPLICON ASSEMBLY WORKFLOW                            #
# ----------------------------------------------------- #

# assemble directly from refined alignments
# -----------------------------------------------------
rule assemble_refined:
    input:
        "results/{preset}/{sample}/alignments.refined.vdjca"
    output:
        clna="results/{preset}/{sample}/assembled.clns",
        txt="results/{preset}/{sample}/assemble.report.txt",
        json="results/{preset}/{sample}/assemble.report.json"
    wildcard_constraints:
        preset="|".join(AMPLICON_PRESETS)
    resources:
        mem_gb=config["mem_gb"]["assemble"]
    log:
        "logs/{preset}.{sample}.assemble.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr assemble -f "
        # "--write-alignments "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.clna} &> {log}"

# ----------------------------------------------------- #
# FRAGMENTED ASSEMBLY WORKFLOW                          #
# ----------------------------------------------------- #

# overlap alignments from same molecule covering CDR3
# -----------------------------------------------------
rule assemble_partial:
    input:
        "results/{preset}/{sample}/alignments.refined.vdjca"
    output:
        vdjca="results/{preset}/{sample}/alignments.recovered.vdjca",
        txt="results/{preset}/{sample}/assemble_partial.report.txt",
        json="results/{preset}/{sample}/assemble_partial.report.json"
    wildcard_constraints:
        preset="|".join(FRAGMENTED_PRESETS)
    resources:
        mem_gb=config["mem_gb"]["assemble"]
    log:
        "logs/{preset}.{sample}.assemble_partial.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr assemblePartial -f "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.vdjca} &> {log}"

# extend edges of CDR3 in recovered alignments
# -----------------------------------------------------
rule extend_recovered:
    input:
        "results/{preset}/{sample}/alignments.recovered.vdjca"
    output:
        vdjca="results/{preset}/{sample}/alignments.extended.vdjca",
        txt="results/{preset}/{sample}/extend_recovered.report.txt",
        json="results/{preset}/{sample}/extend_recovered.report.json"
    wildcard_constraints:
        preset="|".join(FRAGMENTED_PRESETS)
    resources:
        mem_gb=config["mem_gb"]["extend"]
    log:
        "logs/{preset}.{sample}.extend_recovered.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr extend -f "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.vdjca} &> {log}"

# assemble from extended alignents
# -----------------------------------------------------
rule assemble_extended:
    input:
        "results/{preset}/{sample}/alignments.extended.vdjca"
    output:
        clna="results/{preset}/{sample}/extended_assembled.clna",
        txt="results/{preset}/{sample}/assemble_extended.report.txt",
        json="results/{preset}/{sample}/assemble_extended.report.json"
    wildcard_constraints:
        preset="|".join(FRAGMENTED_PRESETS)
    resources:
        mem_gb=config["mem_gb"]["assemble"]
    log:
        "logs/{preset}.{sample}.assemble.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr assemble -f "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.clna} &> {log}"

# assemble contigs
# -----------------------------------------------------
rule assemble_contigs:
    input:
        "results/{preset}/{sample}/extended_assembled.clna"
    output:
        clns="results/{preset}/{sample}/assembled.clns",
        txt="results/{preset}/{sample}/assemble_contigs.report.txt",
        json="results/{preset}/{sample}/assemble_contigs.report.json"
    wildcard_constraints:
        preset="|".join(FRAGMENTED_PRESETS)
    threads:
        config["threads"]["assemble_contigs"]
    resources:
        mem_gb=config["mem_gb"]["assemble"]
    log:
        "logs/{preset}.{sample}.assemble.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr assembleContigs -f "
        "-t {threads} -Xmx{resources.mem_gb}g "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.clns} &> {log}"

# ----------------------------------------------------- #
# COMMON EXPORT RULES                                   #
# ----------------------------------------------------- #

# export clones as tsv from assembled clns file
# -----------------------------------------------------
rule export_clones:
    input:
        "results/{preset}/{sample}/assembled.clns"
    output:
        "results/{preset}/{sample}/clones.tsv"
    resources:
        mem_gb=config["mem_gb"]["export_clones"]
    log:
        "logs/{preset}.{sample}.export_clones.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr exportClones -f "
        "-Xmx{resources.mem_gb}g "
        "--verbose "
        "{input} > {output} "
        "2> {log}"

# export clones in AIRR format from assembled clns file
# -----------------------------------------------------
rule export_airr:
    input:
        "results/{preset}/{sample}/assembled.clns"
    output:
        "results/{preset}/{sample}/airr.tsv"
    log:
        "logs/{preset}.{sample}.export_airr.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr exportAirr -f "
        "--verbose "
        "{input} {output} "
        "&> {log}"

# extend clones from assembled clns file
# -----------------------------------------------------
rule extend_assembled:
    input:
        "results/{preset}/{sample}/assembled.clns"
    output:
        clna="results/{preset}/{sample}/extended.clns",
        txt="results/{preset}/{sample}/extend.report.txt",
        json="results/{preset}/{sample}/extend.report.json"
    threads:
        config["threads"]["extend"]
    resources:
        mem_gb=config["mem_gb"]["extend"]
    log:
        "logs/{preset}.{sample}.extend.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr extend -f "
        "-t {threads} -Xmx{resources.mem_gb}g "
        "--v-anchor L1Begin --j-anchor FR4End "
        "--verbose "
        "--report {output.txt} "
        "--json-report {output.json} "
        "{input} "
        "{output.clna} &> {log}"

# export clones with imputed features in AIRR format from extended clns file
# -----------------------------------------------------
rule export_airr_extended:
    input:
        "results/{preset}/{sample}/extended.clns"
    output:
        "results/{preset}/{sample}/airr_extended.tsv"
    log:
        "logs/{preset}.{sample}.export_airr_extended.log"
    conda:
        config['conda']['mixcr']
    shell:
        "mixcr exportAirr -f "
        "--verbose "
        "{input} {output} "
        "&> {log}"

rule plot_saturation:
    input:
        "results/{preset}/{sample}/airr.tsv"
    output:
        "results/{preset}/{sample}/saturation.png"
    conda:
        config['conda']['plot']
    params:
        colname="consensus_count"
    script:
        "../scripts/saturation.py"