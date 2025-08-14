# ----------------------------------------------------- #
# READ PRE-PROCESSING RULES                             #
# ----------------------------------------------------- #

# nucleotide sequences of 10x library components
R1="GATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
CELLBC="N{16}"
UMI="N{10}"
TSO="TTTCTTATATGGGG"

def get_adapter(wildcards):
    return (R1+CELLBC+UMI+TSO)

def get_adapter_str(wildcards):
    adapter = get_adapter(wildcards)
    length = len(adapter)
    return f"{adapter};min_overlap={length};rightmost"

# filter and orient reads with 10x barcoding structure
# -----------------------------------------------------
rule cutadapt:
    input:
        fastq=get_r1
    output:
        fastq="results/cutadapt/{sample}.fastq.gz"
    log:
        "log/cutadapt.{sample}.log"
    params:
        adapter=get_adapter_str,
        error_rate=config['cutadapt']['error_rate'],
        action=config['cutadapt']['action']
    threads:
        config['threads']['cutadapt']
    conda:
        config['conda']['cutadapt']
    shell:
        "cutadapt -g '{params.adapter}' "
        "--rc -e {params.error_rate} "
        "-o {output.fastq} "
        "--discard-untrimmed "
        "--action={params.action} "
        "--cores {threads} "
        "{input} > {log}"


    
    