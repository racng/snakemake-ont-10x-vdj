#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["savefig.transparent"] = True
plt.rcParams['figure.facecolor'] = (1,1,1,0)
plt.rcParams['axes.facecolor'] = (1,1,1,0)
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["savefig.format"] = 'png'
plt.rcParams["savefig.bbox"] = 'tight'
plt.rcParams["axes.grid"] = False


def read_file(path):
    if path[-3:] == 'tsv':
        return pd.read_table(path)
    elif path[-3:] == 'csv':
        return pd.read_csv(path)

# Calculation method inspired by https://github.com/slowkow/saturation
def sample_saturation(nreads, prob, seed=123):
    # Create a Generator instance
    rng = np.random.default_rng(seed=seed)
    samples = nreads.apply(lambda x:rng.binomial(x, prob)) 
    n_deduped_reads = (samples > 0).sum()
    n_total_reads = samples.sum()
    return (n_total_reads, 1 - n_deduped_reads / n_total_reads)

def saturation_curve(path, colname, min_reads=1000 ,seed=123, label=None,
    save=None):
    df = read_file(path)
    order = np.log10(min_reads/df[colname].sum())
    sampling_freq = np.logspace(order, 0, num=20)
    xy = [sample_saturation(df[colname], p, seed=seed) for p in sampling_freq]
    # df_sat = pd.DataFrame.from_records(xy, 
    #     columns=['Total reads', 'Sequencing saturation'])
    pct = round(xy[-1][1]*100, 1)
    print(f"{pct}%")
    fig = plt.figure()
    plt.plot(*zip(*xy), label=label)
    plt.xscale("log")
    plt.xlabel('Total reads')
    plt.ylabel('Sequencing saturation')
    if save:
        fig.savefig(save)
    plt.close(fig)
    return 


# saturation_curve(
#     "/users/rng/proj/snakemake-ont-10x-vdj/results/10x-ont-cdna/20250425_1_t/airr.tsv", 
#     "consensus_count",1000, label='Nanopore')

# saturation_curve(
#     "/users/rng/proj/tlc/hd_data/10x/cellranger-9.0.0/hg38-t-tropic-virus/20250429_1/outs/multi/vdj_t/all_contig_annotations.csv", 
#     "reads",1000, label='Illumina')

# plt.legend()
# %%
saturation_curve(
    snakemake.input[0],
    snakemake.params['colname'],
    min_reads=1000, 
    seed=0, 
    label=None,
    save=snakemake.output[0]
)