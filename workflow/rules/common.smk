# import basic packages
import pandas as pd
from snakemake.utils import validate

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample_id": str})
    .set_index("sample_id", drop=False)
    .sort_index()
)

# assert sample ids are unique
if samples.index.duplicated().any():
    raise ValueError("Cannot use duplicated sample ids.")

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yaml")
validate(config, schema="../../config/schemas/config.schema.yaml")

# default presets for given chemistry and library type
PRESETS = {
    ('5p_v2', 'cdna'): '10x-ont-cdna',
    ('5p_v2', 'illumina'): '10x-ont-illumina'
}
# preset 10x-ont-cdna-clna doesn't work yet because extend bug 
AMPLICON_PRESETS = ['10x-ont-cdna', '10x-ont-cdna-clna']
FRAGMENTED_PRESETS = ['10x-ont-illumina']
ONT_PRESETS = ['10x-ont-cdna', '10x-ont-cdna-clna', '10x-ont-illumina']
PE_PRESETS = ['10x-sc-xcr-vdj']

# determine mixcr preset:
def map_preset(row):
    c, l = row['10x_chemistry'], row['library_type']
    return PRESETS[(c, l)]
presets = samples.apply(map_preset, axis=1)
samples['preset'] = samples['preset'].fillna(presets)

wildcard_constraints:
    sample="|".join(samples["sample_id"]),
    preset="|".join(samples["preset"])

def get_r1(wildcards):
    return samples.loc[wildcards.sample, 'read1']

def get_r1_r2(wildcards):
    reads = {
        'read1': samples.loc[wildcards.sample, 'read1'],
        'read2': samples.loc[wildcards.sample, 'read2']
    }
    return reads