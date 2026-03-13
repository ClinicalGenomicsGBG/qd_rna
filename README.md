# qd_rna

**qd_rna** is a [cellophane](https://github.com/ClinicalGenomicsGBG/cellophane) based wrapper for analyzing Whole Transcriptome Sequencing (WTS/RNAseq) data, used at Clinical Genomics Göteborg (CGG).

It runs three main processes:
- [nfcore/rnafusion](https://nf-co.re/rnafusion) for detecting gene fusions with Arriba, FusionCatcher and STAR-Fusion. Genome version GRCh38.
- [nfcore/rnaseq](https://nf-co.re/rnaseq) for gene expression data processing and quantification. Genome version GRCh38.
- Subsampling and STAR mapping with hg19 for Qlucore-based tumor classification. 
  - For now the Qlucore interpretation is performed externally, by the clinical geneticists

By default it downsamples the input if the number of reads exceeds 1M (see `properties.subsample.target` in [schema.yaml](schema.yaml))

## Usage

### Minimal manual run

A minimal manual run can be performed, where you specify the fastq file paths in a samples_file. If you do not include slims credentials in [configs/qd_rna.manual.yaml](configs/qd_rna.manual.yaml), the pipeline will skip processes which require slims access. The output will be in `/clinical/data/qdrna/manual` unless specified otherwise in the config.

```bash
# Clone the repo
git clone https://github.com/ClinicalGenomicsGBG/qd_rna.git -b main qd_rna
cd qd_rna

# Create and activate conda environment
module load micromamba
micromamba create -f qdrna_env.yaml --prefix ../qd_rna_env
micromamba activate ../qd_rna_env

# Write or copy samples file, adjust as needed
rsync /clinical/data/qdrna/test_data/samples.test.yaml .

# Start the pipeline in manual mode
python -m qd_rna --config_file configs/qd_rna.manual.yaml --samples_file ./samples.test.yaml
```

#### Samples file format

You can and should specify different things in the samples file depending on your needs.

```bash
- id: "sub_55778845"
  run: "240430_AHGTK7DRX3"
  reads: 1102526
  files:
    - /clinical/data/qdrna/test_data/sub_55778845_240110_AHWV7HDRX3_S1_R1_001.fastq.gz
    - /clinical/data/qdrna/test_data/sub_55778845_240110_AHWV7HDRX3_S1_R2_001.fastq.gz
```

- `id`: [required] Name of the sample to run. It doesn't need to match the name in the fastq name
- `run`: [required] Name of the run. It doesn't need to match the name in the fastq name
- `reads`: [optional] Used by the subsample hook to know which fraction to subsample to. If no reads are specified, QDRNA will subsample the data to a fixed number of reads, which is slower. (see `properties.subsample.target` in [schema.yaml](schema.yaml))
- `files` [optional] Path to the fastq files. **Make sure R1 is the first file, and R2 the second**. If you use only a single fastq, it will be treated as single-end data.

- If your sample is in slims and you added slims credentials, you can specify the `id` and `run` and it will find the rest for you.
- Note: If there are any samples with an identical `id`, it will download and merge them. This also happens if you specify the fastq files in the samples file. To avoid this, change the `id` in the samples file to an `id` that is not in slims.

#### Restarting a run

To restart a run, locate the workdir (in `/clinical/data/qdrna/manual/work` if you used the above config) and copy the directory name (e.g. `260419-154331`). Restart qd_rna with the `--tag` parameter (`--tag 260419-154331`). If you include `clean: true` in the config, the workdir will be removed once qd_rna finishes.

If your fastqs were subsampled in the initial run, you will need to point to those subsampled fastqs in the samples file (e.g. in `260419-154331/subsample/123456-mRNA_R[1|2]_001.subsampled.fq.gz`). You will also need to adjust the `reads` to below the cutoff, or pass `--subsample_target 0` to skip subsampling.


### Routine run

For running the routine pipeline, the process is as above but you specify `--config_file` [configs/qd_rna.routine.yaml](configs/qd_rna.routine.yaml). Additionally, you need slims credentials and email addresses. You do not need samples_files, instead it will start running on any slims samples with `qdrna_reverse` in the raw samplesheet description unless they have a bioinformatics object with `Secondary Analysis State = complete|error|running` or are more than a year old (see `slims.novel.criteria` in [configs/qd_rna.base.yaml](configs/qd_rna.base.yaml).

For the routine run, the results will be on webstore. The working directory is in `/clinical/data/qdrna/cron/work`, but will be cleaned upon successful completion.

If you have your own copy of QDRNA, make sure you also have the `emails.yaml` and the `slims_credentials.yaml` files in the QDRNA `config` directory or edit these parameters to the [configs/qd_rna.routine.yaml](configs/qd_rna.routine.yaml). Prefilled configs are available at `/clinical/exec/qdrna/qdrna/configs/`.

```bash
# Copy the email addresses and slims credentials
rsync /clinical/exec/qdrna/qdrna/configs/emails.yaml configs/
rsync /clinical/exec/qdrna/qdrna/configs/slims_credentials.yaml configs/

# Start the pipeline in routine mode
python -m qd_rna --config_file configs/qd_rna.routine.yaml
```

### Testing

There's a test dataset available:
- `samples file`: `/clinical/exec/qdrna/qdrna/configs/samples/samples.testing.yaml`
- `fastq files`: `/clinical/data/qdrna/test_data/`

To test QDRNA, follow the steps above in [Minimal manual run](#minimal-manual-run)

## Workflow details

The main code for running the pipeline is distributed over the cellophane modules. In the modules, the functions to run are decorated with either `@pre_hook`, `@runner` or `@post_hook`. [cellophane](https://github.com/ClinicalGenomicsGBG/cellophane), the backbone of `qd_rna`, will discover these and resolve in which order they should be run. The `samples` contain the information of what is run through these hooks and runners. The `runners` are decorated with `@outputs`. These outputs are expected to be present after runner completion and will be further processed by the post_hooks. For more information on how cellophane works, see the [cellophane usage](https://github.com/ClinicalGenomicsGBG/cellophane/blob/main/USAGE.md).

The following steps are typically performed in the pipeline (with relevant `modules`):

--- pre_hooks ---
- Fetch nf-core pipelines -- [modules/common.py](modules/common.py)
- Fetch information from slims -- [modules/slims/src/hooks.py](modules/slims/src/hooks.py)
- Find additional samples with the same id -- [modules/common.py](modules/common.py)
- Collect the necessary fastqs from long-term storage -- [modules/hcp/src/hooks.py](modules/hcp/src/hooks.py)
- Unpack the compressed fastqs -- [modules/unpack/src/hooks.py](modules/unpack/src/hooks.py)
- Start email -- [modules/slims/src/hooks.py](modules/slims/src/hooks.py)
- Subsample -- [modules/common.py](modules/common.py)

--- runners ---
- nfcore/rnaseq -- [modules/rnaseq.py](modules/rnaseq.py)
- nfcore/rnafusion -- [modules/rnafusion.py](modules/rnafusion.py)
- STAR mapping for Qlucore -- [modules/qlucore.py](modules/qlucore.py)

--- post_hooks ---
- Copy declared outputs to results directory -- [modules/rsync/src/hooks.py](modules/rsync/src/hooks.py)
- Update slims objects -- [modules/slims/src/hooks.py](modules/slims/src/hooks.py)
- End email -- [modules/slims/src/hooks.py](modules/slims/src/hooks.py)
- Cleanup workdir -- in [cellophane/src/cellophane/cleanup](https://github.com/ClinicalGenomicsGBG/cellophane/tree/main/src/cellophane/cleanup)

## Updates/Development

We currently have the code for all the modules tracked in the qd_rna repo. This means we can make changes as needed to the modules in the repository. At the same time, both [cellophane](https://github.com/ClinicalGenomicsGBG/cellophane) and its [cellophane_modules](https://github.com/ClinicalGenomicsGBG/cellophane_modules) are under active development. The used `cellophane` version is specified in `qd_rna_env.yaml`. We can update this to get the latest features when needed, but may need to update the `modules` as well when there are breaking changes.

In order to update the modules which are available from [cellophane](https://github.com/ClinicalGenomicsGBG/cellophane), activate an environment with the desired `cellophane` version and in the repo root run `cellophane module update`.

When updating the modules, any of the changes to the module in the `qd_rna` repo will be lost. Therefore, if you make changes to a module which you believe are of general interest, consider making a pull request to the [cellophane_modules](https://github.com/ClinicalGenomicsGBG/cellophane_modules). Otherwise, one will have to pull the commits into the updated module after updating.

### Validation

After major updates to the pipeline, we want to validate that the results are as expected. We will try to set up a routine for checking expected variants from a given sample.
