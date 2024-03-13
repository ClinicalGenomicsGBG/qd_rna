# Nextflow module for Cellophane

Module providing convenient access to functionality for launching NextFlow pipelines (particularly nf-core) from runners via SGE.

## Configuration

Option                      | Type | Required | Default  | Description
----------------------------|------|----------|----------|-------------
`nextflow.threads`          | int  |          | 2        | Threads for nextflow manager
`nextflow.config`           | str  |          |          | Nextflow config file
`nextflow.profile`          | str  |          |          | Nextflow profile
`nextflow.workdir`          | str  |          |          | Nextflow workdir
`nextflow.ansi_log`         | bool |          | false    | Enable ANSI log
`nextflow.init`             | str  |          |          | Code to run before running Nextflow (Bash)
`nextflow.env`              | dict |          | {}       | Environment variables that will be passed to the nextflow process

# Functions

```python
def nextflow(
    main: Path,
    *args,
    config: cellophane.Config,
    executor: cellophane.Executor,
    workdir: Path,
    env: dict[str, str] | None = None,
    nxf_config: Path | None = None,
    nxf_work: Path | None = None,
    nxf_profile: str | None = None,
    ansi_log: bool = False,
    resume: bool = False,
    name: str = "nextflow",
    check: bool = True,
    **kwargs
) -> tuple[AsyncResult, UUID]:
```

Launch the NextFlow pipeline specified by `main`  with `*args` as CLI arguments. `config` should be the current configuration. `executor` is the executor instance that will be used to spawn the NextFlow manager process. `workdir` is the working directory for the NextFlow manager. Environment variables can be passed using `env`. `nxf_config`, `nxf_work`, and `nxf_profile` are the NextFlow config file, workdir, and profile, respectively. `ansi_log` enables ANSI log. `resume` enables resuming a previous run. `name` is the name of the process. `check` will raise an exception if the process fails. Any `**kwargs` will be passed to the executor's `submit` method.

## Mixins

`NextflowSamples`

```python
nfcore_samplesheet(self, *_, location: str | Path, **kwargs) -> Path:
```

Write an nf-core compatible sample sheet at `location`. Any `**kwargs` will be added as columns. Values in `**kwargs` may be strings (same will be added to all samples) or a `sample.id -> value` mapping, in which case the value for each sample cn be assigned individually.
