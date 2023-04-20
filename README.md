# Output module for Cellophane

Module for copying output to a specified location using rsync. Uses SGE to run multiple rsync jobs in parallel.

## Configuration

Option             | Type | Required | Default | Description
-------------------|------|----------|---------|-------------
`rsync.sge_queue`  | str  | x        |         | SGE queue for copying
`rsync.sge_pe`     | str  | x        |         | SGE parallel environment for copying
`rsync.sge_slots`  | int  |          | 1       | SGE slots for copying 
`rsync.parallel`   | int  |          | 10      | Maximum number of parallel copy jobs

## Hooks

Name           | When  | Condition | Description
---------------|-------|-----------|-------------
`rsync_output` | Post  | Complete  | Rsync results to output directory

## Classes

`Output`
Output dataclass for samples. Contains information about a single output file or directory as well as the destination directory and name.

```
src: list[Path]
dest_dir: Path
dest_name: Optional[str] = None
```

## Mixins

`OutputSample`

```
output: list[Output]
```
