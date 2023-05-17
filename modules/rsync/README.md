# Output module for Cellophane

Module for copying output to a specified location using rsync. Uses SGE to run multiple rsync jobs in parallel. Uses the data.Output object specified in each sample.

## Configuration

Option                       | Type | Required | Default | Description
-----------------------------|------|----------|---------|-------------
`rsync.skip`                 | bool |          | false   | Don't rsync results
`rsync.base`                 | path | x        |         | Base of directory where results will be stored
`rsync.overwrite`            | bool |          | false   | Overwrite existing results
`rsync.large_file_threshold` | str  |          | "100M"  | Files larger than this will be copied in a separate job (eg. 100M, 1 GB)
`rsync.sge_queue`            | str  | x        |         | SGE queue for copying
`rsync.sge_pe`               | str  | x        |         | SGE parallel environment for copying
`rsync.sge_slots`            | int  |          | 1       | SGE slots for copying

## Hooks

Name           | When  | Condition | Description
---------------|-------|-----------|-------------
`rsync_output` | Post  | Complete  | Rsync results to output directory
