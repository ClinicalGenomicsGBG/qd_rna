# Output module for Cellophane

Module for copying output to a specified location using rsync. Spawns one process for large files, one for small files, and one for directories.

## Configuration

Option                       | Type | Required | Default | Description
-----------------------------|------|----------|---------|-------------
`rsync.overwrite`            | bool |          | false   | Overwrite existing results
`rsync.large_file_threshold` | str  |          | 100M    | Files larger than this will be copied in a separate job (eg. 100M, 1 GB)
`rsync.timeout`              | int  |          | 30      | Timeout (in seconds) for files to be available after rsync is complete

## Hooks

Name           | When  | Condition | Description
---------------|-------|-----------|-------------
`rsync_output` | Post  | Complete  | Rsync results to output directory
