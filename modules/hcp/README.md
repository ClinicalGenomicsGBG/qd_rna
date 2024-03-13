# HCP module for Cellophane

Fetch sample files from HCP backup (using NGPIris). If the backup attribute is set on a `data.Sample` this will be used as the remote keys. Otherwise, the local filenames will be used to search for matching files on the HCP.

## Configuration

Option             | Type      | Required | Default | Description
-------------------|-----------|----------|---------|-------------
`hcp.credentials`  | str       | x        |         | Path to iris credentials file
`hcp.fastq_temp`   | str       | x        |         | Path where fastqs will be stored
`hcp.parallel`     | int       |          | 4       | Number of parallel downloads from HCP

## Hooks

Name        | When | Condition | Description
------------|------|-----------|-------------
`hcp_fetch` | Pre  |           | Fetch sample files from HCP 

## Mixins

`HCPSample`

```
backup: list[str]
```
List of remote keys for files
