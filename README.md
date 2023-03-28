# HCP module for Cellophane

Fetch sample files from HCP backup (using NGPIris). If the backup attribute is set on a `data.Sample` this will be used as the remote keys. Otherwise, the local filenames will be used to search for matching files on the HCP.

## Configuration

type: object
properties:
  iris:
    type: object
    properties:
      credentials:
        type: path
        description: IRIS credentials file
      fastq_temp:
        type: path
        description: Temporary directory for fastq files
      parallel:
        type: integer
        description: Number of parallel HCP conenctions
        default: 4
    required:
      - credentials
      - fastq_temp


Option             | Type      | Required | Default | Description
-------------------|-----------|----------|---------|-------------
`iris.credentials` | str       | x        |         | Path to iris credentials file
`iris.fastq_temp`  | str       | x        |         | Path where fastqs will be stored
`iris.parallel`    | int       |          | 4       | Number of parallel downloads from HCP

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
