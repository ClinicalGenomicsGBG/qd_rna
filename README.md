# Unpacking module for Cellophane

Module for unpacking compressed files. Supports SPRING (.spring) and Petagene (.fasterq) compressed FASTQ via SGE.

## Configuration

Option             | Type | Required | Default | Description
-------------------|------|----------|---------|-------------
`unpack.sge_queue` | str  | x        |         | SGE queue for unpacking
`unpack.sge_pe`    | str  | x        |         | SGE parallel environment for unpacking
`unpack.sge_slots` | int  |          | 40      | SGE slots unpacking 
`unpack.parallel`  | int  |          | 100     | Maximum number of parallel unpack operations

## Hooks

Name     | When | Condition | Description
---------|------|-----------|-------------
`unpack` | Pre  |           | Extract sample files
