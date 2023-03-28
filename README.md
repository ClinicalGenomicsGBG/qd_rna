# SLIMS module for Cellophane

Module for unpacking compressed files. Supports SPRING (.spring) and Petagene (.fasterq) compressed FASTQ via SGE.

## Configuration

  unpack:
    type: object
    properties:
      sge_queue:
        type: string
        description: SGE queue for decompression
      sge_pe:
        type: string
        description: SGE parallel environment for decompression
      sge_slots:
        type: integer
        description: SGE slots (threads) for decompression
        default: 40
      parallel:
        type: integer
        description: Maximum number of decompression jobs to submit at once
        default: 100
    required:
      - sge_queue
      - sge_pe


Option             | Type | Required | Default | Description
-------------------|------|----------|---------|-------------
`unpack.sge_queue` | str  | x        |         | SGE queue for unpacking
`unpack.sge_pe`    | str  | x        |         | SGE parallel environment for unpacking
`unpack.sge_slots` | int  |          | 40      | SGE slots unpacking 
`unpack.parallel`  | int  |          | 100     | Maximum number of parallel unpack operations

## Hooks

Name                    | When | Condition | Description
------------------------|------|-----------|-------------
`petagene_extract`      | Pre  |           | Extract sample files
