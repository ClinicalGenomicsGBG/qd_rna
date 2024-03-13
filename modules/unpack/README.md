# Unpack module for Cellophane

Module for unpacking compressed files. Supports SPRING (.spring) and Petagene (.fasterq) compressed FASTQ via SGE.

## Configuration

Option             | Type | Required | Default | Description
-------------------|------|----------|---------|-------------
`unpack.init`      | str  |          |         | Code to run before unpacking (Bash)
`unpack.exit`      | str  |          |         | Code to run after unpacking (Bash)
`unpack.threads`   | int  |          | 40      | Threads for decompression
`unpack.timeout`   | int  |          | 60      | Timeout (in seconds) to wait for unpacked file to become available

## Hooks

Name     | When | Condition | Description
---------|------|-----------|-------------
`unpack` | Pre  |           | Extract sample files
