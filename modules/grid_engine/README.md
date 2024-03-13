# Grid Engine module for Cellophane

Executor for running jobs on Grid Engine via DRMAA2.

## Configuration

Option              | Type      | Required | Default | Description
--------------------|-----------|----------|---------|-------------
`grid_engine.queue` | str       |          | all.q   | Queue to submit jobs to
`grid_engine.pe`    | str       |          | mpi     | Parallel environment to use

## Executors

Name          | Mem  | CPUs        | Description
--------------|------|-------------|-------------
`grid_engine` | No   | Yes (Slots) | Run jobs on Grid Engine
