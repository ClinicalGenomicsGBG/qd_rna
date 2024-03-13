import logging
from pathlib import Path
from typing import Any
from uuid import UUID

import drmaa2
from attrs import define, field
from cellophane import executors


@define(slots=False)
class GridEngineExecutor(executors.Executor, name="grid_engine"):
    """Executor using multiprocessing."""

    ge_jobs: dict[
        UUID,
        tuple[
            drmaa2.JobSession,
            drmaa2.Job,
        ],
    ] = field(factory=dict, init=False)

    def target(
        self,
        *args: str,
        name: str,
        uuid: UUID,
        workdir: Path,
        env: dict[str, str],
        os_env: bool = True,
        logger: logging.LoggerAdapter,
        cpus: int,
        **kwargs: Any,
    ) -> None:
        del kwargs  # Unused

        _logdir = self.config.logdir / "grid_engine" / uuid.hex
        _logdir.mkdir(exist_ok=True, parents=True)

        try:
            session = drmaa2.JobSession(f"{name}_{uuid.hex}")
            job = session.run_job(
                {
                    "remote_command": args[0],
                    "args": args[1:],
                    "min_slots": cpus,
                    "implementation_specific": {
                        "uge_jt_pe": self.config.grid_engine.pe,
                        "uge_jt_native": (
                            "-l excl=1 "
                            "-S /bin/bash "
                            f"-notify -q {self.config.grid_engine.queue} "
                            f"{'-V' if os_env else ''}"
                        ),
                    },
                    "job_name": f"{name}_{uuid.hex[:8]}",
                    "job_environment": env,
                    "output_path": str(_logdir / f"{name}.out"),
                    "error_path": str(_logdir / f"{name}.err"),
                    "working_directory": str(workdir),
                }
            )
            logger.debug(f"Grid Engine job started ({name=}, {uuid=}, {job.id=})")
            self.ge_jobs[uuid] = (session, job)
        except Exception as exception:  # pylint: disable=broad-except
            logger.error(f"Failed to submit job to Grid Engine ({name=}, {uuid=})")
            with open(
                _logdir / f"{name}.err",
                mode="w",
                encoding="utf-8",
            ) as f:
                f.write(str(exception))
            exit(1)

        session.wait_all_terminated([job])
        job_info = job.get_info()
        session.close()
        session.destroy()
        exit(job_info.exit_status)

    def terminate_hook(self, uuid: UUID, logger: logging.LoggerAdapter) -> int:
        if uuid in self.ge_jobs:
            session, job = self.ge_jobs[uuid]
            logger.warning(f"Terminating SGE job (id={job.id})")
            job.terminate()
            job.wait_terminated()
            session.close()
            session.destroy()
            return 143
