"""Module for fetching files from HCP."""

from pathlib import Path
from uuid import UUID, uuid4
from functools import partial

from cellophane import cfg, executors, Config
from mpire.async_result import AsyncResult
from modules.mail import send_mail

_ROOT = Path(__file__).parent.parent


def nextflow(
    main: Path,
    *args,
    config: cfg.Config,
    executor: executors.Executor,
    workdir: Path,
    env: dict[str, str] | None = None,
    nxf_config: Path | None = None,
    nxf_work: Path | None = None,
    nxf_profile: str | None = None,
    ansi_log: bool = False,
    resume: bool = False,
    name: str = "nextflow",
    check: bool = True,
    **kwargs,
) -> tuple[AsyncResult, UUID]:
    """Submit a Nextflow job to SGE."""

    uuid = uuid4()
    _nxf_log = config.logdir / "nextflow" / f"{name}.{uuid.hex}.log"
    _nxf_config = nxf_config or config.nextflow.get("config")
    _nxf_work = nxf_work or config.nextflow.get("workdir") or workdir / "nxf_work"
    _nxf_launch = workdir / "nxf_launch"
    _nxf_profile = nxf_profile or config.nextflow.get("profile")

    _nxf_log.parent.mkdir(parents=True, exist_ok=True)
    _nxf_work.mkdir(parents=True, exist_ok=True)

    result, uuid = executor.submit(
        str(_ROOT / "scripts" / "nextflow.sh"),
        f"-log {_nxf_log}",
        (f"-config {_nxf_config}" if _nxf_config else ""),
        f"run {main}",
        "-ansi-log false" if not ansi_log or config.nextflow.ansi_log else "",
        f"-work-dir {_nxf_work}",
        "-resume" if resume else "",
        f"-with-report {config.logdir / 'nextflow' / f'{name}.{uuid.hex}.report.html'}",
        (f"-profile {_nxf_profile}" if _nxf_profile else ""),
        *args,
        env={
            "_NXF_INIT": config.nextflow.init,
            **config.nextflow.env,
            **(env or {}),
        },
        workdir=_nxf_launch,
        uuid=uuid,
        name=name,
        cpus=config.nextflow.threads,
        error_callback=partial(
            _nextflow_error_callback,
            uuid=uuid,
            name=name,
            workdir=workdir,
            config=config,
        ),
        **kwargs,
    )

    if check:
        result.get()

    return result, uuid

def send_crash_mail(
    config: Config,
    workdir: Path,
    uuid: UUID,
    name: str,
):
    """Send email if one runner fails.
    
    The function is intented to use as an error_callback,
    and the email should be sent while the pipeline is continuing to run for the remaining samples.
    """
    if not config.mail.send:
        return

    subject = "QD-RNA: Crashing sample"
    body = "\n".join([
    f"<p>QD-RNA has crashed for sample {workdir.name} in UUID {uuid} ({name}).</p>",
    "<p>The pipeline will continue to run if there are other samples in the same run,",
    f"but the results from the above stated tool for {workdir.name} will most likely be broken.</p>",
    f"<p>The workdir was set to the following path: <code>{workdir}</code>.</p>",
    "<p>Please investigate the crash and restart the pipeline for the failed sample if necessary.</p>"
    ])
    
    to = config.mail.to_addr[0] # Do not send crash mail to everyone, only to the first recipient

    send_mail(
        **config.mail.smtp,
        body=body,
        subject=subject,
        to_addr=to,
        from_addr=config.mail.from_addr,
        cc_addr=config.mail.get("cc_addr")
    )
    
def _nextflow_error_callback(
    exception: Exception,
    uuid: UUID,
    name: str,
    workdir: Path,
    config: Config,
) -> None:
    """Error callback for Nextflow jobs."""
    send_crash_mail(config, workdir, uuid, name)
