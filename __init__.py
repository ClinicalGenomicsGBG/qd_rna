from logging import LoggerAdapter
from typing import Optional
from cellophane import modules, data, cfg
from smtplib import SMTP
from email.message import EmailMessage
from mimetypes import guess_type
from pathlib import Path
from jinja2 import Environment


def _send_mail(
    *,
    from_addr: str,
    to_addr: list[str] | str,
    subject: str,
    body: str,
    host: str,
    port: int,
    tls: bool = False,
    cc_addr: Optional[list[str] | str] = None,
    user: Optional[str] = None,
    password: Optional[str] = None,
    attachments: Optional[list[Path]] = None,
    **_,
) -> None:
    conn = SMTP(host, port)
    if tls:
        conn.starttls()
    if user and password:
        conn.login(user, password)
    msg = EmailMessage()
    msg.set_content(body)
    msg["Subject"] = subject
    msg["From"] = from_addr
    msg["To"] = ", ".join(to_addr) if isinstance(to_addr, list) else to_addr
    if cc_addr is not None:
        msg["Cc"] = ", ".join(cc_addr) if isinstance(cc_addr, list) else cc_addr

    for attachment in attachments or []:
        ctype, encoding = guess_type(attachment)
        if ctype is None or encoding is not None:
            ctype = "application/octet-stream"
        maintype, subtype = ctype.split("/", 1)
        with open(attachment, "rb") as fp:
            msg.add_attachment(
                fp.read(),
                maintype=maintype,
                subtype=subtype,
                filename=Path(attachment).name,
            )

    conn.send_message(msg)
    conn.quit()


def _render_mail(subject, body, **kwargs):
    subject = Environment().from_string(subject).render(**kwargs)
    body = Environment().from_string(body).render(**kwargs)

    body_template = Environment().from_string(body)
    subject_template = Environment().from_string(subject)

    subject = subject_template.render(**kwargs)
    body = body_template.render(**kwargs)
    return subject, body


@modules.pre_hook(label="Send start mail", after="all")
def start_mail(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
):
    if not config.mail.skip:
        logger.debug(f"Sending start mail to {config.mail.start.to_addr}")
        subject, body = _render_mail(
            **config.mail,
            **config.mail.start,
            analysis=config.analysis,
            samples=samples,
        )

        cc_addr = config.mail.start.cc_addr if "cc_addr" in config.mail.start else None
        _send_mail(
            **config.mail.smtp,
            body=body,
            subject=subject,
            to_addr=config.mail.start.to_addr,
            from_addr=config.mail.start.from_addr,
            cc_addr=cc_addr,
        )


@modules.post_hook(label="Send end mail", condition="always", after="all")
def end_mail(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
):
    if not config.mail.skip:
        logger.debug(f"Sending end mail to {config.mail.end.to_addr}")
        subject, body = _render_mail(
            **config.mail,
            **config.mail.end,
            analysis=config.analysis,
            samples=samples,
        )
        cc_addr = config.mail.end.cc_addr if "cc_addr" in config.mail.end else None
        _send_mail(
            **config.mail.smtp,
            body=body,
            subject=subject,
            to_addr=config.mail.end.to_addr,
            from_addr=config.mail.end.from_addr,
            cc_addr=cc_addr,
        )
