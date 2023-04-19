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
    cc_addr: list[str] | str,
    subject: str,
    body: str,
    host: str,
    port: int,
    tls: bool = False,
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


def _hook_proto(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    **_,
):
    logger.debug(f"Sending start mail to {config.mail.start.to_addrs}")
    subject, body = _render_mail(**config.mail, **config.mail.start, samples=samples)
    _send_mail(**config.mail, body=body, subject=subject)


start_mail = modules.pre_hook(label="Send start mail", after="all")(_hook_proto)
end_mail = modules.post_hook(label="Send end mail", condition="always")(_hook_proto)
