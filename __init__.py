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
    params: dict[str, str],
    **kwargs,
) -> None:
    conn = SMTP(host, port)
    if tls:
        conn.starttls()
    if user and password:
        conn.login(user, password)
    msg = EmailMessage()
    msg.set_content(body.format(**params))
    msg['Subject'] = subject
    msg['From'] = from_addr
    msg['To'] = ", ".join(to_addr) if isinstance(to_addr, list) else to_addr
    msg['Cc'] = ", ".join(cc_addr) if isinstance(cc_addr, list) else cc_addr

    if isinstance(attachments, str):
        attachments = [attachments]
    for attachment in attachments:
        ctype, encoding = guess_type(attachment)
        if ctype is None or encoding is not None:
            ctype = "application/octet-stream"
        maintype, subtype = ctype.split('/', 1)
        with open(attachment, 'rb') as fp:
            msg.add_attachment(
                fp.read(),
                maintype=maintype,
                subtype=subtype, 
                filename=Path(attachment).name
            )

    conn.send_message(msg)
    conn.quit()


@modules.pre_hook(label="Send end mail", after="all")
def start_mail(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    timestamp: str,
    **kwargs,
):
    logger.debug(f"Sending start mail to {config.mail.start.to_addrs}")
    template = Environment().from_string(config.mail.start.body)
    template.render(
        config=config,
        samples=samples,
        analysis=config.analysis,
    )
    _send_mail(**config.mail, **config.mail.start)


@modules.post_hook(label="Send end mail", condition="always")
def end_mail(
    samples: data.Samples[data.Sample],
    logger: LoggerAdapter,
    config: cfg.Config,
    timestamp: str,
    **kwargs,
):
    logger.debug(f"Sending end mail to {config.mail.end.to_addr}")
    template = Environment().from_string(config.mail.end.body)
    body = template.render(
        config=config,
        samples=samples,
        analysis=config.analysis,
    )
    _send_mail(**config.mail, **config.mail.end)