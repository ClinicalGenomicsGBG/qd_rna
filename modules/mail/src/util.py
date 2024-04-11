"""Utility functions for the mail module."""

from email.message import EmailMessage
from logging import LoggerAdapter
from mimetypes import guess_type
from pathlib import Path
from smtplib import SMTP
from typing import Sequence

from cellophane import Config, Samples
from jinja2 import Environment
from mistletoe import markdown


def send_mail(
    *,
    from_addr: str,
    to_addr: list[str] | str,
    subject: str,
    body: str,
    host: str,
    port: int,
    tls: bool = False,
    cc_addr: list[str] | str | None = None,
    user: str | None = None,
    password: str | None = None,
    attachments: set[Path] | None = None,
    **_,
) -> None:
    """
    Send an email with the specified parameters.

    Args:
        from_addr: The sender's email address.
        to_addr: The recipient's email address(es).
        subject: The subject of the email.
        body: The body of the email (HTML)
        host: The SMTP server host.
        port: The SMTP server port.
        tls: Whether to use TLS (default is False).
        cc_addr: The CC email address(es) (optional).
        user: The SMTP server username (optional).
        password: The SMTP server password (optional).
        attachments: Set of file paths for attachments (optional).

    Returns:
        None
    """
    conn = SMTP(host, port)
    if tls:
        conn.starttls()
    if user and password:
        conn.login(user, password)
    msg = EmailMessage()
    msg.set_content(body, subtype="html")
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


def render_mail(subject, body, **kwargs):
    """
    Render the email subject and body using provided templates and keyword arguments.

    Args:
        subject: Jinja template for the email subject.
        body: Jinja template for the email body (supports Markdown).
        **kwargs: Additional keyword arguments to render the templates.

    Returns:
        Tuple containing the rendered subject and body as HTML.
    """
    body_template = Environment().from_string(body)
    subject_template = Environment().from_string(subject)

    subject_ = subject_template.render(**kwargs)
    body_ = markdown(body_template.render(**kwargs))
    return subject_, body_


def resolve_attachments(
    attachments: Sequence[str],
    logger: LoggerAdapter,
    samples: Samples,
    config: Config,
) -> set[Path]:
    """
    Resolve the attachments based on the provided samples and configuration.

    Args:
        attachments: List of attachment paths with placeholders.
        logger: Logger adapter for logging messages.
        samples: Collection of samples to resolve attachments for.
        config: Configuration settings for resolving attachments.

    Returns:
        Set of resolved attachment paths.
    """

    attachments_: set[Path] = set()
    for sample in samples:
        attachments_ |= {
            Path(
                str(a).format(
                    sample=sample,
                    samples=samples,
                    config=config,
                )
            )
            for a in attachments
        }

    for attachment in attachments_.copy():
        if attachment.is_dir():
            logger.warning(f"Attachment {attachment} is a directory")
            attachments_.remove(attachment)
        elif not attachment.is_file():
            logger.warning(f"Attachment {attachment} is not a file")
            attachments_.remove(attachment)
        elif attachment.is_symlink():
            # Replace the attachment with its resolved path
            attachments_ ^= {attachment, attachment.resolve()}

    return attachments_
