"""Mail module for cellphane."""

from .src.hooks import end_mail, start_mail
from .src.mixins import MailSample, MailSamples
from .src.util import render_mail, resolve_attachments, send_mail

__all__ = [
    "MailSample",
    "MailSamples",
    "end_mail",
    "render_mail",
    "resolve_attachments",
    "send_mail",
    "start_mail",
]
