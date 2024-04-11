"""Mixins for samples with mail attachments"""

from pathlib import Path
from typing import Sequence

from attrs import define, field
from cellophane import Sample, Samples


@define(slots=False)
class MailSample(Sample):
    """A sample with mail attachments"""

    mail_attachments: set[Path] = field(factory=set)


@define(slots=False)
class MailSamples(Samples[MailSample]):
    """A collection of samples with mail attachments"""

    _mail_attachments: set[Path] = field(factory=set, init=False)

    @property
    def mail_attachments(self) -> set[Path]:
        """Get the mail attachments of all samples and the collection itself"""
        return {a for s in self for a in s.mail_attachments} | self._mail_attachments

    @mail_attachments.setter
    def mail_attachments(self, value: Sequence[Path]) -> None:
        self._mail_attachments = {*value}


@Sample.merge.register("mail_attachments")
@Samples.merge.register("_mail_attachments")
def _(this, that):
    return this | that
