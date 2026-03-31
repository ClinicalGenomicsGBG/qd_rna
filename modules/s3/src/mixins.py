"""S3 backup mixin for samples."""

from attrs import Attribute, define, field
from attrs.setters import convert, validate
from cellophane import Sample


@define(slots=False)
class S3Sample(Sample):
    """Sample with S3 backup."""

    s3_remote_keys: list[str] | None = field(
        default=None,
        kw_only=True,
        converter=lambda value: None if value is None else [str(v) for v in value],
        on_setattr=convert,
    )

    s3_bucket: str | None = field(
        default=None,
        kw_only=True,
        on_setattr=validate
    )

    s3_endpoint: str | None = field(
        default=None,
        kw_only=True,
        on_setattr=validate
    )

    @s3_remote_keys.validator
    def _validate_s3_remote_keys(
        self,
        attribute: Attribute,
        value: list[str] | None,
    ) -> None:
        if value and len(value) != len(self.files):
            raise ValueError(
                f"Length mismatch between {attribute.name} and files: "
                f"{len(value)} != {len(self.files)}"
            )


@Sample.merge.register("s3_remote_keys")
def merge_s3_remote_keys(this, that) -> list[str] | None:
    """Merge S3 remote keys."""
    # Use dict.fromkeys to preserve order while removing duplicates
    return None if (this or that) is None else [*(dict.fromkeys(this or []) | dict.fromkeys(that or []))]