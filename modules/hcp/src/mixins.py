"""HCP mixin for samples."""

from attrs import Attribute, define, field
from attrs.setters import convert
from cellophane.src.data import Sample


@define(slots=False)
class HCPSample(Sample):
    """Sample with HCP backup."""

    hcp_remote_keys: set[str] | None = field(
        default=None,
        kw_only=True,
        converter=lambda value: None if value is None else {str(v) for v in value},
        on_setattr=convert,
    )

    @hcp_remote_keys.validator
    def _validate_hcp_remote_keys(
        self,
        attribute: Attribute,
        value: set[str] | None,
    ) -> None:
        if value and len(value) != len(self.files):
            raise ValueError(
                f"Length mismatch between {attribute.name} and files: "
                f"{len(value)} != {len(self.files)}"
            )


@Sample.merge.register("hcp_remote_keys")
def merge_hcp_remote_keys(this, that) -> set[str] | None:
    """Merge HCP remote keys."""
    return None if (this or that) is None else (this or set()) | (that or set())
