from datetime import datetime
from types import SimpleNamespace

from cellophane import Timestamp

from modules.common import resolve_s3_upload_path


class AttrDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


def test_resolve_s3_upload_path():
    fixed_time = datetime(2027, 3, 14).timestamp()
    timestamp = Timestamp(time=fixed_time)

    config = SimpleNamespace(
        s3=SimpleNamespace(
            upload=AttrDict(
                path="s3://qdrna-results/{timestamp[%Y]}",
            )
        )
    )

    resolve_s3_upload_path(config, timestamp)

    assert timestamp["%Y"] == "2027"
    assert config.s3.upload.path == "s3://qdrna-results/2027"
