import time
from pathlib import Path

import drmaa2
from attrs import define


@define
class JobInfoMock:
    state: int = drmaa2.JobState.DONE

    @property
    def exit_status(self) -> int:
        return int(self.state != 7)


@define
class JobMock:
    state: int = drmaa2.JobState.DONE
    id: str = "DUMMY"

    def get_info(self, *args, **kwargs):
        return JobInfoMock(state=self.state)

    def get_state(self, *args, **kwargs):
        return self.state, None

    def terminate(self, *args, **kwargs):
        del args, kwargs

    def wait_terminated(self, *args, **kwargs):
        del args, kwargs


@define
class JobSessionMock:
    state: int = drmaa2.JobState.DONE
    delay: int = 0

    def close(self, *args, **kwargs):
        del args, kwargs  # Unusedc

    def destroy(self, *args, **kwargs):
        del args, kwargs  # Unused

    def run_job(self, *args, **kwargs):
        del args, kwargs  # Unused
        return JobMock(state=self.state)

    def wait_all_terminated(self, *args, **kwargs):
        time.sleep(self.delay)
        del args, kwargs  # Unused
