# from unittest.mock import MagicMock
# from uuid import uuid4

# from pytest import fixture, mark, param, raises
# from pytest_mock import MockerFixture

# from cellophane import sge



import logging
import multiprocessing as mp
from unittest.mock import MagicMock
from uuid import uuid4

import drmaa2
from attrs import define
from grid_engine import GridEngineExecutor
from mpire import WorkerPool
from pytest import FixtureLookupError, fixture, raises
from pytest_mock import MockerFixture


@fixture(scope="function")
def job_session(request):
    try:
        _state = request.getfixturevalue("state")
    except FixtureLookupError:
        _state = drmaa2.JobState.DONE

    class JobInfoMock:
        exit_status: int = int(_state != 7)

    class JobMock(MagicMock):
        get_info: MagicMock = MagicMock(return_value=JobInfoMock())

        def get_state(self):
            return _state, None

    @define
    class JobSessionMock():
        close = MagicMock()
        destroy = MagicMock()
        run_job = MagicMock(return_value=JobMock())
        wait_all_terminated = MagicMock()

    return JobSessionMock()


class Test_GridEngineExecutor:
    def test_submit(self, job_session, tmp_path, mocker: MockerFixture):

        mocker.patch(
            "grid_engine.drmaa2.JobSession",
            return_value=job_session,
        )

        config = MagicMock(workdir=tmp_path)
        logger = logging.getLogger()

        with (
            raises(SystemExit),
            WorkerPool(1) as pool,
        ):
            executor = GridEngineExecutor(config, pool, log_queue=mp.Queue())
            executor.target(
                "SCRIPT",
                "ARG1",
                "ARG2",
                workdir=tmp_path,
                logger=logger,
                cpus=1,
                uuid=uuid4(),
                env={"ENV": "VAR"},
                name="DUMMY",
                wait=True,
            )

# class Test__cleanup:
#     @staticmethod
#     def test__cleanup():
#         _job = MagicMock()
#         _session = MagicMock()
#         _cleanup = sge._cleanup(_job, _session)
#         with raises(SystemExit) as e:
#             _cleanup()
#         _job.terminate.assert_called_once()
#         _job.wait_terminated.assert_called_once()
#         _session.close.assert_called_once()
#         _session.destroy.assert_called_once()
#         assert e.value.code == 1


# class Test_submit:
#     def _run(self):
#         sge.submit(
#             "SCRIPT",
#             "ARG1",
#             "ARG2",
#             name="NAME",
#             # FIXME: Test with minimal config
#             config=MagicMock(
#                 logdir="LOGDIR",
#             ),
#             queue="QUEUE",
#             pe="PE",
#             slots=42,
#             env={"ENV": "VAR"},
#             cwd="CWD",
#             os_env=True,
#             check=True,
#             callback=None,
#             error_callback=None,
#         )

#     @mark.parametrize(
#         "exitcode",
#         [
#             param(0, id="exitcode=0"),
#             param(1, id="exitcode=1"),
#         ],
#     )
#     def test__submit(self, mocker, exitcode):
#         mp_Process_mock = mocker.patch(
#             "cellophane.src.sge.mp.Process",
#             return_value=MagicMock(
#                 exitcode=exitcode,
#             ),
#         )

#         uuid = uuid4()
#         mocker.patch(
#             "cellophane.src.sge.uuid4",
#             return_value=uuid,
#         )

#         if exitcode == 0:
#             self._run()
#         else:
#             assert raises(RuntimeError, self._run)

#         mp_Process_mock.assert_called_once_with(
#             target=sge._run,
#             args=("SCRIPT", "ARG1", "ARG2"),
#             kwargs={
#                 "logdir": "LOGDIR",
#                 "queue": "QUEUE",
#                 "slots": 42,
#                 "pe": "PE",
#                 "name": "NAME",
#                 "uuid": uuid,
#                 "env": {"ENV": "VAR"},
#                 "cwd": "CWD",
#                 "os_env": True,
#                 "callback": None,
#                 "error_callback": None,
#             },
#         )


# class Test__run:
#     uuid = uuid4()

#     def _run(self, tmp_path, mocker, sge_mocks):
#         mocker.patch("cellophane.src.sge.time.sleep")
#         sge._run(
#             "SCRIPT",
#             "ARG1",
#             "ARG2",
#             logdir=tmp_path,
#             uuid=self.uuid,
#             queue="QUEUE",
#             pe="PE",
#             slots=42,
#             name="NAME",
#             env={"ENV": "VAR"},
#             cwd=tmp_path,
#             callback=sge_mocks.callback_mock,
#             error_callback=sge_mocks.error_callback_mock,
#         )

#     @mark.parametrize(
#         "state,code",
#         [
#             param(sge.drmaa2.JobState.DONE, 0, id="DONE"),
#             param(sge.drmaa2.JobState.FAILED, 1, id="FAILED"),
#         ],
#     )
#     def test__run(self, tmp_path, mocker: MockerFixture, sge_mocks, state, code):
#         with raises(SystemExit) as e:
#             self._run(tmp_path, mocker, sge_mocks)

#         # sourcery skip: no-conditionals-in-tests
#         if state == sge.drmaa2.JobState.DONE:
#             sge_mocks.callback_mock.assert_called_once()
#             sge_mocks.error_callback_mock.assert_not_called()
#         else:
#             sge_mocks.callback_mock.assert_not_called()
#             sge_mocks.error_callback_mock.assert_called_once()

#         sge_mocks.run_job_mock.assert_called_once_with(
#             {
#                 "remote_command": "SCRIPT",
#                 "args": ["ARG1", "ARG2"],
#                 "min_slots": 42,
#                 "implementation_specific": {
#                     "uge_jt_pe": "PE",
#                     "uge_jt_native": "-l excl=1 -S /bin/bash -notify -q QUEUE -V",
#                 },
#                 "job_name": f"NAME_{self.uuid.hex[:8]}",
#                 "job_environment": {"ENV": "VAR"},
#                 "output_path": str(tmp_path / "sge" / f"NAME.{self.uuid.hex}.out"),
#                 "error_path": str(tmp_path / "sge" / f"NAME.{self.uuid.hex}.err"),
#                 "working_directory": str(tmp_path),
#             }
#         )

#         assert e.value.code == code

#     def test__run_drmaa2_exception(self, tmp_path, mocker: MockerFixture, sge_mocks):
#         sge_mocks.run_job_mock.side_effect = sge.drmaa2.Drmaa2Exception

#         with raises(SystemExit) as e:
#             self._run(tmp_path, mocker, sge_mocks)

#         assert e.value.code == 1
#         sge_mocks.error_callback_mock.assert_called_once()

#     def test__run_keyboard_interrupt(self, tmp_path, mocker: MockerFixture, sge_mocks):
#         sge_mocks.get_state_mock.side_effect = KeyboardInterrupt

#         with raises(SystemExit) as e:
#             self._run(tmp_path, mocker, sge_mocks)

#         assert e.value.code == 1
#         sge_mocks._cleanup_mock.assert_called_once()
#         sge_mocks.error_callback_mock.assert_called_once()
