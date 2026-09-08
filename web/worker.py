"""Detached worker for one REGOLIT web run: ``python -m web.worker <run_id>``.

Started by the web server in its own process session, so it keeps going when the server reloads,
restarts or loses its client. It reads ``<runs>/<run_id>/request.json`` (written when the run was
queued), runs the model and writes ``summary.json`` (status ``done`` or ``failed``).

``python -m web.worker --finalize <run_id>`` only does the post-processing (layer-stack compaction
and summary) for a run whose model output already exists, for example a run adopted from an older
server process.
"""
import json
import sys
import time

from web.server import RUNS_DIR, execute, fail_run, finalize, load_summary  # noqa: E402


def main(argv):
    finalize_only = "--finalize" in argv
    run_id = [a for a in argv if not a.startswith("--")][0]
    workdir = RUNS_DIR / run_id
    try:
        if finalize_only:
            summary = load_summary(run_id)
            finalize(run_id, summary, summary.get("started_at", time.time()))
        else:
            request = json.loads((workdir / "request.json").read_text())
            execute(run_id, request["overrides"], request["layers"], request.get("presets"))
    except Exception as error:  # noqa: BLE001
        fail_run(run_id, "the model run failed: {}".format(str(error)[:1500]))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
