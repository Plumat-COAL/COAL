"""In-process API loop (Tier C, supplementary mode, not a headline): one interpreter calls the tool's
main(argv) once per alignment. Usage: inproc_loop.py coal|p2n16 JOBS.json HASHES_OUT.json
JOBS.json = [[case, argv], ...]; writes the normalised-output hash of each call so the harness can check it."""
import contextlib
import io
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import grade  # noqa: E402

tool, jobs_path, out_path = sys.argv[1:4]
if tool == "coal":
    from coal.cli import main
else:
    from pal2nal.cli import main

tmp = Path(out_path).with_suffix(".cur.fa")
hashes = {}
for case, argv in json.load(open(jobs_path)):
    with open(tmp, "wb") as raw:
        out = io.TextIOWrapper(raw, encoding="utf-8", write_through=True)
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            try:
                rc = main(argv)
            except SystemExit as e:
                rc = e.code
        out.flush()
        out.detach()
    recs = grade.read_fasta(tmp)
    hashes[case] = grade.records_hash(recs) if (recs and not rc) else f"rc={rc}"
json.dump(hashes, open(out_path, "w"))
