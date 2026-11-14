"""scripts/run_tutorials.py, the tutorials-as-tests runner: how it reads the tutorials (fenced blocks, command vs
output, what is skipped), how it compares outputs, and that every tutorial still parses into something it can
check. Fast and offline; the tutorials themselves run in CI's "Tutorials as tests" job."""
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location("run_tutorials", os.path.join(ROOT, "scripts", "run_tutorials.py"))
rt = importlib.util.module_from_spec(_spec)
sys.modules["run_tutorials"] = rt  # dataclasses look the module up by name
_spec.loader.exec_module(rt)

SAMPLE = """# Sample

## 0. Set up

```
export T=/path/to/COAL/docs/tutorials
mkdir -p coal-tutorials && cd coal-tutorials
```

## 1. Run

```
coal pipeline pairs --query q.fna -o out -j 1
```

```
coal pipeline: 3,142 pairs (rbh)
coal pipeline: summary: decompress_inputs 0.0 s, align_codon_kaks 9.7 s
coal pipeline: summary: results in .../out/dnds.tsv
```

- A list item with its own block:

  ```
  coal a.aln b.cds --strict
  ```
  ```
  coal: error: X: inconsistent
  ...
  exit 1
  ```

```
$ head -c 12 x.paml
   3   1536
```

```r
library(coalr)
```

```
not the output of anything above
```

## Clean up

```
rm -rf out
```
"""


class TestParse(unittest.TestCase):
    def setUp(self):
        self.blocks = rt.parse_blocks(SAMPLE, "01")
        self.run = rt.plan(self.blocks)

    def test_blocks_headings_and_indented_fences(self):
        self.assertEqual([b.heading for b in self.blocks][:3], ["0. Set up", "1. Run", "1. Run"])
        strict = [b for b in self.blocks if "--strict" in b.text][0]
        self.assertEqual(strict.text, "coal a.aln b.cds --strict")  # dedented
        self.assertEqual(strict.outputs[0].text.splitlines()[-1], "exit 1")

    def test_command_and_output_blocks(self):
        kinds = [b.kind for b in self.blocks]
        self.assertEqual(kinds, ["command", "command", "output", "command", "output", "prompt", "other", "output",
                                 "command"])
        # the output after an R block is not attached to the shell block before it
        self.assertEqual([len(b.outputs) for b in self.run], [0, 1, 1, 0, 0])

    def test_output_lines_that_start_like_commands(self):
        for line in ("coal pipeline: 3,142 pairs", "coal: error: x", "coal dnds: 3 pair(s)", "[coal] Newton C",
                     "coal translate: 3 gene(s)", "ok                 coal-tutorial-data/x.gz", "     11 ok"):
            self.assertEqual(rt.classify(rt.Block("01", 0, 1, "", "", line)), "output", line)
        for line in ("coal dnds x.fasta -c 4", "D=../coal-tutorial-data", "python3 - <<'EOF'", "for f in a b; do",
                     "zcat x.gz | awk 1", "$ head -c 120 x"):
            self.assertIn(rt.classify(rt.Block("01", 0, 1, "", "", line)), ("command", "prompt"), line)
        self.assertEqual(rt.classify(rt.Block("01", 0, 1, "", "bash", "anything")), "command")
        self.assertEqual(rt.classify(rt.Block("01", 0, 1, "", "text", "coal x y")), "output")

    def test_prompt_block_and_script(self):
        prompt = [b for b in self.run if b.kind == "prompt"][0]
        self.assertEqual(rt.script_of(prompt, None), "head -c 12 x.paml\n")
        self.assertEqual(rt.documented_lines(prompt), ["   3   1536"])
        setup = self.run[0]
        self.assertTrue(rt.script_of(setup, None).startswith("# (set by the runner) export T="))

    def test_jobs_substitution_only_where_allowed(self):
        b = rt.Block("03", 0, 1, "", "", "coal pipeline clusters --cds $CDS -o x -j 1\ncoal dnds a -j 1 | cut -f1")
        self.assertEqual(rt.script_of(b, 4).count("-j 4"), 2)
        b.tutorial = "02"
        self.assertNotIn("-j 4", rt.script_of(b, 4))

    def test_skips(self):
        rt.apply_skips("01", self.blocks, {"no_pal2nal"})
        self.assertIn("clean-up", self.run[-1].skip)
        b = rt.Block("01", 0, 1, "3. x", "", "pal2nal.pl a b > c")
        rt.apply_skips("01", [b], {"no_pal2nal"})
        self.assertIn("PAL2NAL", b.skip)
        b.skip = ""
        rt.apply_skips("01", [b], set())
        self.assertEqual(b.skip, "")

    def test_exit_statuses(self):
        strict = [b for b in self.run if "--strict" in b.text][0]
        self.assertEqual(rt.expected_exit(strict), {1})
        self.assertTrue(rt.failures_allowed(strict))
        echo = rt.Block("01", 0, 1, "", "", 'coal a b -o c\necho "exit $?"')
        echo.outputs = [rt.Block("01", 1, 2, "", "", "coal: error: x\nexit 1")]
        self.assertEqual(rt.expected_exit(echo), {0})  # the block prints the status itself
        self.assertTrue(rt.failures_allowed(echo))
        self.assertEqual(rt.expected_exit(self.run[1]), {0})
        self.assertFalse(rt.failures_allowed(self.run[1]))


class TestCompare(unittest.TestCase):
    NORM = dict(coal_version="0.7.0", doc_version="0.6.0.dev0", scratch="/tmp/x/coal-tutorials")

    def test_normalisation(self):
        n = lambda s: rt.normalise(s, **self.NORM)
        self.assertEqual(n("coal pipeline: summary: diamond_search 6.0 s, orthologs 0.0 s"),
                         n("coal pipeline: summary: diamond_search 12.3 s, orthologs 0.1 s"))
        self.assertEqual(n("completed 0 0.6.0.dev0 2026-11-14T14:57:54Z"), n("completed 0 0.7.0 2027-01-02T03:04:05Z"))
        self.assertEqual(n("ok (already here)  coal-tutorial-data/a.gz"), n("ok                 coal-tutorial-data/a.gz"))
        self.assertEqual(n("     11 not_converged"), "11 not_converged")
        self.assertEqual(n("results in /tmp/x/coal-tutorials/out/dnds.tsv"), "results in .../out/dnds.tsv")

    def test_missing_lines_and_wildcards(self):
        doc = ["coal pipeline: 3,142 pairs (rbh)", "...", "results in .../out/dnds.tsv",
               "warning: 1,325 of the 4,074 pairs (33%) are saturated: ..."]
        actual = ["coal pipeline: 3,142 pairs (rbh)", "results in /tmp/x/coal-tutorials/out/dnds.tsv",
                  "warning: 1,325 of the 4,074 pairs (33%) are saturated: these are too far apart"]
        self.assertEqual(rt.compare(doc, actual, 0, **self.NORM), [])
        self.assertEqual(rt.compare(["coal pipeline: 3,143 pairs (rbh)"], actual, 0, **self.NORM),
                         ["coal pipeline: 3,143 pairs (rbh)"])

    def test_each_output_line_is_used_once_and_order_is_free(self):
        self.assertEqual(rt.compare(["a", "a"], ["a"], 0), ["a"])
        self.assertEqual(rt.compare(["b", "a"], ["a", "b"], 0), [])

    def test_exit_line_and_ls_columns(self):
        self.assertEqual(rt.compare(["exit 1"], [], 1), [])
        self.assertEqual(rt.compare(["exit 1"], [], 0), ["exit 1"])
        self.assertEqual(rt.compare(["dnds.tsv  failed.tsv  work"], ["dnds.tsv", "failed.tsv", "work"], 0), [])
        self.assertEqual(rt.compare(["dnds.tsv  failed.tsv  work"], ["dnds.tsv", "work"], 0),
                         ["dnds.tsv  failed.tsv  work"])


class TestKnownDrift(unittest.TestCase):
    def test_known_lines_are_split_off_and_entries_marked_used(self):
        known = [{"tutorial": "01", "line": "coal: summary:  1 sequence of 3 was"}, {"tutorial": "02", "line": "x"}]
        used = set()
        new, old = rt.split_known("01", ["coal: summary: 1 sequence of 3 was written ...", "x y"], known, used)
        self.assertEqual((new, len(old), used), (["x y"], 1, {0}))

    def test_the_fixture_is_valid(self):
        entries = rt.load_known_drift()
        for e in entries:
            self.assertIn(e["tutorial"], rt.RUNNABLE)
            self.assertTrue(e["line"].strip())
            self.assertTrue(e.get("since"), e)


class TestRunner(unittest.TestCase):
    @unittest.skipUnless(os.name != "nt" and shutil.which("bash"), "needs a POSIX bash (the job runs on Linux)")
    def test_one_shell_session_statuses_and_failed_commands(self):
        text = ("## 1. x\n\n```\ncd sub\nV=hello\n```\n\n```\necho $V; pwd | sed 's#.*/##'\n```\n\n```\nhello\nsub\n```\n\n"
                "```\ncat /nonexistent/x 2>/dev/null\necho \"exit $?\"\n```\n\n```\nexit 1\n```\n\n```\ngrep -q nothing /dev/null\necho after\n```\n\n"
                "```\nafter\n```\n")
        blocks = rt.plan(rt.parse_blocks(text, "01"))
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "sub"))
            res = rt.run_tutorial("01", blocks, d, os.path.join(d, "logs"), dict(os.environ), None, 60)
        st, out, errs = res[blocks[1].index]
        self.assertEqual((st, out, errs), (0, ["hello", "sub"], []))
        st, out, errs = res[blocks[2].index]
        self.assertEqual((st, out), (0, ["exit 1"]))
        self.assertTrue(errs and rt.failures_allowed(blocks[2]))
        st, out, errs = res[blocks[3].index]
        self.assertEqual(out, ["after"])
        self.assertTrue(errs)                              # a failed command, with nothing documented about it
        self.assertFalse(rt.failures_allowed(blocks[3]))  # ... is a difference


class TestTheTutorials(unittest.TestCase):
    """Every tutorial the job runs still has commands, outputs to check and (where listed) its facts in the text:
    a change of format that the runner cannot read fails here, not silently in CI."""

    def test_every_runnable_tutorial_parses_into_checks(self):
        for t in rt.RUNNABLE:
            with open(os.path.join(rt.TUTORIALS_DIR, rt.TUTORIALS[t]), encoding="utf-8") as fh:
                text = fh.read()
            blocks = rt.plan(rt.parse_blocks(text, t))
            rt.apply_skips(t, blocks, {"not_full"})
            ran = [b for b in blocks if not b.skip]
            self.assertGreater(len(ran), 2, t)
            self.assertGreater(sum(len(rt.documented_lines(b)) for b in ran), 10, t)
            self.assertTrue(rt.doc_coal_version(text), f"{t}: no 'Tested with | COAL x' line")
            flat = re.sub(r"\s+", " ", text)
            for ft, rx, _, what in rt.FACTS:
                if ft == t:
                    self.assertTrue(re.search(rx, flat), f"{rt.TUTORIALS[t]}: the text of the fact '{what}' changed")
            for b in blocks:
                self.assertNotIn(b.kind, ("",), b.where)

    def test_list_mode_runs_offline(self):
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "run_tutorials.py"), "--list", "-t", "05"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("01-codon-alignment.md", r.stdout)  # 05 depends on 01 and 02
        self.assertIn("02-pairs-dnds.md", r.stdout)


if __name__ == "__main__":
    unittest.main()
