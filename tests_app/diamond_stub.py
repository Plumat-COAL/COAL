"""Test stand-in for DIAMOND (NOT DIAMOND): supports just `makedb --in F -d DB` and
`blastp -q Q -d DB -o OUT --outfmt 6 qseqid sseqid bitscore evalue ...`, scoring hits by shared 3-mers."""
import shutil
import sys


def read(path):
    recs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs[name] = ""
        elif name:
            recs[name] += line
    return recs


def opt(args, flag):
    return args[args.index(flag) + 1]


def main(argv):
    cmd, args = argv[0], argv[1:]
    if cmd == "makedb":
        shutil.copy(opt(args, "--in"), opt(args, "-d") + ".dmnd")
        return 0
    if cmd == "blastp":
        q = read(opt(args, "-q"))
        db = read(opt(args, "-d") + ".dmnd")
        kdb = {s: {seq[i:i + 3] for i in range(len(seq) - 2)} for s, seq in db.items()}
        with open(opt(args, "-o"), "w") as out:
            for qid, seq in q.items():
                kq = {seq[i:i + 3] for i in range(len(seq) - 2)}
                hits = sorted(((len(kq & ks), s) for s, ks in kdb.items()), reverse=True)[:5]
                for score, s in hits:
                    if score >= 5:
                        out.write(f"{qid}\t{s}\t{float(score)}\t1e-20\n")
        return 0
    print(f"stub: unsupported command {cmd}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
