from pathlib import Path
from os import SEEK_CUR, SEEK_END


# fast read last line taken from https://stackoverflow.com/a/18603065
def _readlast__bytes(f, sep, size, step):
    # Point cursor 'size' + 'step' bytes away from the end of the file.
    o = f.seek(0 - size - step, SEEK_END)
    # Step 'step' bytes each iteration, halt when 'sep' occurs.
    while f.read(size) != sep:
        f.seek(0 - size - step, SEEK_CUR)

def _readlast__text(f, sep, size, step):
    # Text mode, same principle but without the use of relative offsets.
    o = f.seek(0, SEEK_END)
    o = f.seek(o - size - step)
    while f.read(size) != sep:
        o = f.seek(o - step)

def readlast(f, sep, fixed = False):
    """readlast(f: io.BaseIO, sep: bytes|str, fixed: bool = False) -> bytes|str

    Return the last segment of file `f`, containing data segments separated by
    `sep`.

    Set `fixed` to True when parsing UTF-32 or UTF-16 encoded data (don't forget
    to pass the correct delimiter) in files opened in byte mode.
    """
    size = len(sep)
    step = len(sep) if (fixed is True) else (fixed or 1)
    step = size if fixed else 1
    if not size:
        raise ValueError("Zero-length separator.")
    try:
        if 'b' in f.mode:
            # Process file opened in byte mode.
            _readlast__bytes(f, sep, size, step)
        else:
            # Process file opened in text mode.
            _readlast__text(f, sep, size, step)
    except (OSError, ValueError): 
        # Beginning of file reached.
        f.seek(0, SEEK_SET)
    return f.read()


def check_pdb(outfile: Path):
    return readlast(outfile.open("rb"), b"\n").strip() == b"END"


def check_finished(p1, p2, output_dir):
    outdir = output_dir / f"{p1}-{p2}_results"
    finished = []
    for i in range(1, 6):
        outfile = outdir / f"{p1}-{p2}_{i}" / "unrelaxed_model_1.pdb"
        finished.append(outfile.exists() and check_pdb(outfile))
    return all(finished)


def run(output_dir: Path,
        interactions_file: Path, 
        output_file: Path):
    assert output_dir.is_dir()
    assert interactions_file.is_file()
    interactions = sorted([l.strip().split() for l in interactions_file.open()])    
    print(f"processing {len(interactions)} interactions")
    with output_file.open("w") as of:
        for p1, p2 in interactions:
            if check_finished(p1, p2, output_dir):
                of.write(f"{p1}\t{p2}\n") 


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Extract the list of interactions that finished successfully")
    parser.add_argument("-d", "--output-dir", required=True,
                        help="Path to the output directory")
    parser.add_argument("-i", "--interactions-file", required=True,
                        help="Path to the interactions file")
    parser.add_argument("-o", "--output-file", required=True,
                        help="Path to the output file")
    args = parser.parse_args()

    run(Path(args.output_dir),
        Path(args.interactions_file),
        Path(args.output_file))
