from pathlib import Path
from os import SEEK_CUR, SEEK_END


def separate_pdb(pdb_in: Path, 
                 pdb_a: Path, 
                 pdb_b: Path):
    with pdb_in.open() as i, pdb_a.open("w") as a, pdb_b.open("w") as b:
        curr_residue = -1
        first_chain = True
        for line in i:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                residue_number = int(line[22:26].strip())
                if curr_residue < 0:
                    curr_residue = residue_number
                if residue_number - curr_residue > 10:
                    first_chain = False
                if first_chain:
                    a.write(line)
                else:
                    b.write(line[:21] + "B" + line[22:])    
                curr_residue = residue_number


def separate_results(p1, p2, fd_output_dir, output_dir):
    outdir = fd_output_dir / f"{p1}-{p2}_results"
    for i in range(1, 6):
        outfile = outdir / f"{p1}-{p2}_{i}" / "unrelaxed_model_1.pdb"
        afile = output_dir / f"{p1}-{p2}_{i}_A.pdb"
        bfile = output_dir / f"{p1}-{p2}_{i}_B.pdb"
        separate_pdb(outfile, afile, bfile)


def run(folddock_output_dir: Path,
        interactions_file: Path, 
        output_dir: Path
        ):
    assert folddock_output_dir.is_dir()
    assert output_dir.is_dir()
    assert interactions_file.is_file()
    interactions = sorted([l.strip().split() for l in interactions_file.open()])    
    print(f"processing {len(interactions)} interactions")
    for p1, p2 in interactions:
        separate_results(p1, p2, folddock_output_dir, output_dir)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Extract the list of interactions that finished successfully")
    parser.add_argument("-d", "--folddock-output-dir", required=True,
                        help="Path to the output directory")
    parser.add_argument("-i", "--interactions-file", required=True,
                        help="Path to the interactions file")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="Path to the output file")
    args = parser.parse_args()

    run(Path(args.folddock_output_dir),
        Path(args.interactions_file),
        Path(args.output_dir))
