from pathlib import Path
from typing import Tuple
import numpy as np
from datetime import datetime


def get_time_and_memory(outfile: Path) -> Tuple:
    times = {}
    memory = None
    with outfile.open() as of:
        for line in of:
            if line.startswith("["):
                if "start:" in line or "end:" in line:
                    dt_str, rest = line.split("]")
                    event, key = rest.strip().split(":")
                    event = event.strip()
                    key = key.strip()
                    dt = datetime.strptime(dt_str[1:], "%Y-%m-%d %H:%M:%S")
                    if key not in times:
                        times[key] = {}
                    times[key][event] = dt
                else:
                    memory = [float(i) for i in line.strip()[1:-1].split(",")]
    return memory, times         


def run(folddock_output_dir: Path,
        msa_output_dir: Path, 
        interactions_file: Path, 
        output_file: Path
        ):
    assert folddock_output_dir.is_dir()
    assert interactions_file.is_file()
    interactions = sorted([l.strip().split() for l in interactions_file.open()])    
    print(f"processing {len(interactions)} interactions")
    proteins = set()
    for p1, p2 in interactions:
        proteins.add(p1)
        proteins.add(p2)
    print(f"processing {len(proteins)} proteins")

    protein_to_results = {}
    for outfile in sorted(msa_output_dir.glob("**/*.out")):
        with outfile.open() as of:
            try:
                prot = of.readline().split("=")[-1].strip().split()[-1]
            except (ValueError, IndexError):
                continue
            if prot in proteins and prot not in protein_to_results:
                protein_to_results[prot] = outfile
    print(f"found {len(protein_to_results)} protein result files")
    interaction_to_results = {}
    for outfile in sorted(folddock_output_dir.glob("**/*.out")):
        with outfile.open() as of:
            try:
                p1, p2 = of.readline().split("=")[-1].strip().split()
            except ValueError:
                continue
            if [p1, p2] in interactions and (p1, p2) not in interaction_to_results:
                interaction_to_results[p1,p2] = outfile
    print(f"found {len(interaction_to_results)} result files")
    resources = {}
    for (p1, p2), outfile in interaction_to_results.items():
        memory, times = get_time_and_memory(outfile)
        resources[p1, p2] = {
            "mem_array": memory,
            "times": times,
        }
    print(list(resources[p1,p2]["times"].keys()))
    for prot, outfile in protein_to_results.items():
        memory, times = get_time_and_memory(outfile)
        resources[prot] = {
            "mem_array": memory, 
            "times": times
        }
    stages = ['UNALIGN', 'CDHIT', 'oxmatch', 'fuse_msas', 'AlphaFold 1', 'AlphaFold 2', 'AlphaFold 3', 'AlphaFold 4', 'AlphaFold 5']
    with output_file.open("w") as of:
        of.write("protein1\tprotein2\tmem min\tmem mean\tmem max\ttime UNALIGN\ttime CDHIT\ttime oxmatch\ttime fuse_msas\ttime AlphaFold 1\ttime AlphaFold 2\ttime AlphaFold 3\ttime AlphaFold 4\ttime AlphaFold 5\ttime HHBLITS\ttime features\ttime prediction\ttime total\n")
        for (p1, p2), _ in interaction_to_results.items():
            l = f"{p1}\t{p2}"
            res = resources[p1, p2]
            times = res["times"]
            l += f"\t{min(res['mem_array'])}"
            l += f"\t{np.mean(res['mem_array'])}"
            l += f"\t{max(res['mem_array'])}"
            feature_duration = 0
            prediction_duration = 0
            for stage in stages:
                try:
                    duration = times[stage]["end"] - times[stage]["start"]
                    l += f"\t{duration.seconds}"
                    if stage.startswith("Alpha"):
                        prediction_duration += duration.seconds     
                    else:
                        feature_duration += duration.seconds    
                except KeyError:
                    print(stage, times)
                    l += "\t-1"
            time_hhblits = resources[p1]["times"]["HHBLITS"]
            duration_hhblits = (time_hhblits["end"] - time_hhblits["start"]).seconds
            time_hhblits = resources[p2]["times"]["HHBLITS"]
            duration_hhblits += (time_hhblits["end"] - time_hhblits["start"]).seconds
            feature_duration += duration_hhblits
            l += f"\t{duration_hhblits}\t{feature_duration}\t{prediction_duration}\t{feature_duration + prediction_duration}\n"
            of.write(l)
            

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(
        description="Extract the list of interactions that finished successfully")
    parser.add_argument("-d", "--folddock-output-dir", required=True,
                        help="Path to the output directory")
    parser.add_argument("-m", "--msa-output-dir", required=True,
                        help="Path to the MSA output directory")
    parser.add_argument("-i", "--interactions-file", required=True,
                        help="Path to the interactions file")
    parser.add_argument("-o", "--output-file", required=True,
                        help="Path to the output file")
    args = parser.parse_args()

    run(Path(args.folddock_output_dir),
        Path(args.msa_output_dir),
        Path(args.interactions_file),
        Path(args.output_file))
