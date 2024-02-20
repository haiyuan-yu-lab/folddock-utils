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


def get_resource_usages(resources, p1_res, p2_res):
    memory_max = []
    memory_min = []
    memory_mean = []
    stages = ['UNALIGN', 'CDHIT', 'oxmatch', 'fuse_msas', 'AlphaFold 1',
              'AlphaFold 2', 'AlphaFold 3', 'AlphaFold 4', 'AlphaFold 5']
    times = {stage: [] for stage in stages}
    for res in resources:
        times = res["times"]
        memory_min.append(min(res["mem_array"]))
        memory_mean.append(np.mean(res["mem_array"]))
        memory_max.append(max(res["mem_array"]))
        feature_duration = 0
        for stage in stages:
            try:
                duration = times[stage]["end"] - times[stage]["start"]
                times[stage].append(duration)
                if stage.startswith("Alpha"):
                    prediction_duration += duration.seconds
                else:
                    feature_duration += duration.seconds
            except KeyError:
                print(stage, times)
        time_hhblits = p1_res["times"]["HHBLITS"]
        duration_hhblits = (time_hhblits["end"]
                            - time_hhblits["start"]).seconds
        time_hhblits = p2_res["times"]["HHBLITS"]
        duration_hhblits += (time_hhblits["end"]
                             - time_hhblits["start"]).seconds
        feature_duration += duration_hhblits
        feature_duration += (max(times["UNALIGN"])
                             + max(times["CDHIT"])
                             + max(times["oxmatch"]))
        prediction_duration = (max(times["AlphaFold 1"]),
                               + max(times["AlphaFold 2"])
                               + max(times["AlphaFold 3"])
                               + max(times["AlphaFold 4"])
                               + max(times["AlphaFold 5"]))
    return (
        max(memory_min),
        max(memory_mean),
        max(memory_max),
        max(times["UNALIGN"]),
        max(times["CDHIT"]),
        max(times["oxmatch"]),
        max(times["fuse_msas"]),
        max(times["AlphaFold 1"]),
        max(times["AlphaFold 2"]),
        max(times["AlphaFold 3"]),
        max(times["AlphaFold 4"]),
        max(times["AlphaFold 5"]),
        feature_duration,
        prediction_duration,
        feature_duration + prediction_duration)


def run(folddock_output_dir: Path,
        msa_output_dir: Path,
        interactions_file: Path,
        output_file: Path
        ):
    assert folddock_output_dir.is_dir()
    assert interactions_file.is_file()
    interactions = sorted([line.strip().split()
                           for line in interactions_file.open()])
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
    files_found = 0
    for outfile in sorted(folddock_output_dir.glob("**/*.out")):
        with outfile.open() as of:
            try:
                p1, p2 = of.readline().split("=")[-1].strip().split()
            except ValueError:
                continue
            if (p1, p2) not in interaction_to_results:
                interaction_to_results[p1, p2] = []
            if [p1, p2] in interactions:
                interaction_to_results[p1, p2].append(outfile)
                files_found += 1
    print(f"found {len(files_found)} result files")

    resources = {}
    for (p1, p2), outfiles in interaction_to_results.items():
        for outfile in outfiles:
            memory, times = get_time_and_memory(outfile)
            if (p1, p2) not in resources:
                resources[p1, p2] = {}
            resources[p1, p2].append({
                "mem_array": memory,
                "times": times,
            })
    print(list(resources[p1, p2]["times"].keys()))
    for prot, outfile in protein_to_results.items():
        memory, times = get_time_and_memory(outfile)
        resources[prot] = {
            "mem_array": memory,
            "times": times
        }
    with output_file.open("w") as of:
        of.write("protein1\tprotein2\tmem min\tmem mean\tmem max\t"
                 "time UNALIGN\ttime CDHIT\ttime oxmatch\ttime fuse_msas\t"
                 "time AlphaFold 1\ttime AlphaFold 2\ttime AlphaFold 3\t"
                 "time AlphaFold 4\ttime AlphaFold 5\ttime HHBLITS\t"
                 "time features\ttime prediction\ttime total\n")
        for (p1, p2), _ in interaction_to_results.items():
            (memory_min,
             memory_mean,
             memory_max,
             timesUNALIGN,
             timesCDHIT,
             timesoxmatch,
             timesfuse_msas,
             timesAlphaFold1,
             timesAlphaFold2,
             timesAlphaFold3,
             timesAlphaFold4,
             timesAlphaFold5,
             feature_duration,
             prediction_duration,
             total_duration) = get_resource_usages(resources[p1, p2],
                                                   resources[p1],
                                                   resources[p2])
            of.write(f"{p1}\t{p2}"
                     f"\t{memory_min}"
                     f"\t{memory_mean}"
                     f"\t{memory_max}"
                     f"\t{timesUNALIGN}"
                     f"\t{timesCDHIT}"
                     f"\t{timesoxmatch}"
                     f"\t{timesfuse_msas}"
                     f"\t{timesAlphaFold1}"
                     f"\t{timesAlphaFold2}"
                     f"\t{timesAlphaFold3}"
                     f"\t{timesAlphaFold4}"
                     f"\t{timesAlphaFold5}"
                     f"\t{feature_duration}"
                     f"\t{prediction_duration}"
                     f"\t{total_duration}\n")


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
