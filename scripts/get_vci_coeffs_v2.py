#!/usr/bin/env python3

import glob
import re
import argparse


def parse_vci_file(filename):
    config_pattern = re.compile(
        r"\[\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\s*\]\s*\(([-+0-9.eE]+)\)"
    )

    freq_pattern = re.compile(
        r"\|\s*([-+0-9.eE]+)\s*\|\s*([-+0-9.eE]+)"
    )

    configs = []
    frequency = None

    with open(filename, "r") as f:
        for line in f:
            if frequency is None:
                match = freq_pattern.search(line)
                if match:
                    frequency = float(match.group(2))
                    continue

            match = config_pattern.search(line)
            if match:
                quanta = tuple(int(match.group(i)) for i in range(1, 7))
                ci2 = float(match.group(7))
                configs.append((ci2, quanta))

    configs.sort(reverse=True)
    return frequency, configs


def find_files_by_frequency(files, target_freq, tolerance=None):
    matches = []

    for filename in sorted(files):
        frequency, configs = parse_vci_file(filename)

        if frequency is None:
            continue

        diff = abs(frequency - target_freq)

        if tolerance is None or diff <= tolerance:
            matches.append((diff, filename, frequency, configs))

    matches.sort(key=lambda x: x[0])
    return matches


def make_config_string(quanta):
    return "(" + ",".join(str(x) for x in quanta) + ")"


def latex_rows(label, frequency, configs, nterms=3, threshold=0.01):
    kept = [(ci2, q) for ci2, q in configs if ci2 >= threshold][:nterms]

    if not kept:
        return "% No configurations found above threshold."

    rows = []

    for i, (ci2, quanta) in enumerate(kept):
        config = make_config_string(quanta)

        if i == 0:
            rows.append(
                f"{label} & {frequency:.2f} & ${config}$ & {ci2:.4f} \\\\[2pt]"
            )
        else:
            rows.append(
                f"        &       & ${config}$ & {ci2:.4f} \\\\[2pt]"
            )

    rows.append(r"\hline")
    rows.append(r"\addlinespace")

    return "\n".join(rows)


def print_result(filename, state_label, frequency, configs, args):
    kept = [
        (ci2, quanta)
        for ci2, quanta in configs
        if ci2 >= args.threshold
    ][:args.num]

    print()
    print("=" * 70)
    print(f"{filename}")
    print(f"State = {state_label}")
    print(f"Frequency = {frequency:.6f} cm^-1")
    print("-" * 70)
    print("Leading VCI Contributions")
    print()

    if not kept:
        print("No configurations found above threshold.")
    else:
        for i, (ci2, quanta) in enumerate(kept, start=1):
            config = make_config_string(quanta)
            print(f"{i:2d}.  {config:20s}   {ci2:.8f}")

    print()
    print("LaTeX:")
    print("-" * 70)

    latex = latex_rows(
        state_label,
        frequency,
        configs,
        args.num,
        args.threshold
    )

    print(latex)
    print()

    return latex


def main():
    parser = argparse.ArgumentParser(
        description="Find leading VCI |Ci|^2 configurations and print human-readable and LaTeX table rows."
    )

    parser.add_argument(
        "files",
        nargs="*",
        default=glob.glob("VCI_MODES*.dat"),
        help="VCI files to analyze. Default: VCI_MODES*.dat"
    )

    parser.add_argument(
        "-n",
        "--num",
        type=int,
        default=3,
        help="Number of leading configurations to print. Default: 3"
    )

    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.01,
        help="Minimum |Ci|^2 value to include. Default: 0.01"
    )

    parser.add_argument(
        "--freq",
        type=float,
        default=None,
        help="Target frequency in cm^-1. The script will find the closest matching VCI file."
    )

    parser.add_argument(
        "--tol",
        type=float,
        default=None,
        help="Frequency tolerance in cm^-1. If omitted, the closest match is used."
    )

    parser.add_argument(
        "--label",
        type=str,
        default=None,
        help=r"LaTeX state label, e.g. '$\nu_b$'."
    )

    parser.add_argument(
        "--num-matches",
        type=int,
        default=1,
        help="Number of frequency matches to print. Default: 1"
    )

    parser.add_argument(
        "--no-label-prompt",
        action="store_true",
        help="Use filenames as labels instead of prompting for state labels."
    )

    args = parser.parse_args()

    if not args.files:
        print("No VCI files found.")
        return

    all_latex_rows = []

    if args.freq is not None:
        matches = find_files_by_frequency(
            args.files,
            args.freq,
            tolerance=args.tol
        )

        if not matches:
            print("No matching VCI files found.")
            return

        matches = matches[:args.num_matches]

        print()
        print(f"Target frequency = {args.freq:.6f} cm^-1")

        if args.tol is not None:
            print(f"Tolerance = {args.tol:.6f} cm^-1")

        print("Matches:")
        for diff, filename, frequency, configs in matches:
            print(f"  {filename:25s}  {frequency:12.6f} cm^-1  diff = {diff:.6f}")

        for diff, filename, frequency, configs in matches:
            if args.label is not None and len(matches) == 1:
                state_label = args.label
            elif args.no_label_prompt:
                state_label = filename
            else:
                state_label = input(
                    f"Label for {filename} ({frequency:.2f} cm^-1): "
                )

            latex = print_result(
                filename,
                state_label,
                frequency,
                configs,
                args
            )

            all_latex_rows.append(latex)

    else:
        for filename in sorted(args.files):
            frequency, configs = parse_vci_file(filename)

            if frequency is None:
                print(f"Skipping {filename}: no frequency found.")
                continue

            if args.no_label_prompt:
                state_label = filename
            else:
                state_label = input(
                    f"Label for {filename} ({frequency:.2f} cm^-1): "
                )

            latex = print_result(
                filename,
                state_label,
                frequency,
                configs,
                args
            )

            all_latex_rows.append(latex)

    print()
    print("=" * 70)
    print("Combined LaTeX rows")
    print("=" * 70)
    print("\n".join(all_latex_rows))


if __name__ == "__main__":
    main()
