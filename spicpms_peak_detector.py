#!/usr/bin/env python3
"""
Short-dwell ICP-MS peak/event detector

Purpose:
    Detect nanoparticle or single-cell transient events from short-dwell
    time-resolved ICP-MS data without using Excel formulas.

Main idea:
    1. Read raw counts from a CSV file.
    2. Estimate background from the first N baseline points.
    3. Set a threshold.
    4. Flag bins above threshold.
    5. Group adjacent above-threshold bins into one event.
    6. Export an event summary and optional row-by-row flagged data.

Example:
    python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --count-column Counts

If your CSV has only one numeric column, the program will use that as counts automatically.
"""

import argparse
import csv
import math
import statistics
from pathlib import Path


def parse_float(value):
    """Convert text to float; return None if it cannot be converted."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_csv_numeric(input_path, count_column=None, time_column=None):
    """
    Read CSV data and return:
        times: list of float or None
        counts: list of float
        headers: original column headers
        selected_count_column: column used as counts
        selected_time_column: column used as time, if any
    """
    input_path = Path(input_path)

    with input_path.open("r", newline="", encoding="utf-8-sig") as f:
        sample = f.read(4096)
        f.seek(0)

        # Try to detect delimiter automatically.
        # Single-column CSV files can make Sniffer fail, so use comma as a safe fallback.
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel

        try:
            has_header = csv.Sniffer().has_header(sample)
        except Exception:
            has_header = True

        reader = csv.reader(f, dialect)

        if has_header:
            headers = next(reader)
        else:
            first_row = next(reader)
            headers = [f"Column_{i+1}" for i in range(len(first_row))]
            reader = [first_row] + list(reader)

        rows = list(reader)

    if not rows:
        raise ValueError("The input file has no data rows.")

    # Map header names.
    header_to_index = {h.strip(): i for i, h in enumerate(headers)}

    # Select count column.
    if count_column:
        if count_column not in header_to_index:
            raise ValueError(
                f"Count column '{count_column}' not found. Available columns: {headers}"
            )
        count_idx = header_to_index[count_column]
        selected_count_column = count_column
    else:
        # Auto-detect numeric columns by checking how many values are numeric.
        numeric_scores = []
        for col_idx, header in enumerate(headers):
            numeric_count = 0
            checked = 0
            for row in rows[:200]:
                if col_idx < len(row):
                    checked += 1
                    if parse_float(row[col_idx]) is not None:
                        numeric_count += 1
            numeric_scores.append((numeric_count, col_idx, header))

        numeric_scores.sort(reverse=True)

        if not numeric_scores or numeric_scores[0][0] == 0:
            raise ValueError("No numeric column found. Please specify --count-column.")

        # If there are multiple numeric columns, use the last strongly numeric column.
        best_score = numeric_scores[0][0]
        candidates = [x for x in numeric_scores if x[0] == best_score]
        chosen = max(candidates, key=lambda x: x[1])
        count_idx = chosen[1]
        selected_count_column = chosen[2]

    # Select time column, optional.
    selected_time_column = None
    time_idx = None
    if time_column:
        if time_column not in header_to_index:
            raise ValueError(
                f"Time column '{time_column}' not found. Available columns: {headers}"
            )
        time_idx = header_to_index[time_column]
        selected_time_column = time_column

    counts = []
    times = []

    for row_number, row in enumerate(rows, start=2):
        if count_idx >= len(row):
            continue
        count_value = parse_float(row[count_idx])
        if count_value is None:
            continue

        counts.append(count_value)

        if time_idx is not None and time_idx < len(row):
            time_value = parse_float(row[time_idx])
        else:
            time_value = None
        times.append(time_value)

    if not counts:
        raise ValueError("No valid count values were found.")

    return times, counts, headers, selected_count_column, selected_time_column


def estimate_background(counts, baseline_points, method):
    """Estimate background mean and noise from baseline points."""
    n = min(max(1, baseline_points), len(counts))
    baseline = counts[:n]

    bg_mean = statistics.mean(baseline)

    if len(baseline) > 1:
        bg_std = statistics.stdev(baseline)
    else:
        bg_std = 0.0

    if method == "meanstd":
        noise = bg_std

    elif method == "poisson":
        # Useful for counting data. Add 1 to avoid zero-noise problems.
        noise = math.sqrt(max(bg_mean, 0.0) + 1.0)

    elif method == "robust":
        # Median absolute deviation. Good when baseline has occasional spikes.
        bg_median = statistics.median(baseline)
        deviations = [abs(x - bg_median) for x in baseline]
        mad = statistics.median(deviations)
        bg_mean = bg_median
        noise = 1.4826 * mad

    else:
        raise ValueError("method must be 'meanstd', 'poisson', or 'robust'.")

    # Avoid a zero threshold when the baseline is very flat.
    noise = max(noise, 1e-12)

    return bg_mean, noise, n


def detect_events(
    counts,
    times,
    dwell_us,
    bg_mean,
    threshold,
    min_bins,
    merge_gap_bins,
):
    """
    Detect events by grouping adjacent bins above threshold.
    Allows small gaps between above-threshold bins if merge_gap_bins > 0.
    """
    above = [x >= threshold for x in counts]

    groups = []
    i = 0
    n = len(counts)

    while i < n:
        if not above[i]:
            i += 1
            continue

        start = i
        end = i
        last_above = i
        gap_count = 0
        i += 1

        while i < n:
            if above[i]:
                end = i
                last_above = i
                gap_count = 0
            else:
                gap_count += 1
                if gap_count <= merge_gap_bins:
                    end = i
                else:
                    end = last_above
                    break
            i += 1

        groups.append((start, end))

    event_rows = []
    event_id_for_row = [""] * n

    event_id = 0
    for start, end in groups:
        num_bins = end - start + 1

        # Count how many bins in the group are actually above threshold.
        above_bins = sum(1 for j in range(start, end + 1) if above[j])

        if above_bins < min_bins:
            continue

        event_id += 1

        segment = counts[start : end + 1]
        peak_height = max(segment)
        local_peak_offset = segment.index(peak_height)
        peak_index = start + local_peak_offset

        total_counts = sum(segment)
        area_above_background = sum((x - bg_mean) for x in segment)

        start_time_us = times[start] if times[start] is not None else start * dwell_us
        end_time_us = times[end] if times[end] is not None else (end + 1) * dwell_us
        peak_time_us = times[peak_index] if times[peak_index] is not None else peak_index * dwell_us

        # If time was generated from dwell time, use bin width for duration.
        if times[start] is None or times[end] is None:
            duration_us = num_bins * dwell_us
        else:
            duration_us = end_time_us - start_time_us

        for j in range(start, end + 1):
            event_id_for_row[j] = event_id

        event_rows.append(
            {
                "Event_ID": event_id,
                "Start_Row": start + 1,
                "End_Row": end + 1,
                "Peak_Row": peak_index + 1,
                "Start_Time_us": start_time_us,
                "End_Time_us": end_time_us,
                "Peak_Time_us": peak_time_us,
                "Duration_us": duration_us,
                "Total_Bins": num_bins,
                "Above_Threshold_Bins": above_bins,
                "Peak_Height_Counts": peak_height,
                "Total_Counts": total_counts,
                "Area_Above_Background": area_above_background,
            }
        )

    return above, event_id_for_row, event_rows


def write_event_summary(output_path, event_rows):
    fieldnames = [
        "Event_ID",
        "Start_Row",
        "End_Row",
        "Peak_Row",
        "Start_Time_us",
        "End_Time_us",
        "Peak_Time_us",
        "Duration_us",
        "Total_Bins",
        "Above_Threshold_Bins",
        "Peak_Height_Counts",
        "Total_Counts",
        "Area_Above_Background",
    ]

    with Path(output_path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(event_rows)


def write_flagged_data(output_path, counts, times, dwell_us, above, event_id_for_row):
    with Path(output_path).open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Row", "Time_us", "Counts", "Above_Threshold", "Event_ID"])

        for i, count in enumerate(counts):
            time_us = times[i] if times[i] is not None else i * dwell_us
            writer.writerow([i + 1, time_us, count, int(above[i]), event_id_for_row[i]])


def write_parameters(output_path, params):
    with Path(output_path).open("w", encoding="utf-8") as f:
        for key, value in params.items():
            f.write(f"{key}: {value}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Detect short-dwell ICP-MS transient peaks/events from CSV data."
    )

    parser.add_argument("input_csv", help="Input CSV file containing time-resolved counts.")
    parser.add_argument(
        "--count-column",
        default=None,
        help="Column name containing counts. If omitted, the program auto-detects a numeric column.",
    )
    parser.add_argument(
        "--time-column",
        default=None,
        help="Optional column name containing time. If omitted, time is generated from dwell time.",
    )
    parser.add_argument(
        "--dwell-us",
        type=float,
        default=50.0,
        help="Dwell time in microseconds. Default = 50.",
    )
    parser.add_argument(
        "--baseline-points",
        type=int,
        default=1000,
        help="Number of starting points used to estimate background. Default = 1000.",
    )
    parser.add_argument(
        "--background-method",
        choices=["meanstd", "poisson", "robust"],
        default="poisson",
        help="Background/noise method. Default = poisson.",
    )
    parser.add_argument(
        "--threshold-sigma",
        type=float,
        default=5.0,
        help="Threshold multiplier. Threshold = background + threshold_sigma * noise. Default = 5.",
    )
    parser.add_argument(
        "--threshold-counts",
        type=float,
        default=None,
        help="Optional fixed threshold in counts. If used, this overrides automatic threshold.",
    )
    parser.add_argument(
        "--min-bins",
        type=int,
        default=1,
        help="Minimum number of above-threshold bins required for an event. Default = 1.",
    )
    parser.add_argument(
        "--merge-gap-bins",
        type=int,
        default=0,
        help="Merge events separated by this many below-threshold bins. Default = 0.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output prefix. Default = input filename + '_processed'.",
    )
    parser.add_argument(
        "--no-flagged-output",
        action="store_true",
        help="Do not export row-by-row flagged data. Faster and smaller output.",
    )

    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if args.output_prefix:
        prefix = Path(args.output_prefix)
    else:
        prefix = input_path.with_suffix("")
        prefix = Path(str(prefix) + "_processed")

    times, counts, headers, selected_count_column, selected_time_column = read_csv_numeric(
        input_path,
        count_column=args.count_column,
        time_column=args.time_column,
    )

    bg_mean, noise, actual_baseline_points = estimate_background(
        counts,
        baseline_points=args.baseline_points,
        method=args.background_method,
    )

    if args.threshold_counts is not None:
        threshold = args.threshold_counts
        threshold_source = "fixed user threshold"
    else:
        threshold = bg_mean + args.threshold_sigma * noise
        threshold_source = (
            f"background + {args.threshold_sigma} * noise "
            f"using {args.background_method} method"
        )

    above, event_id_for_row, event_rows = detect_events(
        counts=counts,
        times=times,
        dwell_us=args.dwell_us,
        bg_mean=bg_mean,
        threshold=threshold,
        min_bins=args.min_bins,
        merge_gap_bins=args.merge_gap_bins,
    )

    summary_path = Path(str(prefix) + "_event_summary.csv")
    params_path = Path(str(prefix) + "_parameters.txt")
    flagged_path = Path(str(prefix) + "_flagged_data.csv")

    write_event_summary(summary_path, event_rows)

    if not args.no_flagged_output:
        write_flagged_data(
            flagged_path,
            counts=counts,
            times=times,
            dwell_us=args.dwell_us,
            above=above,
            event_id_for_row=event_id_for_row,
        )

    params = {
        "Input file": input_path,
        "Rows analyzed": len(counts),
        "Selected count column": selected_count_column,
        "Selected time column": selected_time_column if selected_time_column else "Generated from dwell time",
        "Dwell time (us)": args.dwell_us,
        "Background method": args.background_method,
        "Baseline points used": actual_baseline_points,
        "Estimated background": bg_mean,
        "Estimated noise": noise,
        "Threshold source": threshold_source,
        "Threshold counts": threshold,
        "Minimum above-threshold bins": args.min_bins,
        "Merge gap bins": args.merge_gap_bins,
        "Events detected": len(event_rows),
        "Event summary output": summary_path,
        "Flagged row-by-row output": "Not exported" if args.no_flagged_output else flagged_path,
    }
    write_parameters(params_path, params)

    print("Done.")
    print(f"Rows analyzed: {len(counts)}")
    print(f"Background: {bg_mean:.6g}")
    print(f"Noise: {noise:.6g}")
    print(f"Threshold: {threshold:.6g}")
    print(f"Events detected: {len(event_rows)}")
    print(f"Event summary: {summary_path}")
    if not args.no_flagged_output:
        print(f"Flagged data: {flagged_path}")
    print(f"Parameters: {params_path}")


if __name__ == "__main__":
    main()
