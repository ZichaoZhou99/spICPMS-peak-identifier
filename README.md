# spICPMS-peak-identifier
Python based spICPMS peak identifier allow directly convert time resolved data to peak. Suitble for low dwell time. 
# Short-dwell ICP-MS peak detector

This program detects nanoparticle or single-cell ICP-MS transient events without Excel formulas.

## Why this avoids Excel freezing

Excel freezes because it recalculates formulas for every row. This program calculates everything once in Python, then exports finished CSV files that can be opened in Excel.

## Basic usage

Put your raw CSV file in the same folder as the program.

Example:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50
```

If your count column is called `Counts`:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --count-column Counts
```

If your file has a time column:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --time-column Time_us --count-column Counts
```

## Useful settings

Use a stricter threshold:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --threshold-sigma 6
```

Require at least 2 adjacent above-threshold bins:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --min-bins 2
```

Merge peaks separated by 1 below-threshold bin:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --merge-gap-bins 1
```

Skip the large row-by-row output and export only event summary:

```bash
python spicpms_peak_detector.py raw_data.csv --dwell-us 50 --no-flagged-output
```

## Output files

The program creates:

1. `*_event_summary.csv`
   - Event ID
   - Start/end row
   - Start/end time
   - Duration
   - Peak height
   - Total counts
   - Area above background

2. `*_flagged_data.csv`
   - Original row number
   - Time
   - Counts
   - Above-threshold flag
   - Event ID

3. `*_parameters.txt`
   - Background estimate
   - Threshold used
   - Number of events detected
