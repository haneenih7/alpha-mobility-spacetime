# Alpha mobility through space and time

R script for sliding-window analysis of genetic, geographic, and temporal distances, with 200 resampling runs of uncertain sample dates. The supplied script retains its original attribution to Liisa Loog (2017) and implements a SpaceTime-style analysis based on the work cited below.

## Contents

- `spacetime_200run_2times.R` — the analysis script, uploaded unchanged.
- `README.md` — inputs, execution instructions, parameters, and output documentation.

Input datasets and analysis results are not included.

## Requirements

Install R and the `vegan` package:

```r
install.packages("vegan")
```

The script uses base R and `vegan::mantel`. `ggplot2` is commented out and is not required. No particular R or package version is pinned; record `sessionInfo()` with your results.

## Input files

Place these three comma-separated files in the working directory from which you run the script. All files must have a header row.

| File | Expected contents |
| --- | --- |
| `Matrixs.csv` | Numeric, square, symmetric genetic distance matrix. The spelling and capitalization are required. |
| `geo_dist.csv` | Numeric, square, symmetric geographic distance matrix, in kilometres for the configured 500-km cutoff. |
| `time.csv` | One row per individual: column 1 = sample identifier; column 2 = upper numeric date bound; column 3 = lower numeric date bound, in years ago. |

The two matrices must have the same dimensions and sample ordering, matching the rows of `time.csv`. Matching is positional: the script does not join or reorder samples by identifier. Matrix files should contain only the numeric distance columns, without an additional sample-ID or row-number column. Column names in `time.csv` are arbitrary because the script reads columns by position. Require lower bound <= upper bound and finite numeric dates.

Example date-file structure (illustrative only):

```csv
sample_id,upper_years_ago,lower_years_ago
sample_A,3100,2900
sample_B,2700,2500
```

## Run

Put the script and input files together, then run from that directory:

```sh
Rscript spacetime_200run_2times.R
```

Paths are relative to the current working directory, not automatically to the script location. There are no command-line parameters; edit the settings inside the script when needed.

For a reproducible random sequence without changing the script:

```sh
Rscript -e 'set.seed(123); source("spacetime_200run_2times.R"); sessionInfo()' > run.log 2>&1
```

Use the same inputs, settings, R version, and package versions when reproducing results. The script does not set a random seed itself.

## Default analysis settings

| Setting | Value |
| --- | --- |
| Resampling runs (`nr.run`) | 200 |
| Sliding windows (`window.moves`) | 33 |
| First window start (`start.windows.at`) | 4,200 years ago |
| Window width (`gap.window`) | 2,000 years |
| Window stride (`move.window.by`) | 100 years |
| Geographic exclusion | Pairs with distance < 500 km |
| Angle grid (`nsim`) | 200 angles from 0 to just below 90 degrees |
| Correlation method | Spearman Mantel statistic |

The first window spans 4,200 to 2,200 years ago; the final window spans 1,000 to -1,000 years ago. These are the literal bounds produced by the current settings. Adjust the window configuration if a different temporal range is intended.

For each run, the script draws one date uniformly between each sample's date bounds and calculates pairwise temporal distances. In each window, it evaluates the correlation of genetic distance with a combined distance:

```text
combined distance = sqrt(geographic distance^2 + (temporal distance * tan(angle))^2)
```

`maxratio` records the angle in degrees at the highest observed correlation. It is not a direct estimate in km/year. Each window also receives one null comparison per run, obtained by independently permuting the geographic and temporal sample orderings.

## Output

The script writes `alpha_500km_w2000.txt` to the working directory, overwriting an existing file with that name. A complete run produces 6,600 rows (33 windows x 200 runs).

The file is whitespace-delimited, despite its `.txt` extension. It is written with the default `write.table()` row names and quoting. Read it with:

```r
results <- read.table("alpha_500km_w2000.txt", header = TRUE, row.names = 1,
                      check.names = FALSE)
```

| Column | Meaning |
| --- | --- |
| `window_number` | Window index, 1–33 |
| `time_from`, `time_to` | Window bounds in years ago |
| `nsamples` | Number of individuals selected in the window |
| `maxratio` | Angle at maximum observed correlation |
| `maxcor` | Maximum observed correlation over the angle grid |
| `peak` | Maximum correlation minus the larger endpoint correlation |
| `geo_corr`, `time_corr` | Correlations at the geographic and near-temporal endpoints |
| `p_maxratio`, `p_maxcor`, `p_peak`, `p_geo_corr`, `p_time_corr` | Corresponding statistics from the permuted comparison |
| `NA` | Unused placeholder column, filled with missing values |
| `run` | Resampling run index, 1–200 |

The `p_` prefix denotes **permuted statistics, not calculated P-values**. The Mantel calls use `permutations = 0`; the script handles its own matrix permutations and does not calculate final empirical P-values. Full observed angle-correlation vectors are retained in the in-memory `corvec.array` but are not written to disk.

## Implementation notes

- Window membership uses the original lower date bound (column 3 of `time.csv`), with inclusive bounds. It does **not** use the randomly drawn date or full date-range overlap; resampled dates affect temporal distances.
- Pairs closer than 500 km are set to missing in the observed matrices. The null calculation uses the masked genetic matrix with separately permuted, original geographic and temporal matrices.
- Missing correlation statistics are replaced with zero. If several angles tie for the maximum, `maxratio` is set to zero. Interpret these cases carefully.
- No automatic input validation or handling of empty, very small, or degenerate windows is provided. Check that enough informative pairs remain after filtering.
- The script was checked for R syntax before upload. The complete analysis was not executed as part of this upload; scientific outputs have not been validated here.

## Attribution and reference

The script header credits **Liisa Loog, 2017**. That attribution is preserved.

Loog, L., Mirazón Lahr, M., Kovacevic, M., Manica, A., Eriksson, A., and Thomas, M. G. (2017). *Estimating mobility using sparse data: Application to human genetic variation*. Proceedings of the National Academy of Sciences, 114(46), 12213–12218. [https://doi.org/10.1073/pnas.1703642114](https://doi.org/10.1073/pnas.1703642114).

The source header contains a placeholder DOI; the complete reference is supplied here. No new software license is assigned by this upload.
