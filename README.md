# Net_Dyns_Project_EpiRank
project for Networks and Dynamics 
Reorganized code for the EpiRank paper:

Huang CY, Chin WC, Wen TH, Fu YH, Tsai YS. 2019. **EpiRank: Modeling Bidirectional Disease Spread in Asymmetric Commuting Networks.** *Scientific Reports*. https://doi.org/10.1038/s41598-019-41719-8

The original project code is based on https://github.com/wcchin/EpiRank. This repository keeps the package-style script layout while adding the workbook-based data pipeline from the GUI implementation.

## Repository Layout

- `data/`: default Taiwan township workbook dataset (`bs.xlsx`, `cn.xlsx`, `Flu.xlsx`, `ev.xlsx`, `SARS.xlsx`). Same-name CSV files are supported as a fallback.
- `data_flow/`: legacy example origin-destination commuting flow data.
- `data_sars/`: SARS case-study data, including commuting flows, SARS counts, point coordinates, and Taiwan GeoJSON files.
- `scripts/EpiRank/`: reorganized EpiRank implementation and analysis helpers.
- `scripts/example.py`: simple command-line example using the default xlsx dataset.
- `scripts/example_notebook.ipynb`: notebook version of the simple example.
- `scripts/run_taiwan_dataset.py`: full xlsx-first Taiwan run that writes rankings, correlations, and classification breaks.
- `scripts/sars_run/run_for_sars.py`: SARS replication runner.
- `scripts/sars_run/run_for_sars.ipynb`: notebook version of the SARS replication workflow.
- `scripts/sars_run/plot_sars_results.py`: plotting script for SARS CSV outputs.
- `scripts/sars_run/plot_sars_results.ipynb`: notebook version of the plotting workflow.
- `results/`: default output location for generated CSV result tables.
- `requirements.txt`: Python package requirements.

## Environment

Use Python 3.10 or newer with NetworkX 3.x.

Create a conda environment:

```bash
conda create -n epirank python=3.10
conda activate epirank
pip install -r requirements.txt
```

The core requirements are:

```text
networkx>=3.0
numpy>=1.24
pandas>=2.0
openpyxl>=3.1
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
```

For notebooks, also install Jupyter if your environment does not already include it:

```bash
pip install notebook nbconvert ipykernel
```

## Run The Example

From the repository root:

```bash
python scripts/example.py
```

This loads the workbook files under `data/`, builds a weighted directed commuting graph, runs EpiRank with `daytime=0.5` and `d=0.95`, and prints the top-ranked nodes plus baseline metric names.

You can also open:

```text
scripts/example_notebook.ipynb
```

## Run The Taiwan Workbook Dataset

```bash
python scripts/run_taiwan_dataset.py
```

By default, this reads xlsx files from `data/` and writes:

- `results/taiwan/ranking.csv`
- `results/taiwan/correlations.csv`
- `results/taiwan/epirank_breaks.csv`
- `results/taiwan/figures/correlations_table.png`
- `results/taiwan/figures/epirank_scores.png`
- `results/taiwan/figures/epirank_levels.png`
- `results/taiwan/figures/level_counts.png`

If an xlsx workbook is unavailable, place a same-name CSV file in `data/` and the loader will use it as a backup.

To run the Flu/EV/SARS daytime-by-damping parameter grid used for the sensitivity analysis:

```bash
python scripts/run_taiwan_dataset.py --parameter-grid
```

This adds `results/taiwan/parameter_grid.csv` and `results/taiwan/figures/parameter_grid.png`. The grid heatmaps mark the best cell in each panel.

To check whether the high-risk locations stay stable under parameter changes:

```bash
python scripts/run_taiwan_dataset.py --stability
```

This adds `results/taiwan/parameter_stability.csv` and `results/taiwan/figures/parameter_stability.png`. The stability outputs compare each setting against the reference configuration given by `--daytime` and `--damping`, including top-10 overlap, top-20 overlap, Core-I overlap, rank Spearman, rank Kendall, and mean rank shift for the reference top towns.

## Generate Paper Figures

The GUI comments map plots to paper figures. The same figure-oriented workflow is available without GUI code:

```bash
python scripts/plot_paper_figures.py
```

This writes Figure 2, Figure 3, Figure 4, Figure 6, Figure 7, Figure 8, Figure 9, and Figure 10 to `results/paper_figures/`. Add `--grid` to also generate Figure 11 and its parameter-grid CSV.

## Run The SARS Replication

Quick check with a coarse sensitivity grid:

```bash
python scripts/sars_run/run_for_sars.py --daytime-step 0.5 --damping-step 0.5 --top-n 3
```

Full sensitivity grid used by the notebook:

```bash
python scripts/sars_run/run_for_sars.py --daytime-step 0.05 --damping-step 0.05 --top-n 10
```

By default, the SARS runner writes CSV tables to:

```text
results/sars/
```

Generated files:

- `movement_comparison.csv`: backward-only, bidirectional, and forward-only EpiRank correlations.
- `baseline_comparison.csv`: PageRank, HITS hub, and HITS authority baseline correlations.
- `sensitivity.csv`: correlation results across daytime and damping parameters.
- `top_scores.csv`: top EpiRank locations for the movement settings.

To write results somewhere else:

```bash
python scripts/sars_run/run_for_sars.py --output-dir /tmp/epirank-sars-results
```

The notebook equivalent is:

```text
scripts/sars_run/run_for_sars.ipynb
```

## Generate Plots

After the SARS CSV tables exist in `results/sars/`, generate figures with:

```bash
python scripts/sars_run/plot_sars_results.py
```

This writes PNG files to:

```text
results/sars/figures/
```

Generated figures:

- `correlation_spearman.png`
- `correlation_pearson.png`
- `sensitivity_spearman.png`
- `sensitivity_pearson.png`
- `top_scores.png`

To read result tables from another location or write figures elsewhere:

```bash
python scripts/sars_run/plot_sars_results.py \
  --results-dir /tmp/epirank-sars-results \
  --figures-dir /tmp/epirank-sars-figures
```

The notebook equivalent is:

```text
scripts/sars_run/plot_sars_results.ipynb
```

## Notes On Imports

The scripts and notebooks locate the repository root from their file location or current notebook directory, then add `scripts/` to Python's import path. This allows imports such as:

```python
from EpiRank import epirank
from EpiRank import additional_analysis
```

Because paths are resolved from the repository layout, the scripts can be launched from the repository root or by absolute path from another working directory.
