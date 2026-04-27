# Net_Dyns_Project_EpiRank
project for Networks and Dynamics 
Modernized code for the EpiRank paper:

Huang CY, Chin WC, Wen TH, Fu YH, Tsai YS. 2019. **EpiRank: Modeling Bidirectional Disease Spread in Asymmetric Commuting Networks.** *Scientific Reports*. https://doi.org/10.1038/s41598-019-41719-8

The original project code is based on https://github.com/wcchin/EpiRank. This repository reorganizes the data and scripts for a modern Python environment and NetworkX 3.x.

## Repository Layout

- `data_flow/`: example origin-destination commuting flow data.
- `data_sars/`: SARS case-study data, including commuting flows, SARS counts, point coordinates, and Taiwan GeoJSON files.
- `scripts/EpiRank/`: modernized EpiRank implementation and analysis helpers.
- `scripts/example_modern.py`: simple command-line example using `data_flow/od_flow_data.csv`.
- `scripts/example_notebook_modern.ipynb`: notebook version of the simple example.
- `scripts/sars_run/run_for_sars_modern.py`: SARS replication runner.
- `scripts/sars_run/run_for_sars_modern.ipynb`: notebook version of the SARS replication workflow.
- `scripts/sars_run/plot_sars_results_modern.py`: plotting script for SARS CSV outputs.
- `scripts/sars_run/plot_sars_results_modern.ipynb`: notebook version of the plotting workflow.
- `results/`: default output location for generated CSV result tables.
- `requirements.txt`: Python package requirements for the modernized code.

## Environment

Use Python 3.10. The code has been tested with a modern stack compatible with NetworkX 3.x.

Create a conda environment:

```bash
conda create -n epirank-modern python=3.10
conda activate epirank-modern
pip install -r requirements.txt
```

The core requirements are:

```text
networkx>=3.0
numpy>=1.24
pandas>=2.0
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
python scripts/example_modern.py
```

This loads `data_flow/od_flow_data.csv`, builds a weighted directed commuting graph, runs EpiRank with `daytime=0.5` and `d=0.95`, and prints the top-ranked nodes plus baseline metric names.

You can also open:

```text
scripts/example_notebook_modern.ipynb
```

## Run The SARS Replication

Quick check with a coarse sensitivity grid:

```bash
python scripts/sars_run/run_for_sars_modern.py --daytime-step 0.5 --damping-step 0.5 --top-n 3
```

Full sensitivity grid used by the notebook:

```bash
python scripts/sars_run/run_for_sars_modern.py --daytime-step 0.05 --damping-step 0.05 --top-n 10
```

By default, the SARS runner writes CSV tables to:

```text
results/sars/
```

Generated files:

- `movement_comparison_modern.csv`: backward-only, bidirectional, and forward-only EpiRank correlations.
- `baseline_comparison_modern.csv`: PageRank, HITS hub, and HITS authority baseline correlations.
- `sensitivity_modern.csv`: correlation results across daytime and damping parameters.
- `top_scores_modern.csv`: top EpiRank locations for the movement settings.

To write results somewhere else:

```bash
python scripts/sars_run/run_for_sars_modern.py --output-dir /tmp/epirank-sars-results
```

The notebook equivalent is:

```text
scripts/sars_run/run_for_sars_modern.ipynb
```

## Generate Plots

After the SARS CSV tables exist in `results/sars/`, generate figures with:

```bash
python scripts/sars_run/plot_sars_results_modern.py
```

This writes PNG files to:

```text
results/sars/figures/
```

Generated figures:

- `correlation_spearman_modern.png`
- `correlation_pearson_modern.png`
- `sensitivity_spearman_modern.png`
- `sensitivity_pearson_modern.png`
- `top_scores_modern.png`

To read result tables from another location or write figures elsewhere:

```bash
python scripts/sars_run/plot_sars_results_modern.py \
  --results-dir /tmp/epirank-sars-results \
  --figures-dir /tmp/epirank-sars-figures
```

The notebook equivalent is:

```text
scripts/sars_run/plot_sars_results_modern.ipynb
```

## Notes On Imports

The scripts and notebooks locate the repository root from their file location or current notebook directory, then add `scripts/` to Python's import path. This allows imports such as:

```python
from EpiRank import epirank_modern
from EpiRank import additional_analysis_modern
```

Because paths are resolved from the repository layout, the scripts can be launched from the repository root or by absolute path from another working directory.
