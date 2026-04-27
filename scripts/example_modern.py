# -*- encoding: utf-8 -*-

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import additional_analysis_modern as aa
from EpiRank import epirank_modern as epirank

import pandas as pd


df = pd.read_csv(REPO_ROOT / "data_flow" / "od_flow_data.csv", index_col=0)
print(df.head())

g = epirank.make_DiGraph(
    df,
    origin_col="origin",
    destination_col="destination",
    flow_col="flow",
    largest_connected_component=False,
    exclude_selfloop=False,
)

epi_vals05 = epirank.run_epirank(g, daytime=0.5, d=0.95)
baselines = aa.calculate_metrices(g, d=0.95)

print("done")
print("top 5 epirank values:")
for node, value in sorted(epi_vals05.items(), key=lambda item: item[1], reverse=True)[:5]:
    print(node, value)

print("baseline metrics:", [name for name, _ in baselines])
