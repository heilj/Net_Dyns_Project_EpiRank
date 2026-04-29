# -*- encoding: utf-8 -*-

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import additional_analysis as aa
from EpiRank import data as epirank_data
from EpiRank import epirank


g, town_data = epirank_data.load_taiwan_dataset(
    REPO_ROOT / "data",
    exclude_selfloop=False,
)
print("loaded xlsx dataset:", len(town_data), "towns")

epi_vals05 = epirank.run_epirank(
    g,
    daytime=0.5,
    d=0.95,
    exfac=epirank_data.make_population_exfac(g, town_data),
)
baselines = aa.calculate_metrices(g, d=0.95)

print("done")
print("top 5 epirank values:")
for node, value in sorted(epi_vals05.items(), key=lambda item: item[1], reverse=True)[:5]:
    print(node, value)

print("baseline metrics:", [name for name, _ in baselines])
