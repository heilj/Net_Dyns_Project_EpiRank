# -*- coding: utf-8 -*-

"""Data loading helpers for the Taiwan EpiRank workbook dataset.

The default source is the xlsx dataset in the repository ``data/`` directory.
CSV files with matching stems can be used as a backup when a workbook is not
available.
"""

from pathlib import Path

import networkx as nx
import pandas as pd


KEY_POST_CODE = "post_code"
KEY_DB_ID = "db_ID"
KEY_COUNTY = "county"
KEY_TOWN = "town"
KEY_POS_XY = "pos_xy"
KEY_POPULATION = "population"
KEY_AREA = "area"
KEY_DENSITY = "density"
KEY_NORMALIZED_DENSITY = "normalized_density"
KEY_AGE_0_14 = "age_0_14"
KEY_AGE_15_64 = "age_15_64"
KEY_AGE_65 = "age_65"
KEY_LOCAL_COMMUTER_TYPE1 = "local_commuter_type1"
KEY_OUT_COMMUTER_TYPE1 = "out_commuter_type1"
KEY_IN_COMMUTER_TYPE1 = "in_commuter_type1"
KEY_COMMUTER_TYPE1 = "commuter_type1"
KEY_RAILROAD_ZONE = "railroad_zone"
KEY_FLU_TOTAL_CASES = "flu_total_cases"
KEY_EV_AVERAGE_CASES = "EV_average_cases"
KEY_SARS_TOTAL_CASES = "sars_total_cases"

GTAIPEI_DB_IDS = set(range(0, 29)) | set(range(303, 310)) | set(range(330, 397))


def default_data_dir():
    return Path(__file__).resolve().parents[2] / "data"


def _resolve_source(data_dir, stem, preferred="xlsx"):
    data_dir = Path(data_dir)
    candidates = (
        [data_dir / f"{stem}.xlsx", data_dir / f"{stem}.csv"]
        if preferred == "xlsx"
        else [data_dir / f"{stem}.csv", data_dir / f"{stem}.xlsx"]
    )
    for path in candidates:
        if path.exists():
            return path
    names = ", ".join(path.name for path in candidates)
    raise FileNotFoundError(f"expected one of {names} under {data_dir}")


def _read_table(data_dir, stem, sheet_name, header=None, preferred="xlsx"):
    path = _resolve_source(data_dir, stem, preferred=preferred)
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path, sheet_name=sheet_name, header=header)
    return pd.read_csv(path, header=header)


def sorted_map(mapping):
    return sorted(mapping.items(), key=lambda kv: (-kv[1], kv[0]))


def build_basic_table_of_towns(
    town_data,
    data_dir=None,
    sheet="town_data",
    number_of_sub_towns=409,
    row_base=2,
):
    """Load township metadata from ``bs.xlsx`` or ``bs.csv``."""
    data_dir = default_data_dir() if data_dir is None else data_dir
    df = _read_table(data_dir, "bs", sheet, header=None)
    old_town = None

    for row_idx in range(number_of_sub_towns):
        row = df.iloc[row_base - 1 + row_idx]
        if pd.isna(row.iloc[0]):
            continue

        db_id = int(row.iloc[0])
        county = row.iloc[1]
        town = row.iloc[2]
        raw_population = 0 if pd.isna(row.iloc[8]) else row.iloc[8]
        sub_percentage = 0.0 if pd.isna(row.iloc[9]) else float(row.iloc[9])
        population = (
            float(raw_population / sub_percentage)
            if sub_percentage > 0
            else float(raw_population)
        )

        if town != old_town:
            town_data[db_id] = {
                KEY_COUNTY: county,
                KEY_TOWN: town,
                KEY_POS_XY: (
                    round(float(0 if pd.isna(row.iloc[6]) else row.iloc[6]), 2),
                    round(float(0 if pd.isna(row.iloc[7]) else row.iloc[7]), 2),
                ),
                KEY_POPULATION: population,
                KEY_AREA: float(0 if pd.isna(row.iloc[11]) else row.iloc[11]),
                KEY_DENSITY: float(0 if pd.isna(row.iloc[12]) else row.iloc[12]),
                KEY_NORMALIZED_DENSITY: float(0 if pd.isna(row.iloc[13]) else row.iloc[13]),
                KEY_AGE_0_14: float(0 if pd.isna(row.iloc[14]) else row.iloc[14]),
                KEY_AGE_15_64: float(0 if pd.isna(row.iloc[15]) else row.iloc[15]),
                KEY_AGE_65: float(0 if pd.isna(row.iloc[16]) else row.iloc[16]),
                KEY_LOCAL_COMMUTER_TYPE1: 0,
                KEY_OUT_COMMUTER_TYPE1: 0,
                KEY_IN_COMMUTER_TYPE1: 0,
                KEY_RAILROAD_ZONE: 0,
            }
        old_town = town
    return town_data


def _town_lookup(town_data):
    return {
        (town[KEY_COUNTY], town[KEY_TOWN]): db_id
        for db_id, town in town_data.items()
    }


def _build_reported_cases(
    town_data,
    data_dir,
    stem,
    sheet,
    target_key,
    value_type,
    number_of_towns=353,
    row_base=2,
):
    data_dir = default_data_dir() if data_dir is None else data_dir
    df = _read_table(data_dir, stem, sheet, header=None)
    lookup = _town_lookup(town_data)
    for row_idx in range(number_of_towns):
        row = df.iloc[row_base - 1 + row_idx]
        db_id = lookup.get((row.iloc[0], row.iloc[1]))
        if db_id is None:
            continue
        raw = row.iloc[2]
        town_data[db_id][target_key] = value_type(0 if pd.isna(raw) else raw)
    return town_data


def build_flu_reported_cases(town_data, data_dir=None, sheet="2009", **kwargs):
    return _build_reported_cases(
        town_data, data_dir, "Flu", sheet, KEY_FLU_TOTAL_CASES, int, **kwargs
    )


def build_ev_reported_cases(town_data, data_dir=None, sheet="2000_2008", **kwargs):
    return _build_reported_cases(
        town_data, data_dir, "ev", sheet, KEY_EV_AVERAGE_CASES, float, **kwargs
    )


def build_sars_reported_cases(town_data, data_dir=None, sheet="2003", **kwargs):
    return _build_reported_cases(
        town_data, data_dir, "SARS", sheet, KEY_SARS_TOTAL_CASES, int, **kwargs
    )


def build_commuting_network(
    g,
    town_data,
    data_dir=None,
    sheet="353C",
    number_of_towns=353,
    row_base=6,
    col_base=6,
    exclude_selfloop=False,
):
    """Build the 353-node directed commuting graph from ``cn.xlsx``/``cn.csv``."""
    data_dir = default_data_dir() if data_dir is None else data_dir
    df = _read_table(data_dir, "cn", sheet, header=None)
    cache = {}

    for row_idx in range(number_of_towns):
        r = row_base - 1 + row_idx
        row_seq_raw = df.iat[r, 0]
        if pd.isna(row_seq_raw):
            continue
        row_seq_no = int(row_seq_raw)
        if row_seq_no in cache:
            row_code, row_db_id = cache[row_seq_no]
        else:
            row_code = "" if pd.isna(df.iat[r, 1]) else str(df.iat[r, 1])
            if pd.isna(df.iat[r, 2]):
                continue
            row_db_id = int(df.iat[r, 2])
            cache[row_seq_no] = (row_code, row_db_id)

        for col_idx in range(number_of_towns):
            c = col_base - 1 + col_idx
            col_seq_raw = df.iat[0, c]
            if pd.isna(col_seq_raw):
                continue
            col_seq_no = int(col_seq_raw)
            if col_seq_no in cache:
                col_code, col_db_id = cache[col_seq_no]
            else:
                col_code = "" if pd.isna(df.iat[1, c]) else str(df.iat[1, c])
                if pd.isna(df.iat[2, c]):
                    continue
                col_db_id = int(df.iat[2, c])
                cache[col_seq_no] = (col_code, col_db_id)

            raw_commuters = df.iat[r, c]
            commuters = int(0 if pd.isna(raw_commuters) else raw_commuters)
            if commuters <= 0 or row_db_id not in town_data or col_db_id not in town_data:
                continue
            if exclude_selfloop and row_seq_no == col_seq_no:
                continue

            if not g.has_node(row_seq_no):
                g.add_node(
                    row_seq_no,
                    post_code=row_code,
                    db_ID=row_db_id,
                    posx=town_data[row_db_id][KEY_POS_XY][0],
                    posy=town_data[row_db_id][KEY_POS_XY][1],
                    population=town_data[row_db_id][KEY_POPULATION],
                )
            if not g.has_node(col_seq_no):
                g.add_node(
                    col_seq_no,
                    post_code=col_code,
                    db_ID=col_db_id,
                    posx=town_data[col_db_id][KEY_POS_XY][0],
                    posy=town_data[col_db_id][KEY_POS_XY][1],
                    population=town_data[col_db_id][KEY_POPULATION],
                )

            g.add_edge(
                row_seq_no,
                col_seq_no,
                weight=float(commuters),
                commuter_type1=float(commuters),
            )
            town_data[row_db_id][KEY_OUT_COMMUTER_TYPE1] += commuters
            town_data[col_db_id][KEY_IN_COMMUTER_TYPE1] += commuters
    return g


def load_taiwan_dataset(data_dir=None, include_sars=True, exclude_selfloop=False):
    """Load xlsx-first Taiwan township data and return ``(graph, town_data)``."""
    data_dir = default_data_dir() if data_dir is None else data_dir
    town_data = {}
    graph = nx.DiGraph()
    build_basic_table_of_towns(town_data, data_dir=data_dir)
    build_flu_reported_cases(town_data, data_dir=data_dir)
    build_ev_reported_cases(town_data, data_dir=data_dir)
    if include_sars:
        build_sars_reported_cases(town_data, data_dir=data_dir)
    build_commuting_network(
        graph,
        town_data,
        data_dir=data_dir,
        exclude_selfloop=exclude_selfloop,
    )
    return graph, town_data


def rank_from_town_data(graph, town_data, key, default=0):
    return {
        node: town_data[attrs[KEY_DB_ID]].get(key, default)
        for node, attrs in graph.nodes(data=True)
    }


def make_population_exfac(graph, town_data):
    """Return a node-keyed population vector for use as ``exfac``."""
    return {
        node: town_data[attrs[KEY_DB_ID]][KEY_POPULATION]
        for node, attrs in graph.nodes(data=True)
    }
