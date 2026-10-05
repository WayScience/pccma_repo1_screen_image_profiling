"""Helpers for processing pooled bulk (well-level) profiles in the style of the JUMP profiling recipe.

The recipe (https://github.com/broadinstitute/jump-profiling-recipe) processes the compound screens as

    pool wells -> drop NaN/inf features -> per-plate control stats -> variant filter ->
    MAD normalize (per plate) -> rank inverse normal transform -> pooled feature selection

We add one step the recipe doesn't have, between dropping non-finite features and the variant filter: a
correction for a plate-wide, left-to-right gradient found in both screens (see `correct_column_gradient`).

These helpers implement the steps that pycytominer does not provide (pooling, dropping non-finite
features, the per-plate variant filter, the position correction) and thin wrappers that apply pycytominer's
`normalize` and `feature_select` with the settings used here. Nothing in this module is specific to one screen.
"""

import pathlib
import warnings
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from pycytominer import feature_select, normalize
from scipy.stats import median_abs_deviation

PLATE_COLUMN = "Metadata_Plate"
WELL_COLUMN = "Metadata_Well"
BATCH_ID_COLUMN = "Metadata_Batch_Id"
N_PLATE_COLUMNS = 24  # both screens use 384-well (16 x 24) plates


def pool_plate_profiles(
    paths: Iterable[pathlib.Path],
    plate_column: str = PLATE_COLUMN,
    replace_spaces_in_plate_names: bool = True,
) -> pd.DataFrame:
    """Stack the per-plate profile files (one parquet per plate) into one table, one row per well.

    Every file must have the same metadata columns. Feature columns may differ between plates: a feature that is
    missing from some plates is NaN there, is reported with a warning, and is dropped by `drop_nonfinite_features`.
    Spaces in plate names are replaced by underscores by default, so that "Assay Plate_1_3" and "Assay_Plate_1_3"
    (both spellings occur in CHP-134) are the same plate.
    """
    paths = sorted(pathlib.Path(path) for path in paths)
    if not paths:
        raise FileNotFoundError("no per-plate profile files to pool")

    frames = [pd.read_parquet(path) for path in paths]
    metadata = [c for c in frames[0].columns if c.startswith("Metadata_")]
    for path, frame in zip(paths, frames):
        if {c for c in frame.columns if c.startswith("Metadata_")} != set(metadata):
            raise ValueError(f"{path.name} has different metadata columns than {paths[0].name}")

    in_every_plate = set.intersection(*(set(frame.columns) for frame in frames))
    in_any_plate = set().union(*(set(frame.columns) for frame in frames))
    if in_any_plate != in_every_plate:
        warnings.warn(
            f"{len(in_any_plate - in_every_plate)} columns are missing from some plates and will be NaN there: "
            f"{sorted(in_any_plate - in_every_plate)[:5]}"
        )

    pooled = pd.concat(frames, ignore_index=True)
    if replace_spaces_in_plate_names:
        pooled[plate_column] = pooled[plate_column].str.replace(" ", "_", regex=False)
    return pooled


def split_feature_columns(profiles: pd.DataFrame, features: Sequence[str]) -> list[str]:
    """The metadata columns: everything that is not one of ``features``, in their original order."""
    feature_set = set(features)
    return [column for column in profiles.columns if column not in feature_set]


def drop_nonfinite_features(profiles: pd.DataFrame, features: Sequence[str]) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Drop every feature that has a NaN or an infinite value in any well (as the JUMP recipe does).

    Returns the table without those columns, the features that remain and the features that were dropped.
    """
    features = list(features)
    finite = np.isfinite(profiles[features].to_numpy(dtype=float)).all(axis=0)
    kept = [feature for feature, keep in zip(features, finite) if keep]
    dropped = [feature for feature, keep in zip(features, finite) if not keep]
    return profiles.drop(columns=dropped), kept, dropped


def control_mask(profiles: pd.DataFrame, control_query: str) -> np.ndarray:
    """Boolean array marking the control wells that match a pandas query (the same kind pycytominer takes as ``samples``)."""
    return profiles.index.isin(profiles.query(control_query).index)


def drop_wells(profiles: pd.DataFrame, query: str) -> tuple[pd.DataFrame, int]:
    """Remove the wells that match a pandas query; returns the remaining wells and how many were removed."""
    drop = control_mask(profiles, query)
    return profiles.loc[~drop].reset_index(drop=True), int(drop.sum())


# ---- position correction ----
# Every plate has the same left-to-right gradient: DMSO controls differ between plate columns 2 and 24 (where
# they sit), and compound wells drift smoothly between the two across the columns in between. The gradient shows
# up in the controls, which have no compound at all, so it looks like a plate-wide technical effect rather than
# biology. These functions estimate that gradient from the compound wells (so it doesn't just fit noise from the
# ~30 controls per plate) and subtract it from every well before the variant filter and MAD normalization, so
# those steps see a plate baseline that isn't stretched by the two control columns sitting at opposite ends of it.


def split_well(wells: pd.Series) -> pd.DataFrame:
    """Row letter and column number of well names such as A01."""
    parts = wells.str.extract(r"^([A-Za-z]+)(\d+)$")
    return pd.DataFrame({"well_row": parts[0], "well_column": pd.to_numeric(parts[1])})


def extract_compound_id(batch_ids: pd.Series, pattern: str = r"^(BRD-[A-Za-z]\d+)") -> pd.Series:
    """The compound identifier inside a batch identifier (e.g. BRD-K75699339 from BRD-K75699339-001-12-9)."""
    return batch_ids.str.extract(pattern)[0]


def platemap_groups(
    profiles: pd.DataFrame, plate_column: str, compound_column: str, min_jaccard: float = 0.9
) -> pd.Series:
    """Label each well with its plate map: the group of plates that carry (nearly) the same compounds.

    Plates are compared by the set of compounds they hold (``compound_column``, missing for control and empty
    wells); plates whose sets overlap by at least ``min_jaccard`` (intersection over union) are joined, and the
    connected groups are the plate maps. A plate with no compounds forms a group of its own. Groups are named
    ``platemap_01``, ``platemap_02``, ... in order of first appearance.
    """
    sets = profiles.loc[profiles[compound_column].notna()].groupby(plate_column)[compound_column].agg(frozenset)
    plates = list(pd.unique(profiles[plate_column]))
    parent = {plate: plate for plate in plates}

    def find(plate):
        while parent[plate] != plate:
            parent[plate] = parent[parent[plate]]
            plate = parent[plate]
        return plate

    held = [plate for plate in plates if plate in sets.index]
    for i, first in enumerate(held):
        for second in held[i + 1 :]:
            union = len(sets[first] | sets[second])
            if union and len(sets[first] & sets[second]) / union >= min_jaccard:
                parent[find(second)] = find(first)
    roots = list(dict.fromkeys(find(plate) for plate in plates))
    names = {root: f"platemap_{i + 1:02d}" for i, root in enumerate(roots)}
    return profiles[plate_column].map({plate: names[find(plate)] for plate in plates})


def column_curve(
    profiles: pd.DataFrame,
    features: Sequence[str],
    compound_mask: np.ndarray,
    group: np.ndarray,
    well_column: str = WELL_COLUMN,
    n_columns: int = N_PLATE_COLUMNS,
) -> pd.DataFrame:
    """Each feature's plate-column gradient: for every column, the median across plate-map groups of that
    column's mean value among the compound wells of that group (a table of columns 1..``n_columns`` x features).

    Using the median across groups, rather than pooling every compound well together, keeps one unusually
    represented layout from dominating the curve. Columns with no compound wells (here, the two control columns)
    reuse the nearest compound column's value, since the curve can't be estimated there directly. The curve is
    centered so its average over columns is 0: it describes the gradient's shape, not its overall level.
    """
    features = list(features)
    col = split_well(profiles[well_column])["well_column"].to_numpy()
    index = pd.MultiIndex.from_arrays([group[compound_mask], col[compound_mask]], names=["group", "col"])
    group_column_means = (
        pd.DataFrame(profiles.loc[compound_mask, features].to_numpy(dtype=float), index=index, columns=features)
        .groupby(level=["group", "col"])
        .mean()
    )
    curve = group_column_means.groupby(level="col").median()
    curve = curve.reindex(range(1, n_columns + 1)).ffill().bfill()
    return curve - curve.mean(axis=0)


def correct_column_gradient(
    profiles: pd.DataFrame,
    features: Sequence[str],
    plate_column: str = PLATE_COLUMN,
    well_column: str = WELL_COLUMN,
    batch_id_column: str = BATCH_ID_COLUMN,
    min_jaccard: float = 0.9,
) -> pd.DataFrame:
    """Subtract each feature's plate-column gradient (see `column_curve`) from every well, controls included.

    The gradient is estimated once per plate-map group (plates sharing the same compounds; see `platemap_groups`)
    from that group's compound wells, then applied to every well at that column on every plate, in every group.
    """
    features = list(features)
    compound_id = extract_compound_id(profiles[batch_id_column])
    compound_mask = compound_id.notna().to_numpy()
    group = platemap_groups(
        profiles.assign(**{"_position_correction_cid": compound_id}), plate_column, "_position_correction_cid", min_jaccard
    ).to_numpy()

    curve = column_curve(profiles, features, compound_mask, group, well_column)
    col = split_well(profiles[well_column])["well_column"].to_numpy()
    corrected = profiles.copy()
    corrected[features] = profiles[features].to_numpy(dtype=float) - curve.reindex(col).to_numpy()
    return corrected


def plate_control_stats(
    profiles: pd.DataFrame,
    features: Sequence[str],
    is_control: np.ndarray,
    plate_column: str = PLATE_COLUMN,
    min_controls: int = 3,
) -> dict[str, pd.DataFrame]:
    """Median, MAD and |MAD / median| of every feature over each plate's control wells (tables of plates x features).

    The MAD is the raw median absolute deviation (no 1.4826 factor). ``abs_coef_var`` is 0 wherever the median is 0
    or the ratio is not finite, so a feature with a zero median never counts as variable.
    """
    features = list(features)
    controls = profiles.loc[is_control]
    counts = controls.groupby(plate_column).size().reindex(profiles[plate_column].unique(), fill_value=0)
    too_few = counts[counts < min_controls]
    if len(too_few):
        raise ValueError(f"plates with fewer than {min_controls} control wells: {too_few.to_dict()}")

    medians, mads = {}, {}
    for plate, plate_controls in controls.groupby(plate_column, sort=False):
        values = plate_controls[features].to_numpy(dtype=float)
        medians[plate] = np.nanmedian(values, axis=0)
        mads[plate] = median_abs_deviation(values, axis=0, nan_policy="omit")
    median = pd.DataFrame.from_dict(medians, orient="index", columns=features)
    mad = pd.DataFrame.from_dict(mads, orient="index", columns=features)

    with np.errstate(divide="ignore", invalid="ignore"):
        abs_coef_var = (mad / median).abs()
    abs_coef_var = abs_coef_var.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return {"median": median, "mad": mad, "abs_coef_var": abs_coef_var}


def variant_features(stats: dict[str, pd.DataFrame], min_abs_coef_var: float = 1e-3) -> tuple[list[str], list[str]]:
    """Features that vary in the control wells of every plate: MAD is not 0 and |MAD / median| is above ``min_abs_coef_var``.

    A feature that fails on even one plate is dropped. Returns (kept features, dropped features), in column order.
    """
    passes = (stats["mad"] != 0) & (stats["abs_coef_var"] > min_abs_coef_var)
    keep = passes.all(axis=0)
    return keep.index[keep].tolist(), keep.index[~keep].tolist()


def normalize_by_plate(
    profiles: pd.DataFrame,
    features: Sequence[str],
    control_query: str,
    plate_column: str = PLATE_COLUMN,
) -> pd.DataFrame:
    """MAD-normalize every plate on its own control wells with pycytominer (median and MAD from that plate's controls)."""
    features = list(features)
    meta_features = split_feature_columns(profiles, features)

    normalized = []
    for _, plate_profiles in profiles.groupby(plate_column, sort=False):
        normalized.append(
            normalize(
                profiles=plate_profiles.reset_index(drop=True),
                features=features,
                meta_features=meta_features,
                samples=control_query,
                method="mad_robustize",
            )
        )
    return pd.concat(normalized, ignore_index=True)


def inverse_normal_transform(profiles: pd.DataFrame, features: Sequence[str], random_state: int = 0) -> pd.DataFrame:
    """Rank-based inverse normal transform of every feature over all wells, with pycytominer's quantile method.

    pycytominer builds the quantile landmarks from at most 10,000 randomly chosen rows (scikit-learn's default), which
    matters for screens with more wells than that; seeding numpy's global generator makes that choice reproducible.
    """
    features = list(features)
    np.random.seed(random_state)
    return normalize(
        profiles=profiles,
        features=features,
        meta_features=split_feature_columns(profiles, features),
        samples="all",
        method="inverse_normal",
    )


def select_features(
    profiles: pd.DataFrame,
    features: Sequence[str],
    correlation_threshold: float = 0.9,
    freq_cut: float = 0.05,
    unique_cut: float = 0.01,
) -> pd.DataFrame:
    """Pooled feature selection with pycytominer, in the order the JUMP recipe uses.

    (1) frequency filter (the recipe's "variance" filter), (2) correlation filter, (3) pycytominer's default blocklist
    (the same list the recipe ships), (4) drop any feature with a missing value.
    """
    return feature_select(
        profiles=profiles,
        features=list(features),
        operation=["frequency_threshold", "correlation_threshold", "blocklist", "drop_na_columns"],
        freq_cut=freq_cut,
        unique_cut=unique_cut,
        corr_threshold=correlation_threshold,
        na_cutoff=0,
    )


def spherize(
    profiles: pd.DataFrame,
    control_query: str,
    epsilon: float = 1e-6,
    method: str = "ZCA-cor",
    freq_cut: float = 0.05,
    unique_cut: float = 0.01,
) -> pd.DataFrame:
    """Whiten the pooled profiles with a transform fitted on the pooled control wells (the sphering step used so far).

    Features with too little variation among the controls are removed first, since they cannot be whitened.
    """
    varying = feature_select(
        profiles=profiles,
        operation="frequency_threshold",
        freq_cut=freq_cut,
        unique_cut=unique_cut,
        samples=control_query,
    )
    return normalize(
        profiles=varying,
        method="spherize",
        samples=control_query,
        spherize_center=True,
        spherize_method=method,
        spherize_epsilon=epsilon,
    )
