#!/usr/bin/env python3
# =============================================================================
# 08_extract_osp_profiles.py
# Extract midazolam oral PK profiles from OSP observed-data database
# =============================================================================
#
# Input:  data/literature_osp_midazolam/raw/ObsDataPK_OSP.xlsx
# Output: data/literature_osp_midazolam/extracted/*.csv
#         data/literature_osp_midazolam/provenance.md
#         data/literature_osp_midazolam/transforms.md
#
# Usage:  python3 analysis/08_extract_osp_profiles.py

import os
import sys
import csv
import math
import shutil
from datetime import datetime

try:
    import openpyxl
except ImportError:
    print("ERROR: openpyxl required. Install with: pip install openpyxl")
    sys.exit(1)

# --- Paths ------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_ROOT = os.path.join(PROJECT_ROOT, "data", "literature_osp_midazolam")
RAW_DIR = os.path.join(DATA_ROOT, "raw")
OUT_DIR = os.path.join(DATA_ROOT, "extracted")
XLSX_PATH = os.path.join(RAW_DIR, "ObsDataPK_OSP.xlsx")

os.makedirs(OUT_DIR, exist_ok=True)

# --- Constants --------------------------------------------------------------
MW_MIDAZOLAM = 325.78  # g/mol

# Unit conversion factors to ng/mL
UNIT_FACTORS = {
    "ng/mL": 1.0,
    "ng/ml": 1.0,
    "µg/L":  1.0,        # 1 µg/L = 1 ng/mL
    "ug/L":  1.0,
    "mg/L":  1000.0,     # 1 mg/L = 1000 ng/mL
    "pg/mL": 0.001,      # 1 pg/mL = 0.001 ng/mL
    "nmol/l": MW_MIDAZOLAM / 1000.0,  # nmol/L * MW/1000 = ng/mL
    "nmol/L": MW_MIDAZOLAM / 1000.0,
    "ng/L":  0.001,      # 1 ng/L = 0.001 ng/mL
}

# --- Logging ----------------------------------------------------------------
log_lines = []


def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    log_lines.append(line)


# --- Load workbook ----------------------------------------------------------
log(f"Loading {XLSX_PATH}")
wb = openpyxl.load_workbook(XLSX_PATH, data_only=True)

# --- Parse Studies sheet ----------------------------------------------------
log("Parsing Studies sheet...")
ws_studies = wb["Studies"]
studies = {}
for row_idx in range(3, ws_studies.max_row + 1):
    sid = ws_studies.cell(row=row_idx, column=1).value
    if sid is None:
        continue
    studies[sid] = {
        "study_id": sid,
        "study": ws_studies.cell(row=row_idx, column=2).value,
        "reference": ws_studies.cell(row=row_idx, column=3).value,
        "grouping": ws_studies.cell(row=row_idx, column=4).value,
        "compound": ws_studies.cell(row=row_idx, column=5).value,
        "compartment": ws_studies.cell(row=row_idx, column=6).value,
        "data_type": ws_studies.cell(row=row_idx, column=7).value,
        "source": ws_studies.cell(row=row_idx, column=8).value,
        "dose_str": ws_studies.cell(row=row_idx, column=10).value,
        "dose_unit": ws_studies.cell(row=row_idx, column=12).value,
        "route": ws_studies.cell(row=row_idx, column=13).value,
        "fasted_fed": ws_studies.cell(row=row_idx, column=22).value,
        "species": ws_studies.cell(row=row_idx, column=26).value,
        "n": ws_studies.cell(row=row_idx, column=27).value,
        "n_female": ws_studies.cell(row=row_idx, column=28).value,
        "ethnicity": ws_studies.cell(row=row_idx, column=69).value,
    }

log(f"  Total studies parsed: {len(studies)}")

# --- Filter Studies: Midazolam, PO, Human -----------------------------------
mdz_study_ids = set()
for sid, s in studies.items():
    compound = str(s["compound"] or "").lower()
    route = str(s["route"] or "").lower()
    species = str(s["species"] or "").lower()
    if "midazolam" in compound and "po" in route and "human" in species:
        mdz_study_ids.add(sid)

log(f"  Midazolam PO Human studies: {len(mdz_study_ids)}")

# --- Parse dose strings -----------------------------------------------------
def parse_dose_mg(dose_str, dose_unit):
    """Parse dose string to numeric mg. Returns None if unparseable or mg/kg."""
    if dose_str is None:
        return None
    ds = str(dose_str).strip()
    du = str(dose_unit or "mg").strip().lower()

    # Exclude mg/kg doses
    if "mg/kg" in du or "mg/kg" in ds.lower():
        return None

    # Handle parenthetical notes like "15 (actually 7.5)"
    if "(" in ds:
        ds = ds.split("(")[0].strip()

    # Handle "X mg" embedded
    ds = ds.replace("mg", "").strip()

    # Handle ranges like "0.075 - 40"
    if "-" in ds and not ds.startswith("-"):
        parts = ds.split("-")
        try:
            return float(parts[0].strip())
        except ValueError:
            return None

    try:
        val = float(ds)
        if "µg" in str(dose_unit or "") or "ug" in str(dose_unit or ""):
            val = val / 1000.0
        return val
    except ValueError:
        return None


# --- Parse PK-Profiles sheet ------------------------------------------------
log("Parsing PK-Profiles sheet...")
ws_profiles = wb["PK-Profiles"]

profiles_raw = []
for row_idx in range(2, ws_profiles.max_row + 1):
    sid = ws_profiles.cell(row=row_idx, column=1).value  # ID
    if sid not in mdz_study_ids:
        continue

    analyte = ws_profiles.cell(row=row_idx, column=5).value
    if not analyte or "midazolam" not in str(analyte).lower():
        continue

    compartment = ws_profiles.cell(row=row_idx, column=6).value

    profiles_raw.append({
        "study_id": sid,
        "study_label": ws_profiles.cell(row=row_idx, column=2).value,
        "reference": ws_profiles.cell(row=row_idx, column=3).value,
        "grouping": ws_profiles.cell(row=row_idx, column=4).value,
        "analyte": analyte,
        "compartment": str(compartment or ""),
        "time": ws_profiles.cell(row=row_idx, column=7).value,
        "time_unit": ws_profiles.cell(row=row_idx, column=8).value,
        "avg": ws_profiles.cell(row=row_idx, column=9).value,
        "avg_unit": ws_profiles.cell(row=row_idx, column=10).value,
        "avg_type": ws_profiles.cell(row=row_idx, column=11).value,
        "var": ws_profiles.cell(row=row_idx, column=12).value,
        "var_unit": ws_profiles.cell(row=row_idx, column=13).value,
        "var_type": ws_profiles.cell(row=row_idx, column=14).value,
        "lloq": ws_profiles.cell(row=row_idx, column=15).value,
        "comment": ws_profiles.cell(row=row_idx, column=16).value,
    })

log(f"  Midazolam PO profile rows (all compartments): {len(profiles_raw)}")

# --- Filter to Plasma only --------------------------------------------------
profiles_plasma = [p for p in profiles_raw
                   if "plasma" in p["compartment"].lower()]
log(f"  Plasma only: {len(profiles_plasma)}")

# Exclude non-plasma (whole blood, urine, etc.)
excluded_compartments = set(p["compartment"] for p in profiles_raw
                           if "plasma" not in p["compartment"].lower())
if excluded_compartments:
    log(f"  Excluded compartments: {excluded_compartments}")

# --- Normalize units to ng/mL ----------------------------------------------
skipped_units = set()
profiles_norm = []
for p in profiles_plasma:
    unit = str(p["avg_unit"] or "").strip()
    factor = UNIT_FACTORS.get(unit)

    if factor is None:
        # Skip non-concentration units (%, µg mass, etc.)
        skipped_units.add(unit)
        continue

    avg_val = p["avg"]
    if avg_val is None or not isinstance(avg_val, (int, float)):
        continue

    conc_ngml = float(avg_val) * factor

    # Convert variability to SD in ng/mL
    sd_ngml = None
    var_val = p["var"]
    var_type = str(p["var_type"] or "").strip().lower()
    var_unit = str(p["var_unit"] or "").strip()
    var_factor = UNIT_FACTORS.get(var_unit, factor)

    study_info = studies.get(p["study_id"], {})
    n = study_info.get("n")
    if n is not None:
        try:
            n = int(n)
        except (ValueError, TypeError):
            n = None

    if var_val is not None and isinstance(var_val, (int, float)):
        var_converted = float(var_val) * var_factor
        if "sd" in var_type:
            sd_ngml = var_converted
        elif "sem" in var_type:
            if n and n > 0:
                sd_ngml = var_converted * math.sqrt(n)
        # Skip 95th CI and geom. SD for now (not directly convertible)

    # Parse dose
    dose_mg = parse_dose_mg(study_info.get("dose_str"), study_info.get("dose_unit"))

    time_val = p["time"]
    if time_val is None or not isinstance(time_val, (int, float)):
        continue

    profiles_norm.append({
        "study_id": p["study_id"],
        "study_label": p["study_label"],
        "grouping": p["grouping"] or "",
        "time_h": float(time_val),
        "conc_ngml": round(conc_ngml, 4),
        "sd_ngml": round(sd_ngml, 4) if sd_ngml is not None else "",
        "n": n or "",
        "dose_mg": dose_mg or "",
        "data_type": study_info.get("data_type", ""),
    })

if skipped_units:
    log(f"  Skipped non-concentration units: {skipped_units}")
log(f"  Normalized plasma profiles: {len(profiles_norm)}")

# --- Separate Heizmann 1983 individual curves --------------------------------
heizmann_rows = [p for p in profiles_norm
                 if "heizmann" in str(p["study_label"]).lower()
                 and p["data_type"] == "Individual"]
mean_rows = [p for p in profiles_norm
             if p not in heizmann_rows]

log(f"  Heizmann 1983 individual rows: {len(heizmann_rows)}")
log(f"  Mean profile rows: {len(mean_rows)}")

# --- Write mean_profiles_midazolam_po.csv -----------------------------------
mean_csv_path = os.path.join(OUT_DIR, "mean_profiles_midazolam_po.csv")
fieldnames_mean = ["study_id", "study_label", "grouping", "time_h",
                   "conc_ngml", "sd_ngml", "n", "dose_mg"]
with open(mean_csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames_mean)
    writer.writeheader()
    for row in sorted(mean_rows, key=lambda x: (x["study_id"], x["time_h"])):
        writer.writerow({k: row[k] for k in fieldnames_mean})

log(f"  Wrote {mean_csv_path} ({len(mean_rows)} rows)")

# --- Write study_table_midazolam_po.csv -------------------------------------
study_table_path = os.path.join(OUT_DIR, "study_table_midazolam_po.csv")

# Aggregate per study_id
study_agg = {}
for row in mean_rows:
    sid = row["study_id"]
    if sid not in study_agg:
        info = studies.get(sid, {})
        study_agg[sid] = {
            "study_id": sid,
            "study_label": row["study_label"],
            "reference": info.get("reference", ""),
            "grouping": row["grouping"],
            "n": info.get("n", ""),
            "dose_mg": row["dose_mg"],
            "n_timepoints": 0,
            "time_min_h": float("inf"),
            "time_max_h": float("-inf"),
            "data_type": info.get("data_type", ""),
            "fasted_fed": info.get("fasted_fed", ""),
        }
    study_agg[sid]["n_timepoints"] += 1
    t = row["time_h"]
    study_agg[sid]["time_min_h"] = min(study_agg[sid]["time_min_h"], t)
    study_agg[sid]["time_max_h"] = max(study_agg[sid]["time_max_h"], t)

study_table_fields = ["study_id", "study_label", "reference", "grouping", "n",
                      "dose_mg", "n_timepoints", "time_range_h", "data_type",
                      "fasted_fed"]
with open(study_table_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=study_table_fields)
    writer.writeheader()
    for sid in sorted(study_agg.keys()):
        s = study_agg[sid]
        s["time_range_h"] = f"{s['time_min_h']:.1f}-{s['time_max_h']:.1f}"
        writer.writerow({k: s[k] for k in study_table_fields})

log(f"  Wrote {study_table_path} ({len(study_agg)} studies)")

# --- Write heizmann_1983_individual.csv -------------------------------------
heiz_csv_path = os.path.join(OUT_DIR, "heizmann_1983_individual.csv")
heiz_fields = ["study_id", "study_label", "grouping", "time_h", "conc_ngml"]
with open(heiz_csv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=heiz_fields)
    writer.writeheader()
    for row in sorted(heizmann_rows, key=lambda x: (str(x["grouping"]), x["time_h"])):
        writer.writerow({k: row[k] for k in heiz_fields})

log(f"  Wrote {heiz_csv_path} ({len(heizmann_rows)} rows)")

# --- Write provenance.md ---------------------------------------------------
n_unique_studies = len(set(r["study_id"] for r in profiles_norm))
provenance_path = os.path.join(DATA_ROOT, "provenance.md")
with open(provenance_path, "w") as f:
    f.write("# Data Provenance\n\n")
    f.write("## Source\n\n")
    f.write("- **Database:** Open-Systems-Pharmacology/Database-for-observed-data\n")
    f.write("- **File:** ObsDataPK_OSP.xlsx\n")
    f.write(f"- **Download date:** {datetime.now().strftime('%Y-%m-%d')}\n")
    f.write("- **URL:** https://github.com/Open-Systems-Pharmacology/Database-for-observed-data\n\n")
    f.write("## Extraction Pipeline\n\n")
    f.write("Script: `analysis/08_extract_osp_profiles.py`\n\n")
    f.write("### Filtering Steps\n\n")
    f.write(f"| Step | Filter | Rows Remaining |\n")
    f.write(f"|------|--------|----------------|\n")
    f.write(f"| 0 | All PK-Profiles rows | {ws_profiles.max_row - 1} |\n")
    f.write(f"| 1 | Study ID in midazolam PO human studies | {len(profiles_raw)} |\n")
    f.write(f"| 2 | Compartment = Plasma | {len(profiles_plasma)} |\n")
    f.write(f"| 3 | Valid concentration units (ng/mL convertible) | {len(profiles_norm)} |\n")
    f.write(f"| 4a | Mean profiles (excluding Heizmann individual) | {len(mean_rows)} |\n")
    f.write(f"| 4b | Heizmann 1983 individual curves | {len(heizmann_rows)} |\n\n")
    f.write(f"### Summary\n\n")
    f.write(f"- Unique studies with extracted profiles: {n_unique_studies}\n")
    f.write(f"- Total mean profile data points: {len(mean_rows)}\n")
    f.write(f"- Heizmann 1983 individual data points: {len(heizmann_rows)}\n")
    if skipped_units:
        f.write(f"- Skipped non-concentration units: {', '.join(sorted(skipped_units))}\n")
    if excluded_compartments:
        f.write(f"- Excluded non-plasma compartments: {', '.join(sorted(excluded_compartments))}\n")

log(f"  Wrote {provenance_path}")

# --- Write transforms.md ---------------------------------------------------
transforms_path = os.path.join(DATA_ROOT, "transforms.md")
with open(transforms_path, "w") as f:
    f.write("# Unit Conversion and Data Transforms\n\n")
    f.write("## Concentration Unit Normalization\n\n")
    f.write("All concentrations are normalized to ng/mL using the following factors:\n\n")
    f.write("| Source Unit | Conversion Factor | Formula |\n")
    f.write("|------------|-------------------|----------|\n")
    for unit, factor in sorted(UNIT_FACTORS.items()):
        if unit == "nmol/l" or unit == "nmol/L":
            f.write(f"| {unit} | {factor:.4f} | conc_nmol/L * MW / 1000 (MW = {MW_MIDAZOLAM}) |\n")
        else:
            f.write(f"| {unit} | {factor} | conc * {factor} |\n")
    f.write(f"\nMW (Midazolam) = {MW_MIDAZOLAM} g/mol\n\n")
    f.write("## Variability Conversion\n\n")
    f.write("| Source VarType | Conversion to SD |\n")
    f.write("|---------------|------------------|\n")
    f.write("| arith. SD | Used directly (with unit conversion) |\n")
    f.write("| arith. SEM | SD = SEM * sqrt(N), where N from Studies sheet |\n")
    f.write("| geom. SD | Not converted (geometric SD not directly comparable) |\n")
    f.write("| 95th CI | Not converted |\n\n")
    f.write("## Dose Parsing\n\n")
    f.write("- Dose strings parsed from Studies sheet column 10\n")
    f.write("- Parenthetical notes stripped: '15 (actually 7.5)' → 15\n")
    f.write("- mg/kg doses excluded (cannot convert without individual body weight)\n")
    f.write("- Dose unit 'µg' converted to mg by dividing by 1000\n")

log(f"  Wrote {transforms_path}")
log("Extraction complete.")
