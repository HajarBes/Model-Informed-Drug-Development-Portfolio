from pathlib import Path
import re
import csv

try:
    from pypdf import PdfReader
except ImportError:
    raise SystemExit("Install dependency first: pip install pypdf")

REPO = Path(__file__).resolve().parents[1]

FILES = {
    "lumakras_label": REPO / "data/raw/pk/sotorasib_LUMAKRAS_FDA_label_2025.pdf",
    "nejm_cb100": REPO / "data/raw/clinical/NEJM_CodeBreaK100_sotorasib_NSCLC_2021.pdf",
    "nejm_cb300": REPO / "data/raw/clinical/NEJM_CodeBreaK300_sotorasib_panitumumab_mCRC_2023.pdf",
    "vectibix_label": REPO / "data/raw/pk/panitumumab_VECTIBIX_FDA_label.pdf",
}

def pdf_text(path: Path) -> str:
    reader = PdfReader(str(path))
    parts = []
    for page in reader.pages:
        t = page.extract_text() or ""
        parts.append(t)
    return "\n".join(parts)

def find_first(pattern: str, text: str, flags=0):
    m = re.search(pattern, text, flags)
    return m.group(1).strip() if m else None

def main():
    # --- Extract from LUMAKRAS label ---
    label = pdf_text(FILES["lumakras_label"])
    half_life_h = find_first(r"mean terminal elimination half-life is\s+(\d+\s*hours?)", label, re.I)

    # --- Extract from VECTIBIX label ---
    vect = pdf_text(FILES["vectibix_label"])
    panitumumab_dose = find_first(r"recommended dosage.*?is\s+(6\s*mg/kg)\s+every\s+(14)\s+days", vect, re.I | re.S)
    if panitumumab_dose:
        panitumumab_sched = f"{panitumumab_dose} IV q{find_first(r'every\s+(14)\s+days', vect, re.I)}d"
    else:
        panitumumab_sched = None
    panitumumab_dose = find_first(r"recommended dosage.*?is\s+(6\s*mg/kg)\s+every\s+(14)\s+days", vect, re.I | re.S)
    if panitumumab_dose:
    # keep it simple and stable for portfolio
        panitumumab_sched = "6 mg/kg IV q14d"
    else:
        panitumumab_sched = None

    # --- Extract from NEJM CodeBreaK 100 (NSCLC) ---
    cb100 = pdf_text(FILES["nejm_cb100"])
    orr_cb100 = find_first(r"Objective response\s+— %.*?\)\s+(\d+\.?\d*)", cb100, re.I)
    pfs_cb100 = find_first(r"median progression-free survival.*?was\s+(\d+\.?\d*)\s+months", cb100, re.I)

    # --- Extract from NEJM CodeBreaK 300 (mCRC combo) ---
    cb300 = pdf_text(FILES["nejm_cb300"])
    orr_cb300 = find_first(r"objective response.*?was\s+(\d+\.?\d*)%\s", cb300, re.I)

    out_dir = REPO / "data/processed"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / "clinical_targets.csv"

    rows = [
        {
            "indication": "NSCLC",
            "therapy": "sotorasib (mono)",
            "dose": "960 mg QD",
            "ORR_percent": orr_cb100,
            "median_PFS_months": pfs_cb100,
            "source_pdf": str(FILES["nejm_cb100"].relative_to(REPO)),
        },
        {
            "indication": "mCRC",
            "therapy": "sotorasib + panitumumab (combo)",
            "dose": "sotorasib 960 mg QD + panitumumab 6 mg/kg q2w",
            "ORR_percent": orr_cb300,
            "median_PFS_months": "5.6",  # from label table; we keep fixed here
            "source_pdf": str(FILES["lumakras_label"].relative_to(REPO)),
        },
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)

    print("Wrote:", out_csv)
    print("Half-life (from label):", half_life_h)
    print("Panitumumab schedule (from label):", panitumumab_sched)

if __name__ == "__main__":
    main()
