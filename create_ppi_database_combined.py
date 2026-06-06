#!/usr/bin/env python3
"""
Robust local ETL script for Protein-Protein Docking Benchmark v5.5 -> SQLite.

Uses local files only:
  - Table_BM5.5.xlsx
  - benchmark5.5/  (unzipped folder; can contain structures/ subfolder)

Matches ER diagram column names:
  - Protein.protein_sequence
  - Chain.chain_sequence

Run:
  python create_ppi_database_local_v5.py
or:
  python create_ppi_database_local_v5.py --xlsx Table_BM5.5.xlsx --pdb-dir benchmark5.5 --db protein_interaction_benchmark.db
"""

import argparse
import gzip
import re
import sqlite3
from pathlib import Path
import pandas as pd

TYPE_DESC = {
    "EI": "Enzyme-Inhibitor",
    "ES": "Enzyme-Substrate",
    "ER": "Enzyme complex with a regulatory or accessory chain",
    "AA": "Antibody-Antigen",
    "AS": "Antigen - Single domain Antibody",
    "OG": "Others, G-protein containing",
    "OR": "Others, Receptor containing",
    "OX": "Others, miscellaneous",
}

DIFFICULTIES = {
    "rigid": "Rigid-body",
    "medium": "Medium Difficulty",
    "difficult": "Difficult",
}

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O", "MSE": "M", "ASX": "B", "GLX": "Z",
}


def clean_text(value):
    if pd.isna(value):
        return None
    text = str(value).strip()
    if text == "" or text.lower() == "nan":
        return None
    return text


def try_float(value):
    value = clean_text(value)
    if value is None:
        return None
    value = value.replace(",", "")
    value = re.sub(r"[^0-9eE+\-.]", "", value)
    try:
        return float(value)
    except Exception:
        return None


def parse_complex_id(complex_id):
    clean = str(complex_id).replace("*", "").strip()
    m = re.match(r"^([0-9A-Za-z]{4})_([^:]+):(.+)$", clean)
    if not m:
        return clean[:4].upper(), "", ""
    return m.group(1).upper(), m.group(2).replace("_", ""), m.group(3).replace("_", "")


def parse_structure_field(value):
    value = clean_text(value)
    if value is None:
        return None, ""
    s = value.replace("*", "").strip()
    m = re.match(r"^([0-9A-Za-z]{4})_?([^\s]*)", s)
    if not m:
        return s[:4].upper(), ""
    pdb_id = m.group(1).upper()
    chains = re.sub(r"\([^)]*\)", "", m.group(2)).replace("_", "")
    return pdb_id, chains


def infer_difficulty(row_index):
    if row_index < 162:
        return "rigid"
    if row_index < 222:
        return "medium"
    return "difficult"


def open_text_maybe_gz(path):
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "rt", errors="ignore")


def parse_atom_line_fallback(line):
    """Return (chain, residue_name, residue_id) from a PDB ATOM/HETATM line, allowing whitespace fallback."""
    # Fixed-column PDB first
    if len(line) >= 27:
        resname = line[17:20].strip().upper()
        chain = line[21].strip() or "_"
        resid = line[22:27].strip()
        if resname:
            return chain, resname, resid

    # Whitespace fallback for nonstandard cleaned files
    parts = line.split()
    # Typical: ATOM serial atom res chain resseq ...
    if len(parts) >= 6:
        resname = parts[3].upper()
        chain = parts[4] if len(parts[4]) == 1 else "_"
        resid = parts[5]
        return chain, resname, resid
    return "_", None, None


def pdb_sequences_from_file(pdb_file):
    """Extract chain sequences from SEQRES, falling back to unique ATOM/HETATM residues."""
    seqres = {}
    atom_order = []
    atom_seen = set()

    with open_text_maybe_gz(pdb_file) as handle:
        for line in handle:
            record = line[:6].strip().upper()
            if record == "SEQRES":
                chain = line[11].strip() or "_"
                residues = line[19:70].split()
                for res in residues:
                    aa = AA3.get(res.upper())
                    if aa:
                        seqres.setdefault(chain, []).append(aa)
            elif record in {"ATOM", "HETATM"}:
                chain, resname, resid = parse_atom_line_fallback(line)
                if resname not in AA3:
                    continue
                key = (chain, resid)
                if key not in atom_seen:
                    atom_seen.add(key)
                    atom_order.append((chain, AA3[resname]))

    if seqres:
        return {chain: "".join(residues) for chain, residues in seqres.items() if residues}

    by_chain = {}
    for chain, aa in atom_order:
        by_chain.setdefault(chain, []).append(aa)
    return {chain: "".join(residues) for chain, residues in by_chain.items() if residues}


def load_benchmark_records(xlsx_path):
    sheets = pd.read_excel(xlsx_path, sheet_name=None, header=None, dtype=str, engine="openpyxl")
    records = []
    current_difficulty = None
    complex_pattern = re.compile(r"^[0-9A-Za-z]{4}[_:][A-Za-z0-9_]+:[A-Za-z0-9_]+")

    for _, raw in sheets.items():
        for _, row in raw.iterrows():
            vals = [clean_text(x) for x in row.tolist()]
            row_text = " ".join(v for v in vals if v).lower()

            if "rigid" in row_text:
                current_difficulty = "rigid"
                continue
            if "medium" in row_text:
                current_difficulty = "medium"
                continue
            if "difficult" in row_text:
                current_difficulty = "difficult"
                continue

            hit = next((i for i, v in enumerate(vals) if v and complex_pattern.search(v.replace(" ", ""))), None)
            if hit is None:
                continue

            rec = vals[hit: hit + 11]
            while len(rec) < 11:
                rec.append(None)
            rec.append(current_difficulty)
            records.append(rec)
    return records


def possible_file_keys(path):
    stem = path.stem.upper()
    if stem.endswith(".PDB") or stem.endswith(".ENT"):
        stem = Path(stem).stem.upper()
    keys = {stem, stem[:4]}
    for token in re.split(r"[^0-9A-Z]+", stem):
        if len(token) == 4 and re.match(r"^[0-9A-Z]{4}$", token):
            keys.add(token)
    return keys


def build_sequence_cache(pdb_dir):
    seq_cache = {}
    pdb_files = []
    allowed_suffixes = {".pdb", ".ent", ".brk", ".pdb.gz", ".ent.gz"}
    for path in pdb_dir.rglob("*"):
        if not path.is_file():
            continue
        name = path.name.lower()
        if any(name.endswith(suf) for suf in allowed_suffixes):
            pdb_files.append(path)

    parsed_files = 0

    def add_to_cache(key, seqs):
        key = key.upper()
        entry = seq_cache.setdefault(key, {"chains": {}, "all_sequences": []})
        for chain, seq in seqs.items():
            if seq:
                entry["chains"].setdefault(chain, seq)
                entry["chains"].setdefault(chain.upper(), seq)
                entry["chains"].setdefault(chain.lower(), seq)
                if seq not in entry["all_sequences"]:
                    entry["all_sequences"].append(seq)

    for pdb_path in pdb_files:
        seqs = pdb_sequences_from_file(pdb_path)
        if not seqs:
            continue
        parsed_files += 1
        for key in possible_file_keys(pdb_path):
            add_to_cache(key, seqs)

    return seq_cache, len(pdb_files), parsed_files


def get_sequence(seq_cache, pdb_id, chain_id, chain_index=0):
    if pdb_id is None:
        return None
    entry = seq_cache.get(pdb_id.upper())
    if not entry:
        return None

    chains = entry["chains"]
    for candidate in [chain_id, str(chain_id).upper(), str(chain_id).lower(), "_"]:
        if candidate in chains:
            return chains[candidate]

    unique_sequences = entry["all_sequences"]
    if unique_sequences:
        return unique_sequences[min(chain_index, len(unique_sequences) - 1)]
    return None


def create_schema(conn):
    conn.executescript("""
    PRAGMA foreign_keys = ON;

    DROP TABLE IF EXISTS Complex_to_Chain;
    DROP TABLE IF EXISTS Interaction;
    DROP TABLE IF EXISTS Complex;
    DROP TABLE IF EXISTS Chain;
    DROP TABLE IF EXISTS Protein;
    DROP TABLE IF EXISTS Structure;
    DROP TABLE IF EXISTS Difficulty;
    DROP TABLE IF EXISTS ComplexType;

    CREATE TABLE ComplexType (
        complex_type_id INTEGER PRIMARY KEY,
        type_code TEXT NOT NULL UNIQUE,
        type_description TEXT NOT NULL
    );

    CREATE TABLE Difficulty (
        difficulty_id INTEGER PRIMARY KEY,
        difficulty_level TEXT NOT NULL UNIQUE,
        difficulty_description TEXT
    );

    CREATE TABLE Structure (
        structure_id INTEGER PRIMARY KEY,
        pdb_id TEXT NOT NULL,
        structure_type TEXT NOT NULL CHECK(structure_type IN ('bound', 'unbound')),
        UNIQUE(pdb_id, structure_type)
    );

    CREATE TABLE Protein (
        protein_id INTEGER PRIMARY KEY,
        protein_name TEXT NOT NULL,
        protein_sequence TEXT,
        UNIQUE(protein_name, protein_sequence)
    );

    CREATE TABLE Chain (
        chain_id INTEGER PRIMARY KEY,
        structure_id INTEGER NOT NULL,
        protein_id INTEGER NOT NULL,
        pdb_chain_id TEXT NOT NULL,
        chain_sequence TEXT,
        role TEXT,
        FOREIGN KEY (structure_id) REFERENCES Structure(structure_id),
        FOREIGN KEY (protein_id) REFERENCES Protein(protein_id),
        UNIQUE(structure_id, pdb_chain_id)
    );

    CREATE TABLE Complex (
        complex_id TEXT PRIMARY KEY,
        benchmark_name TEXT,
        complex_type_id INTEGER NOT NULL,
        difficulty_id INTEGER NOT NULL,
        bound_structure_id INTEGER NOT NULL,
        multimer_comment TEXT,
        FOREIGN KEY (complex_type_id) REFERENCES ComplexType(complex_type_id),
        FOREIGN KEY (difficulty_id) REFERENCES Difficulty(difficulty_id),
        FOREIGN KEY (bound_structure_id) REFERENCES Structure(structure_id)
    );

    CREATE TABLE Interaction (
        interaction_id INTEGER PRIMARY KEY,
        complex_id TEXT NOT NULL UNIQUE,
        interaction_type TEXT,
        rmsd REAL,
        hetatms_protein_1 TEXT,
        hetatms_protein_2 TEXT,
        interface_surface_area REAL,
        FOREIGN KEY (complex_id) REFERENCES Complex(complex_id)
    );

    CREATE TABLE Complex_to_Chain (
        complex_id TEXT NOT NULL,
        chain_id INTEGER NOT NULL,
        complex_chain_role TEXT,
        PRIMARY KEY (complex_id, chain_id),
        FOREIGN KEY (complex_id) REFERENCES Complex(complex_id),
        FOREIGN KEY (chain_id) REFERENCES Chain(chain_id)
    );
    """)



def clean_complex_id_for_update(x):
    return str(x).replace("*", "").replace(" ", "").strip() if x is not None else None

def looks_like_complex_for_update(value):
    if value is None:
        return False
    compact = str(value).replace(" ", "").replace("*", "").strip()
    return bool(re.match(r"^[0-9A-Za-z]{4}[_:][A-Za-z0-9_]+:[A-Za-z0-9_]+", compact))

def infer_interface_area_from_row(vals, hit):
    """Infer iSA/DASA from a benchmark row by choosing a plausible large numeric value."""
    candidates = []
    for offset, cell in enumerate(vals[hit + 1:], start=hit + 1):
        val = try_float(cell)
        if val is None:
            continue
        if val >= 100:
            candidates.append((offset, val, cell))
    if candidates:
        return max(candidates, key=lambda x: x[1])[1]

    fallback = []
    for offset, cell in enumerate(vals[hit + 1:], start=hit + 1):
        val = try_float(cell)
        if val is not None and val > 10:
            fallback.append((offset, val, cell))
    if fallback:
        return max(fallback, key=lambda x: x[1])[1]
    return None

def extract_interface_areas_from_excel(xlsx_path):
    """Extract interface surface area values from Table_BM5.5.xlsx using robust row scanning."""
    sheets = pd.read_excel(xlsx_path, sheet_name=None, header=None, dtype=str, engine="openpyxl")
    records = {}
    matched_rows = 0
    for _, raw in sheets.items():
        for _, row in raw.iterrows():
            vals = [clean_text(v) for v in row.tolist()]
            hit = next((i for i, v in enumerate(vals) if looks_like_complex_for_update(v)), None)
            if hit is None:
                continue
            complex_id = clean_complex_id_for_update(vals[hit])
            if not complex_id:
                continue
            matched_rows += 1
            isa = infer_interface_area_from_row(vals, hit)
            if isa is not None:
                records[complex_id] = isa
    return records, matched_rows

def get_sequence_with_fallback(seq_cache, keys, chain_id, chain_index=0):
    """Try multiple cache keys and then chain-specific/fallback sequence lookup."""
    seen = set()
    for key in keys:
        if not key:
            continue
        key = str(key).upper()
        if key in seen:
            continue
        seen.add(key)
        entry = seq_cache.get(key)
        if not entry:
            continue
        chains = entry.get("chains", {})
        for candidate in [chain_id, str(chain_id).upper(), str(chain_id).lower(), "_"]:
            if candidate in chains:
                return chains[candidate]
        all_sequences = entry.get("all_sequences", [])
        if all_sequences:
            return all_sequences[min(chain_index, len(all_sequences) - 1)]
    return None

def fix_sequences_after_load(conn, seq_cache):
    """Populate Chain.chain_sequence using role-aware benchmark filename fallbacks."""
    rows = conn.execute("""
        SELECT c.complex_id, s.pdb_id, ch.chain_id, ch.pdb_chain_id, ch.role, p.protein_id
        FROM Complex c
        JOIN Complex_to_Chain cc ON c.complex_id = cc.complex_id
        JOIN Chain ch ON cc.chain_id = ch.chain_id
        JOIN Structure s ON ch.structure_id = s.structure_id
        JOIN Protein p ON ch.protein_id = p.protein_id
        ORDER BY c.complex_id, ch.role, ch.pdb_chain_id
    """).fetchall()

    updated_chains = 0
    updated_proteins = 0
    missing = 0
    role_counter = {}

    for complex_id, structure_pdb, chain_id, pdb_chain_id, role, protein_id in rows:
        bound_pdb, _, _ = parse_complex_id(complex_id)
        if role == "protein_1_unbound":
            keys = [
                structure_pdb,
                f"{bound_pdb}_R_U", f"{bound_pdb}_RU", f"{bound_pdb}_R",
                f"{bound_pdb}_REC_U", f"{bound_pdb}_RECEPTOR_U",
                bound_pdb,
            ]
        elif role == "protein_2_unbound":
            keys = [
                structure_pdb,
                f"{bound_pdb}_L_U", f"{bound_pdb}_LU", f"{bound_pdb}_L",
                f"{bound_pdb}_LIG_U", f"{bound_pdb}_LIGAND_U",
                bound_pdb,
            ]
        else:
            keys = [structure_pdb, bound_pdb]

        role_counter[role] = role_counter.get(role, 0) + 1
        seq = get_sequence_with_fallback(seq_cache, keys, pdb_chain_id, role_counter[role] - 1)
        if not seq:
            missing += 1
            continue

        conn.execute("UPDATE Chain SET chain_sequence = ? WHERE chain_id = ?", (seq, chain_id))
        updated_chains += 1
        try:
            conn.execute(
                "UPDATE Protein SET protein_sequence = COALESCE(protein_sequence, ?) WHERE protein_id = ?",
                (seq, protein_id),
            )
            updated_proteins += 1
        except sqlite3.IntegrityError:
            pass

    conn.commit()
    return updated_chains, updated_proteins, missing

def fix_interface_area_after_load(conn, xlsx_path):
    """Populate Interaction.interface_surface_area from the Excel metadata table."""
    interface_records, matched_rows = extract_interface_areas_from_excel(xlsx_path)
    updated = 0
    for complex_id, isa in interface_records.items():
        cur = conn.execute(
            "UPDATE Interaction SET interface_surface_area = ? WHERE complex_id = ?",
            (isa, complex_id),
        )
        updated += cur.rowcount
    conn.commit()
    return matched_rows, len(interface_records), updated


def main():
    parser = argparse.ArgumentParser(description="Load Protein-Protein Docking Benchmark v5.5 into SQLite.")
    parser.add_argument("--xlsx", default="Table_BM5.5.xlsx", help="Path to downloaded Table_BM5.5.xlsx")
    parser.add_argument("--pdb-dir", default="benchmark5.5", help="Path to unzipped benchmark5.5 folder")
    parser.add_argument("--db", default="protein_interaction_benchmark.db", help="Output SQLite database path")
    args = parser.parse_args()

    xlsx_path = Path(args.xlsx)
    pdb_dir = Path(args.pdb_dir)
    db_path = Path(args.db)

    if not xlsx_path.exists():
        raise FileNotFoundError(f"Could not find Excel table: {xlsx_path.resolve()}")
    if not pdb_dir.exists() or not pdb_dir.is_dir():
        raise FileNotFoundError(f"Could not find unzipped PDB folder: {pdb_dir.resolve()}")

    records = load_benchmark_records(xlsx_path)
    if not records:
        raise ValueError("No benchmark records were found. Check that the Excel file is the BM5.5 table.")

    seq_cache, pdb_file_count, parsed_pdb_count = build_sequence_cache(pdb_dir)
    print(f"Found {len(records)} benchmark records in {xlsx_path}")
    print(f"Found {pdb_file_count} PDB-like files in {pdb_dir}")
    print(f"Parsed sequences from {parsed_pdb_count} PDB-like files")
    print(f"Built sequence cache for {len(seq_cache)} PDB/stem keys")
    if parsed_pdb_count == 0:
        print("WARNING: No sequences were parsed from PDB files. Check file extensions and PDB content.")

    conn = sqlite3.connect(db_path)
    create_schema(conn)

    for code, desc in TYPE_DESC.items():
        conn.execute("INSERT INTO ComplexType(type_code, type_description) VALUES (?, ?)", (code, desc))
    for level, desc in DIFFICULTIES.items():
        conn.execute("INSERT INTO Difficulty(difficulty_level, difficulty_description) VALUES (?, ?)", (level, desc))
    conn.commit()

    type_id = dict(conn.execute("SELECT type_code, complex_type_id FROM ComplexType"))
    diff_id = dict(conn.execute("SELECT difficulty_level, difficulty_id FROM Difficulty"))

    def upsert_structure(pdb_id, structure_type):
        if pdb_id is None:
            return None
        conn.execute("INSERT OR IGNORE INTO Structure(pdb_id, structure_type) VALUES (?, ?)", (pdb_id, structure_type))
        return conn.execute(
            "SELECT structure_id FROM Structure WHERE pdb_id = ? AND structure_type = ?",
            (pdb_id, structure_type),
        ).fetchone()[0]

    def upsert_protein(name, protein_sequence):
        name = clean_text(name) or "Unknown protein"
        conn.execute(
            "INSERT OR IGNORE INTO Protein(protein_name, protein_sequence) VALUES (?, ?)",
            (name, protein_sequence),
        )
        return conn.execute(
            """
            SELECT protein_id
            FROM Protein
            WHERE protein_name = ?
              AND COALESCE(protein_sequence, '') = COALESCE(?, '')
            """,
            (name, protein_sequence),
        ).fetchone()[0]

    def upsert_chain(structure_id, protein_id, pdb_chain_id, chain_sequence, role):
        conn.execute(
            """
            INSERT INTO Chain(structure_id, protein_id, pdb_chain_id, chain_sequence, role)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(structure_id, pdb_chain_id) DO UPDATE SET
                chain_sequence = COALESCE(excluded.chain_sequence, Chain.chain_sequence),
                role = COALESCE(Chain.role, excluded.role)
            """,
            (structure_id, protein_id, pdb_chain_id, chain_sequence, role),
        )
        return conn.execute(
            "SELECT chain_id FROM Chain WHERE structure_id = ? AND pdb_chain_id = ?",
            (structure_id, pdb_chain_id),
        ).fetchone()[0]

    loaded = 0
    skipped = 0
    missing_sequences = 0
    unmatched_pdb_ids = set()

    for idx, rec in enumerate(records):
        detected_difficulty = rec[11] if len(rec) > 11 else None
        rec = rec[:11] + [None] * max(0, 11 - len(rec))
        complex_id, cat, pdbid1, protein1, het1, pdbid2, protein2, het2, rmsd, interface_area, multimer = rec[:11]

        complex_id_clean = str(complex_id).replace("*", "").strip()
        cat = clean_text(cat)
        if cat not in type_id:
            skipped += 1
            continue

        difficulty = detected_difficulty or infer_difficulty(idx)
        bound_pdb, _, _ = parse_complex_id(complex_id_clean)
        p1_pdb, p1_chains = parse_structure_field(pdbid1)
        p2_pdb, p2_chains = parse_structure_field(pdbid2)

        bound_sid = upsert_structure(bound_pdb, "bound")
        p1_sid = upsert_structure(p1_pdb, "unbound")
        p2_sid = upsert_structure(p2_pdb, "unbound")

        conn.execute(
            """
            INSERT OR REPLACE INTO Complex(
                complex_id, benchmark_name, complex_type_id, difficulty_id, bound_structure_id, multimer_comment
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (complex_id_clean, complex_id_clean, type_id[cat], diff_id[difficulty], bound_sid, clean_text(multimer)),
        )

        conn.execute(
            """
            INSERT OR REPLACE INTO Interaction(
                complex_id, interaction_type, rmsd, hetatms_protein_1, hetatms_protein_2, interface_surface_area
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            (complex_id_clean, cat, try_float(rmsd), clean_text(het1), clean_text(het2), try_float(interface_area)),
        )

        for role, pdb_id, structure_id, protein_name, chain_string in [
            ("protein_1_unbound", p1_pdb, p1_sid, protein1, p1_chains),
            ("protein_2_unbound", p2_pdb, p2_sid, protein2, p2_chains),
        ]:
            if pdb_id is None or structure_id is None:
                continue
            chain_string = chain_string or "_"
            if pdb_id.upper() not in seq_cache:
                unmatched_pdb_ids.add(pdb_id.upper())

            for chain_index, chain_label in enumerate(chain_string):
                chain_sequence = get_sequence(seq_cache, pdb_id, chain_label, chain_index)
                if chain_sequence is None:
                    missing_sequences += 1
                protein_id = upsert_protein(protein_name, chain_sequence)
                db_chain_id = upsert_chain(structure_id, protein_id, chain_label, chain_sequence, role)
                conn.execute(
                    """
                    INSERT OR IGNORE INTO Complex_to_Chain(complex_id, chain_id, complex_chain_role)
                    VALUES (?, ?, ?)
                    """,
                    (complex_id_clean, db_chain_id, role),
                )

        loaded += 1

    conn.commit()

    # Post-load cleanup: robustly populate PDB chain sequences and interface surface area.
    # These steps are included here so the full ETL process runs from one script.
    seq_updated, prot_updated, seq_missing_after_fix = fix_sequences_after_load(conn, seq_cache)
    isa_rows_seen, isa_rows_extracted, isa_rows_updated = fix_interface_area_after_load(conn, xlsx_path)

    print(f"Created {db_path}")
    print(f"Post-load sequence rows updated: {seq_updated}")
    print(f"Post-load sequence lookups still missing: {seq_missing_after_fix}")
    print(f"Benchmark rows scanned for interface surface area: {isa_rows_seen}")
    print(f"Interface surface area values extracted: {isa_rows_extracted}")
    print(f"Interaction rows updated with interface surface area: {isa_rows_updated}")
    print(f"Loaded complexes: {loaded}")
    print(f"Skipped records: {skipped}")
    print(f"Requested chain records with missing extracted sequence: {missing_sequences}")
    if unmatched_pdb_ids:
        sample = ", ".join(sorted(list(unmatched_pdb_ids))[:20])
        print(f"PDB IDs from table not found in sequence cache: {len(unmatched_pdb_ids)}")
        print(f"Sample unmatched PDB IDs: {sample}")

    for table in ["ComplexType", "Difficulty", "Structure", "Protein", "Chain", "Complex", "Interaction", "Complex_to_Chain"]:
        count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
        print(f"{table}: {count}")

    chain_nulls = conn.execute("SELECT COUNT(*) FROM Chain WHERE chain_sequence IS NULL").fetchone()[0]
    protein_nulls = conn.execute("SELECT COUNT(*) FROM Protein WHERE protein_sequence IS NULL").fetchone()[0]
    chain_nonnulls = conn.execute("SELECT COUNT(*) FROM Chain WHERE chain_sequence IS NOT NULL").fetchone()[0]
    protein_nonnulls = conn.execute("SELECT COUNT(*) FROM Protein WHERE protein_sequence IS NOT NULL").fetchone()[0]
    print(f"Chain rows with sequence: {chain_nonnulls}")
    print(f"Chain rows still missing sequence: {chain_nulls}")
    print(f"Protein rows with sequence: {protein_nonnulls}")
    print(f"Protein rows still missing sequence: {protein_nulls}")
    isa_nonnulls = conn.execute("SELECT COUNT(*) FROM Interaction WHERE interface_surface_area IS NOT NULL").fetchone()[0]
    isa_nulls = conn.execute("SELECT COUNT(*) FROM Interaction WHERE interface_surface_area IS NULL").fetchone()[0]
    print(f"Interaction rows with interface_surface_area: {isa_nonnulls}")
    print(f"Interaction rows still missing interface_surface_area: {isa_nulls}")

    examples = conn.execute(
        """
        SELECT c.complex_id, us.pdb_id, ch.pdb_chain_id, LENGTH(ch.chain_sequence) AS seq_len,
               SUBSTR(ch.chain_sequence, 1, 20) AS seq_start
        FROM Complex c
        JOIN Complex_to_Chain cc ON c.complex_id = cc.complex_id
        JOIN Chain ch ON cc.chain_id = ch.chain_id
        JOIN Structure us ON ch.structure_id = us.structure_id
        WHERE ch.chain_sequence IS NOT NULL
        LIMIT 5
        """
    ).fetchall()
    print("Example non-NULL chain sequences:")
    for row in examples:
        print(row)

    conn.close()


if __name__ == "__main__":
    main()
