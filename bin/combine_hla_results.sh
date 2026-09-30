#!/usr/bin/env bash
set -euo pipefail

# INPUT:
#   Immediate subfolders under INPUT_BASE.
#
# MHCflurry files:
#   mhcflurry_*.out
#
# Each MHCflurry line must contain exactly:
#   peptide allele predicted_binding_affinity affinity_percentile_rank
#
# netMHC files:
#   ALL files whose basename starts with "netMHC" and ends EXACTLY with "out"
#
#
# OUTPUT:
#   One row per UNIQUE MHCflurry peptide.
#   For each HLA allele, columns contain both MHCflurry outputs and all
#   standard netMHCpan fields:
#     Pos, Core, Of, Gp, Gl, Ip, Il, Icore, Identity,
#     Score_EL, %Rank_EL, BindLevel
#


INPUT_BASE="hla_pred_out"
OUTPUT_DIR=""
OUTPUT_NAME="combined_hla_predictions_long.csv"

usage() {
    cat <<EOF
Usage:
  $0 [-i INPUT_DIR] [-o OUTPUT_DIR] [-n OUTPUT_NAME]

Options:
  -i INPUT_DIR      Directory containing prediction subfolders.
                    Default:
                    hla_pred_out

  -o OUTPUT_DIR     Directory for final CSV files.
                    If omitted, each CSV goes in its source folder.

  -n OUTPUT_NAME    Base output filename.
                    Default: combined_hla_predictions.csv

  -h                Show this help.

Notes:
  * All files starting with "netMHC" and ending EXACTLY with "out" are combined.
  * Existing final CSV files are skipped.
  * Only peptides with MHCflurry percentile rank < 3 OR netMHCpan
    %Rank_EL < 3 OR %Rank_BA < 3 (for any allele) are written.
  * Intermediate data are stored on disk to avoid excessive RAM use.
EOF
}

while getopts ":i:o:n:h" opt; do
    case "$opt" in
        i) INPUT_BASE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        n) OUTPUT_NAME="$OPTARG" ;;
        h)
            usage
            exit 0
            ;;
        :)
            echo "ERROR: -$OPTARG requires an argument." >&2
            usage >&2
            exit 1
            ;;
        \?)
            echo "ERROR: unknown option -$OPTARG" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -d "$INPUT_BASE" ]]; then
    echo "ERROR: input directory does not exist: $INPUT_BASE" >&2
    exit 1
fi

if [[ "$OUTPUT_NAME" == */* ]]; then
    echo "ERROR: -n must be a filename only. Use -o for the directory." >&2
    exit 1
fi

INPUT_BASE="$(cd "$INPUT_BASE" && pwd)"

if [[ -n "$OUTPUT_DIR" ]]; then
    mkdir -p "$OUTPUT_DIR"
    OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"
fi

python3 - "$INPUT_BASE" "$OUTPUT_DIR" "$OUTPUT_NAME" <<'PY'
import csv
import glob
import os
import re
import sqlite3
import sys

INPUT_BASE = sys.argv[1]
OUTPUT_DIR = sys.argv[2]
OUTPUT_NAME = sys.argv[3]

NET_RE = re.compile(r"^netMHC.*out$")


def normalize_allele(value):
    value = value.strip().upper().replace(" ", "")

    x = value[4:] if value.startswith("HLA-") else value
    m = re.fullmatch(r"([A-Z0-9]+)\*?(\d{2}):?(\d{2})", x)

    if m:
        return f"HLA-{m.group(1)}*{m.group(2)}:{m.group(3)}"

    return value


def allele_column_prefix(allele):
    x = allele.replace("HLA-", "HLA_")
    x = x.replace("*", "_")
    x = x.replace(":", "_")
    x = x.replace("-", "_")
    return x


def mhcflurry_sort_key(path):
    name = os.path.basename(path)
    m = re.search(r"mhcflurry_(\d+)\.out$", name)
    if m:
        return int(m.group(1))
    return 10**18


def find_netmhc_files(folder):
    result = []

    for path in glob.glob(os.path.join(folder, "netMHC*")):
        if (
            os.path.isfile(path)
            and NET_RE.fullmatch(os.path.basename(path))
        ):
            result.append(path)

    return sorted(result)


def output_path(folder):
    folder_name = os.path.basename(folder.rstrip(os.sep))

    if OUTPUT_DIR:
        return os.path.join(
            OUTPUT_DIR,
            f"{folder_name}_{OUTPUT_NAME}"
        )

    return os.path.join(folder, OUTPUT_NAME)


def db_path_for(folder, out_path):
    folder_name = os.path.basename(folder.rstrip(os.sep))
    parent = os.path.dirname(out_path)
    return os.path.join(parent, f".{folder_name}.hla_combine.sqlite")


def create_db(db_path):
    if os.path.exists(db_path):
        os.remove(db_path)

    con = sqlite3.connect(db_path)

    # Fast, disk-backed settings suitable for a disposable intermediate DB.
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA temp_store=FILE")
    con.execute("PRAGMA cache_size=-65536")  # about 64 MB cache

    con.executescript(
        """
        CREATE TABLE peptides (
            ord INTEGER PRIMARY KEY AUTOINCREMENT,
            peptide TEXT NOT NULL UNIQUE
        );

        CREATE TABLE alleles (
            allele TEXT PRIMARY KEY
        );

        CREATE TABLE mhc (
            peptide TEXT NOT NULL,
            allele TEXT NOT NULL,
            affinity TEXT,
            percentile TEXT,
            PRIMARY KEY (peptide, allele)
        ) WITHOUT ROWID;

        CREATE TABLE net (
            peptide TEXT NOT NULL,
            allele TEXT NOT NULL,
            pos TEXT,
            core TEXT,
            of_value TEXT,
            gp TEXT,
            gl TEXT,
            ip TEXT,
            il TEXT,
            icore TEXT,
            identity TEXT,
            score_el TEXT,
            rank_el TEXT,
            score_ba TEXT,
            rank_ba TEXT,
            aff_nm TEXT,
            bind_level TEXT,
            PRIMARY KEY (peptide, allele)
        ) WITHOUT ROWID;

        CREATE INDEX idx_mhc_peptide ON mhc(peptide);
        CREATE INDEX idx_net_peptide ON net(peptide);
        """
    )

    return con


def parse_mhcflurry(con, files):
    total = 0
    duplicates = 0

    insert_peptide = """
        INSERT OR IGNORE INTO peptides(peptide)
        VALUES (?)
    """
    insert_allele = """
        INSERT OR IGNORE INTO alleles(allele)
        VALUES (?)
    """
    insert_mhc = """
        INSERT OR IGNORE INTO mhc(
            peptide, allele, affinity, percentile
        )
        VALUES (?, ?, ?, ?)
    """

    for file_i, path in enumerate(
        sorted(files, key=mhcflurry_sort_key), 1
    ):
        print(
            f"        MHCflurry {file_i}/{len(files)}: "
            f"{os.path.basename(path)}",
            flush=True,
        )

        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line_number, line in enumerate(fh, 1):
                s = line.strip()

                if not s:
                    continue

                fields = s.split()

                if len(fields) != 4:
                    raise ValueError(
                        f"{path}:{line_number}: expected exactly 4 columns, "
                        f"found {len(fields)}"
                    )

                peptide, allele, affinity, percentile = fields
                allele = normalize_allele(allele)

                # Validate numeric fields.
                float(affinity)
                float(percentile)

                con.execute(insert_peptide, (peptide,))
                con.execute(insert_allele, (allele,))

                cur = con.execute(
                    insert_mhc,
                    (peptide, allele, affinity, percentile),
                )

                if cur.rowcount == 0:
                    duplicates += 1

                total += 1

        # Commit file-by-file so the transaction does not become enormous.
        con.commit()

    return total, duplicates


def parse_netmhc(con, files):
    """
    Parse netMHCpan output using the TABLE HEADER rather than fixed column
    positions. This supports both:

      EL-only:
        Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity
        Score_EL %Rank_EL BindLevel

      EL + BA:
        Pos MHC Peptide Core Of Gp Gl Ip Il Icore Identity
        Score_EL %Rank_EL Score_BA %Rank_BA Aff(nM) BindLevel

    BindLevel may be absent for a non-binder or may contain two tokens
    such as '<= SB' or '<= WB'.
    """
    total = 0
    duplicates = 0
    files_with_ba = 0
    files_without_ba = 0

    insert_allele = """
        INSERT OR IGNORE INTO alleles(allele)
        VALUES (?)
    """

    insert_net = """
        INSERT OR IGNORE INTO net(
            peptide, allele,
            pos, core, of_value, gp, gl, ip, il, icore,
            identity,
            score_el, rank_el,
            score_ba, rank_ba, aff_nm,
            bind_level
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    header_to_internal = {
        "Pos": "pos",
        "MHC": "mhc",
        "Peptide": "peptide",
        "Core": "core",
        "Of": "of_value",
        "Gp": "gp",
        "Gl": "gl",
        "Ip": "ip",
        "Il": "il",
        "Icore": "icore",
        "Identity": "identity",
        "Score_EL": "score_el",
        "%Rank_EL": "rank_el",
        "Score_BA": "score_ba",
        "%Rank_BA": "rank_ba",
        "Aff(nM)": "aff_nm",
        "BindLevel": "bind_level",
    }

    for file_i, path in enumerate(files, 1):
        if file_i == 1 or file_i % 25 == 0 or file_i == len(files):
            print(
                f"        netMHC {file_i}/{len(files)}: "
                f"{os.path.basename(path)}",
                flush=True,
            )

        header = None
        this_file_has_ba = False

        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line_number, line in enumerate(fh, 1):
                fields = line.split()

                if not fields:
                    continue

                # Detect the prediction-table header.
                if (
                    len(fields) >= 3
                    and fields[0] == "Pos"
                    and fields[1] == "MHC"
                    and fields[2] == "Peptide"
                ):
                    header = fields
                    this_file_has_ba = (
                        "Score_BA" in header
                        or "%Rank_BA" in header
                        or "Aff(nM)" in header
                    )
                    continue

                if header is None:
                    continue

                # Prediction rows begin with an integer position and HLA allele.
                if (
                    not fields[0].isdigit()
                    or len(fields) < 3
                    or not fields[1].upper().startswith("HLA-")
                ):
                    continue

                # Everything before BindLevel is positional. BindLevel itself
                # can be blank or two tokens (e.g. '<= SB'), so join the rest.
                if "BindLevel" in header:
                    bind_idx = header.index("BindLevel")
                    fixed_header = header[:bind_idx]

                    if len(fields) < len(fixed_header):
                        continue

                    raw = dict(zip(fixed_header, fields[:len(fixed_header)]))
                    raw["BindLevel"] = " ".join(fields[len(fixed_header):])
                else:
                    if len(fields) < len(header):
                        continue
                    raw = dict(zip(header, fields[:len(header)]))

                row = {}
                for h, v in raw.items():
                    internal = header_to_internal.get(h)
                    if internal is not None:
                        row[internal] = v

                if "peptide" not in row or "mhc" not in row:
                    continue

                peptide = row["peptide"]
                allele = normalize_allele(row["mhc"])

                con.execute(insert_allele, (allele,))

                cur = con.execute(
                    insert_net,
                    (
                        peptide,
                        allele,
                        row.get("pos", ""),
                        row.get("core", ""),
                        row.get("of_value", ""),
                        row.get("gp", ""),
                        row.get("gl", ""),
                        row.get("ip", ""),
                        row.get("il", ""),
                        row.get("icore", ""),
                        row.get("identity", ""),
                        row.get("score_el", ""),
                        row.get("rank_el", ""),
                        row.get("score_ba", ""),
                        row.get("rank_ba", ""),
                        row.get("aff_nm", ""),
                        row.get("bind_level", ""),
                    ),
                )

                if cur.rowcount == 0:
                    duplicates += 1

                total += 1

        if header is not None:
            if this_file_has_ba:
                files_with_ba += 1
            else:
                files_without_ba += 1

        con.commit()

    return total, duplicates, files_with_ba, files_without_ba


def write_csv(con, out_path):
    """
    Write ONE row per unique peptide + HLA allele combination.

    A row is retained only if that SAME peptide+allele combination has:
        MHCflurry percentile rank < 3
        OR netMHCpan %Rank_EL < 3
        OR netMHCpan %Rank_BA < 3
    """
    fieldnames = [
        "peptide",
        "allele",
        "mhcflurry_affinity",
        "mhcflurry_percentile_rank",
        "netmhcpan_pos",
        "netmhcpan_core",
        "netmhcpan_of",
        "netmhcpan_gp",
        "netmhcpan_gl",
        "netmhcpan_ip",
        "netmhcpan_il",
        "netmhcpan_icore",
        "netmhcpan_identity",
        "netmhcpan_score_el",
        "netmhcpan_rank_el",
        "netmhcpan_score_ba",
        "netmhcpan_rank_ba",
        "netmhcpan_aff_nm",
        "netmhcpan_bind_level",
    ]

    query = """
        WITH keys AS (
            SELECT peptide, allele FROM mhc
            UNION
            SELECT peptide, allele FROM net
        )
        SELECT
            k.peptide,
            k.allele,
            m.affinity,
            m.percentile,
            n.pos,
            n.core,
            n.of_value,
            n.gp,
            n.gl,
            n.ip,
            n.il,
            n.icore,
            n.identity,
            n.score_el,
            n.rank_el,
            n.score_ba,
            n.rank_ba,
            n.aff_nm,
            n.bind_level
        FROM keys AS k
        LEFT JOIN mhc AS m
          ON m.peptide = k.peptide
         AND m.allele = k.allele
        LEFT JOIN net AS n
          ON n.peptide = k.peptide
         AND n.allele = k.allele
        WHERE
            (
                m.percentile IS NOT NULL
                AND TRIM(m.percentile) <> ''
                AND CAST(m.percentile AS REAL) < 3.0
            )
            OR
            (
                n.rank_el IS NOT NULL
                AND TRIM(n.rank_el) <> ''
                AND CAST(n.rank_el AS REAL) < 3.0
            )
            OR
            (
                n.rank_ba IS NOT NULL
                AND TRIM(n.rank_ba) <> ''
                AND CAST(n.rank_ba AS REAL) < 3.0
            )
        ORDER BY k.peptide, k.allele
    """

    tmp_path = out_path + ".tmp"
    row_count = 0

    try:
        with open(tmp_path, "w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out)
            writer.writerow(fieldnames)

            for row in con.execute(query):
                writer.writerow([
                    "" if value is None else value
                    for value in row
                ])
                row_count += 1

        os.replace(tmp_path, out_path)

    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise

    return row_count


def process_folder(folder):
    name = os.path.basename(folder.rstrip(os.sep))
    out_path = output_path(folder)

    # Skip immediately if the final result already exists.
    if os.path.isfile(out_path):
        return {
            "status": "exists",
            "output": out_path,
        }

    mhc_files = glob.glob(
        os.path.join(folder, "mhcflurry_*.out")
    )
    net_files = find_netmhc_files(folder)

    if not mhc_files or not net_files:
        return {
            "status": "missing",
            "mhc_files": len(mhc_files),
            "net_files": len(net_files),
        }

    db_path = db_path_for(folder, out_path)

    print(
        f"START {name}: "
        f"MHCflurry files={len(mhc_files)}, "
        f"netMHC files={len(net_files)}",
        flush=True,
    )

    con = None

    try:
        con = create_db(db_path)

        mhc_total, mhc_duplicates = parse_mhcflurry(
            con, mhc_files
        )

        (
            net_total,
            net_duplicates,
            net_files_with_ba,
            net_files_without_ba,
        ) = parse_netmhc(con, net_files)

        output_rows = write_csv(
            con, out_path
        )

        return {
            "status": "created",
            "output": out_path,
            "output_rows": output_rows,
            "mhc_rows": mhc_total,
            "net_rows": net_total,
            "mhc_duplicates": mhc_duplicates,
            "net_duplicates": net_duplicates,
            "mhc_files": len(mhc_files),
            "net_files": len(net_files),
            "net_files_with_ba": net_files_with_ba,
            "net_files_without_ba": net_files_without_ba,
        }

    finally:
        if con is not None:
            con.close()

        # Intermediate database is disposable.
        if os.path.exists(db_path):
            os.remove(db_path)


folders = sorted(
    p
    for p in glob.glob(os.path.join(INPUT_BASE, "*"))
    if os.path.isdir(p)
)

created = 0
existing = 0
missing = 0
failed = 0

for folder in folders:
    name = os.path.basename(folder)

    try:
        result = process_folder(folder)

        if result["status"] == "exists":
            print(f"SKIP  {name}: output already exists")
            print(f"      -> {result['output']}")
            existing += 1
            continue

        if result["status"] == "missing":
            print(
                f"SKIP  {name}: "
                f"MHCflurry files={result['mhc_files']}, "
                f"netMHC files={result['net_files']}"
            )
            missing += 1
            continue

        print(
            f"OK    {name}: "
            f"qualifying peptide+allele rows={result['output_rows']}, "
            f"MHCflurry rows read={result['mhc_rows']}, "
            f"netMHC rows read={result['net_rows']}"
        )

        print(
            f"      netMHC files with BA columns="
            f"{result['net_files_with_ba']}, "
            f"without BA columns={result['net_files_without_ba']}"
        )

        if result["mhc_duplicates"] or result["net_duplicates"]:
            print(
                f"      duplicate peptide+allele rows ignored: "
                f"MHCflurry={result['mhc_duplicates']}, "
                f"netMHC={result['net_duplicates']}"
            )

        print(f"      -> {result['output']}")
        created += 1

    except Exception as exc:
        print(f"ERROR {name}: {exc}", file=sys.stderr)
        failed += 1

print()
print(
    f"Finished: created={created}, "
    f"already existing={existing}, "
    f"missing inputs={missing}, "
    f"failed={failed}"
)

if failed:
    sys.exit(2)
PY
