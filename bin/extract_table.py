#!/bin/env python3

import sys
import csv
import subprocess
import fnmatch
from pathlib import Path
from copy import deepcopy
from contextlib import ExitStack
import posixpath
from urllib.parse import unquote, urlsplit
import re
import os


FULL_MATCH = 0
ONE_MISMATCH = 1
TWO_MISMATCHES = 2

OUTPUT_DIR = Path(".")

BLAST_COLUMNS = [
    "qaccver", "qlen", "saccver", "slen", "sseq", "pident", "nident", "length",
    "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "staxids", "sscinames", "scomnames", "qcovs", "sskingdom",
    "sskingdoms", "stitle"
]

FULL_MATCH_COLUMNS = [
    "qaccver", "qlen", "saccver", "slen", "sseq", "pident", "nident", "length",
    "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue",
    "bitscore", "staxids", "sscinames", "scomnames", "qcovs", "sskingdom",
    "sskingdoms", "stitle", "organism_in_human", "microbe_location",
    "variant_type", "variant_tissue", "mutation", "HLA", "MHCf_rank", "net4_aff", "net4_rank"
]

MATCH_TABLE_COLUMNS = [
    "Peptide Sequence", 
    "Peptide length",
    "% Identity",
    "Query coverage",
    "Subject",
    "Peptide (BLAST results)",
    "MHC Class",
    "HLA",
    "Taxid",
    "Organism name",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Mutation",
    "MHCf_rank",
    "net4_aff",
    "net4_rank",
    "Antigen type",
    "Variant type",
    "Cancer",
    "Sample Tissue",
    "Gene",
    "WT Peptide",
    "AA Change",
    "Source",
    "Microbe location",
    "Title"
]

caatlas_dict = {
    8: "Breast Cancer",
    9: "Colon Carcinoma",
    10: "Glioblastoma",
    11: "Leukemia",
    12: "Lung Cancer",
    13: "Lymphoma",
    14: "Melanoma",
    15: "Meningioma"
}

ise_cancer_dict = {
    "AML":   "Acute Myeloid lLeukemia",
    "BALL":  "B-Acute Lymphoblastic Leukemia",
    "BRCA":  "Breast Cancer",
    "CCRC":  "Clear cell renal cell carcinoma",
    "CLL":   "Chronic lymphocytic leukemia",
    "COAD":  "Colon adenocarcinoma",
    "GBM":   "Glioblastoma",
    "MG":    "Meningioma",
    "NB":    "Neuroblastoma",
    "NHL":   "Non-Hodgkin lymphoma",
    "NSCLC": "Non-small cell lung cancer",
    "OV":    "Ovarian",
    "SKCM":  "Skin cutaneous melanoma",
    "TALL":  "T-cell acute lymphoblastic leukemia"
}

ise_tissue_dict = {
    "AML":   "bone marrow, blood",
    "BALL":  "bone marrow, blood",
    "BRCA":  "breast",
    "CCRC":  "kidney",
    "CLL":   "bone marrow, blood, lymph nodes, spleen",
    "COAD":  "colon",
    "GBM":   "brain",
    "MG":    "brain",
    "NB":    "sympathetic nervous system",
    "NHL":   "lymphatic tissue",
    "NSCLC": "lung" ,
    "OV":    "ovary",
    "SKCM":  "skin",
    "TALL":  "bone marow, blood" 
}

#    full_match_path = OUTPUT_DIR / "_full_matches.csv"
#    match_table_path = OUTPUT_DIR / "_match_tables.csv"

def get_filename_from_url(url):
    path = urlsplit(url).path
    decoded_path = unquote(path)
    return posixpath.basename(decoded_path)
    
def parse_pair(s):
    match = re.fullmatch(r'(\d+):(\d+)', s)
    return (int(match.group(1)), int(match.group(2))) if match else None
    
def main():
    if len(sys.argv) < 5:
        print(f"Usage: {sys.argv[0]} <prefix> <peptide_file> <taxid_file> <blast_output1> [blast_output2 ...]", file=sys.stderr)
        sys.exit(1)
     
    prefix        = sys.argv[1]
    peptide_file  = sys.argv[2]
    taxid_file    = sys.argv[3]
    blast_files   = sys.argv[4:]
    
    start_batch = 0
    end_batch = len(blast_files)
    batch_str = os.getenv("PBS_ARRAY_INDEX")
    if batch_str is not None:
        batch = int(batch_str)
        start_batch = (batch - 1) * 1000 if batch == 1 else ((batch - 1) * 1000) + 1
        start_batch = min(start_batch, len(blast_files))
        end_batch   = batch * 1000
        end_batch = min(end_batch, len(blast_files))
        prefix = prefix + "_" + str(start_batch) + "-" + str(end_batch)
    print("Start of batch: " + str(start_batch), flush = True)
    print("End of batch: " + str(end_batch), flush = True)
 
    pep_lookup = load_peptides(peptide_file)         # e.g. ref/SNV-derived.txt
    taxid_lookups = load_taxids(taxid_file)          # e.g. ref/HMRGD_txids_lineages2.csv

    file_paths = []
    file_paths.append(OUTPUT_DIR / (prefix + "_full_matches.csv"))
    for idx in range(len(taxid_lookups)):
        file_paths.append(OUTPUT_DIR / f"{prefix}_match_tables{idx}.csv")

    with ExitStack() as stack:
        files = [stack.enter_context(open(path, "w", newline="")) for path in file_paths]
        full_writer = csv.DictWriter(
            files[0],
            fieldnames=FULL_MATCH_COLUMNS,
            extrasaction="ignore"
        )
        full_writer.writeheader()

        table_writers = []
        for fidx in range(len(taxid_lookups)):
            table_writers.append(csv.DictWriter(
                files[fidx+1],
                fieldnames=MATCH_TABLE_COLUMNS,
                extrasaction="ignore"
            ))
            table_writers[fidx].writeheader()

        for blastp_output in blast_files[start_batch:end_batch]:
            full_matches = extract_full_matches(
                blastp_output,
                pep_lookup
            )
            full_writer.writerows(full_matches)

            match_tables = build_match_tables(
                full_matches,
                taxid_lookups
            )

            for fidx in range(len(match_tables)):
                table_writers[fidx].writerows(match_tables[fidx])

            # Free memory before next file
            del full_matches
            del match_tables

 
def load_peptides(peptide_file):
    """
    Read the peptide/SNV file once and store ALL rows by qaccver.
    Multiple rows with the same qaccver are stored as a list of dicts.
    """
    pep_lookup = {}
    peptide_base = get_filename_from_url(peptide_file)
    print(peptide_base, flush = True)

    with open(peptide_file, "r", encoding = "utf-8", newline="", errors="replace") as f:

        # Skip header
        next(f)

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= 2:
                parts = line.strip().split(",")
                    
            parts = [p.strip().strip('"') for p in parts]   # Remove surrounding quotes
            entry_added = False
            qaccver = ""
            
            if fnmatch.fnmatch(peptide_base, "*-derived.txt"):
                if len(parts) >= 13:
                    qaccver = parts[4]
                    entry = create_tsnadb_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "*tsnadb*"):
                if len(parts) >= 14:
                    qaccver = parts[8]
                    entry = create_tsnadb_valid_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "*neoantigens.txt"):
                if len(parts) >= 11: 
                    qaccver = parts[3]
                    entry = create_dbpepneo_entry(parts)
            elif peptide_base == "ITSNdb.csv":
                if len(parts) >= 11:
                    qaccver = parts[3]
                    entry = create_itsndb_entry(parts)
            elif peptide_base == "Valds.csv":
                if len(parts) >= 5:
                    qaccver = parts[1]
                    entry = create_itsndb_valids_entry(parts)
            elif peptide_base == "neodb_all.csv":
                if len(parts) >= 11:
                    qaccver = parts[1]
                    entry = create_neodb_entry(parts)
            elif peptide_base == "PTMAntigenList_SiteLevel.txt":
                if len(parts) >= 17:
                    qaccver = parts[0]
                    for idx in range(8,16):
                        if int(parts[idx]) > 0:
                            cancer_name = caatlas_dict.get(idx, "")
                            entry = create_caatlas_entry(parts, cancer_name)
                            pep_lookup.setdefault(qaccver, []).append(entry)
                            entry_added = True
                    if not entry_added:
                        entry = create_caatlas_entry(parts, "")
            elif peptide_base == "sequences_list.csv":
                if len(parts) >= 11:
                    qaccver = parts[1]
                    entry = create_caatlas_extra_entry(parts)
            elif peptide_base == "data_one_row_per_epitope.txt":
                if len(parts) >= 9:
                    qaccver = parts[3]
                    entry = create_cadbio_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "caped*.txt"):
                if len(parts) >= 9:
                    qaccver = parts[4]
                    entry = create_caped_entry(parts)
            elif peptide_base == "tcell_export.csv":
                if len(parts) >= 5:
                    qaccver = parts[0]
                    entry = create_cedar_entry(parts)
            elif peptide_base == "tcell_sub.csv":
                if len(parts) >= 8:
                    qaccver = parts[2]
                    entry = create_cedar_extra_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "*immunogenicNeo-peptideDataset.txt"):
                if len(parts) >= 14:
                    qaccver = parts[12]
                    entry = create_tumoragdb_entry(parts)
            elif peptide_base == "ISE_Normal.csv":
                if len(parts) >= 7:
                    qaccver = parts[1]
                    entry = create_ise_normal_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "ISE*.csv"):
                if len(parts) >= 5:
                    qaccver = parts[0]
                    cidx = peptide_base.find("_") + 1
                    cridx = peptide_base.rfind(".")
                    cancer_name = ise_cancer_dict.get(peptide_base[cidx:cridx], "")
                    tissue      = ise_tissue_dict.get(peptide_base[cidx:cridx], "")
                    entry = create_ise_entry(parts, cancer_name, tissue)
            elif peptide_base == "combined.csv":
                if len(parts) >= 18:
                    qaccver = parts[4]
                    entry = create_ligand_entry(parts)
            elif peptide_base == "tcellepitopes.csv":
                if len(parts) >= 3:
                    qaccver = parts[0]
                    entry = create_tantigen_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "bedran*.csv"):
                if len(parts) >= 6:
                    qaccver = parts[1]
                    entry = create_bedran_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "beijer*.csv"):
                if len(parts) >= 12:
                    qaccver = parts[2]
                    entry = create_beijer_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "Loffier.csv"):
                if len(parts) >= 4:
                    qaccver = parts[3].strip()
                    entry = create_loffier_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "Mizukoshi.csv"):
                if len(parts) >= 7:
                    qaccver = parts[4]
                    entry = create_mizukoshi_entry(parts)
            elif fnmatch.fnmatch(peptide_base, "reparaz*.csv"):
                if len(parts) >= 7:
                    qaccver = parts[3]
                    entry = create_reparaz_entry(parts)

            # Store ALL matches for this peptide
            if not entry_added:
                if qaccver == "":
                    print("Mismatched length:")
                    print(parts)
                else:
                    pep_lookup.setdefault(qaccver, []).append(entry)

    return pep_lookup

def create_tsnadb_entry(parts):
    mutation = parts[2].strip()
    mut_parts = mutation.split('_', 1)
    if len(mut_parts) == 2:
        gene = mut_parts[0]
        aa_change = mut_parts[1]
    else:
        gene = mutation
        aa_change = ""
 
    entry = {
        "source": "TSNAdb",
        "antigen_type": "neoantigen", 
        "variant_type": parts[0].strip(),
        "variant_tissue": parts[1].strip(),
        "HLA": parts[3].strip(),
        "peptide": parts[4].strip(),
        "MHCf_rank": parts[7].strip(),
        "net4_aff": parts[8].strip(),
        "net4_rank": parts[9].strip(),
        "gene": gene,
        "aa_change": aa_change
    }
    return entry
    
    
def create_tsnadb_valid_entry(parts):
    mutation = parts[6].strip()
    mut_parts = mutation.split('_', 1)
    if len(mut_parts) == 2:
        aa_change = mut_parts[1]
    else:
        aa_change = ""
 
    entry = {
        "source": "TSNAdb",
        "antigen_type": "neoantigen", 
        "variant_tissue": parts[2].strip(),
        "variant_type": parts[4].strip(),
        "gene": parts[5].strip(),
        "peptide": parts[8].strip(),
        "HLA": parts[9].strip(),
        "aa_change": aa_change
    }
    return entry
    
def create_dbpepneo_entry(parts):   # Done
    entry = {
        "source": "dbPepNeo",
        "antigen_type": "neoantigen", 
        "cancer": parts[0].strip(),
        "gene": parts[1].strip(),
        "HLA": parts[2].strip(),
        "peptide": parts[3].strip(),
        "net4_aff": parts[4].strip(),
        "net4_rank": parts[5].strip(),
        "binding_level": parts[6].strip(),
        "wt_peptide": parts[9].strip()
    }
    return entry
    
def create_itsndb_entry(parts):         # Done
    entry = {
        "source": "ITSNdb",
        "antigen_type": "neoantigen", 
        "cancer": parts[0].strip(),
        "wt_peptide": parts[4].strip(),
        "HLA": parts[5].strip(),
        "gene": parts[9].strip()
    }
    return entry
    
def create_itsndb_valids_entry(parts):         # Done
    entry = {
        "source": "ITSNdb",
        "antigen_type": "neoantigen", 
        "HLA": parts[2].strip(),
        "variant_type": parts[4].strip()
    }
    return entry

def create_neodb_entry(parts):          # Done
    entry = {
        "source": "NeoDB",
        "antigen_type": "neoantigen",
        "peptide": parts[1].strip(),
        "aa_change": parts[2].strip(),
        "HLA": parts[3].strip(),
        "wt_peptide": parts[4].strip(),
        "gene": parts[5].strip(),
        "cancer": parts[9].strip(),
        "variant_type": parts[10].strip()
    }
    return entry

def create_caatlas_entry(parts, cancer_name):        # Done
    entry = {
        "source": "CaAtlas",
        "antigen_type": "TAA",
        "peptide": parts[0].strip(),
        "aa_change": parts[3].strip(),
        "gene": parts[4].strip(),
        "cancer": cancer_name
     }
    return entry  

def create_caatlas_extra_entry(parts):        # Done
    entry = {
        "source": "CaAtlas",
        "antigen_type": "TAA",
        "peptide": parts[1].strip(),
        "HLA": parts[9].strip(),
        "cancer": parts[10].strip()
     }
    return entry 
    
def create_cadbio_entry(parts):         # Done
    entry = {
        "source": "CADBIO",
        "antigen_type": "TAA",
        "peptide": parts[3].strip(),
        "HLA": parts[6].strip(),
        "gene": parts[7].strip()
    }
    return entry  

def create_caped_entry(parts):          # Done
    entry = {
        "source": "CAPED",
        "antigen_type": "TAA", 
        "gene": parts[0].strip(),
#       "variant_tissue": parts[1],    # Some of these have the type of cancer
        "HLA": parts[2].strip(),       
        "peptide": parts[4].strip()
    }
    return entry  
    
def create_cedar_entry(parts):          # Done
    entry = {
        "source": "CEDAR",
        "antigen_type": "TAA", 
        "peptide": parts[0].strip(),
        "cancer": parts[1].strip(),
        "HLA": parts[2].strip()
    }
    return entry   
    
def create_cedar_extra_entry(parts):          # Done
    entry = {
        "source": "CEDAR",
        "antigen_type": "TAA", 
        "peptide": parts[2].strip(),
        "cancer": parts[3].strip(),
        "gene": parts[6].strip(),
        "variant_tissue": parts[7].strip()
    }
    return entry   

def create_ise_entry(parts, cancer_name, tissue_type):          # Done
    entry = {
        "source": "ISE",
        "antigen_type": "TAA", 
        "peptide": parts[0],
        "cancer": cancer_name,
        "variant_tissue": tissue_type
    }
    return entry   

def create_ise_normal_entry(parts):          # Done
    entry = {
        "source": "ISE",
        "antigen_type": "TAA", 
        "cancer": "Normal",
        "peptide": parts[1],
        "gene": parts[2],
        "variant_tissue": parts[4]
    }
    return entry
    
def create_ligand_entry(parts):         # Remove row numbers in 1st column
    entry = {
        "source": "LigandMHCatlas",
        "antigen_type": "TAA", 
        "HLA": parts[2].strip(),
        "peptide": parts[4].strip(),
        "binding_level": parts[5].strip(),
        "net4_aff": parts[6].strip(),
        "cancer": parts[12].strip()
    }
    return entry   
    
def create_tumoragdb_entry(parts):      # What to do with aa_wt
    entry = {
        "source": "TumorAgDB",
        "antigen_type": "neoantigen", 
        "variant_tissue": parts[2].strip(),
        "cancer": parts[3].strip(),
        "gene": parts[8].strip(),
        "aa_change": parts[11].strip() + parts[9].strip() + parts[10].strip(),
        "peptide": parts[12].strip(),
        "wt_peptide": parts[13].strip()     
    }
    return entry 
    
def create_tantigen_entry(parts):      # Done 
    entry = {
        "source": "TANTIGEN",
        "antigen_type": "TAA", 
        "peptide": parts[0].strip(),
        "HLA": parts[1].strip()
    }
    return entry     
    
def create_bedran_entry(parts):         # Done    
    entry = {
        "source": "Paper - Bedran et al.",
        "antigen_type": "TAA",
        "cancer": "Liver cancer",
        "peptide": parts[1].strip(),
        "gene": parts[2].strip()
     }
    return entry   
    
def create_beijer_entry(parts):         # Done
    entry = {
        "source": "Paper - Beijer et al.",
        "antigen_type": "TAA", 
        "cancer": "Liver cancer",
        "gene": parts[0].strip(),
        "peptide": parts[2].strip(),
        "variant_tissue": parts[7].strip()
    }
    return entry   
    
def create_loffier_entry(parts):        # Done
    entry = {
        "source": "Paper - Loffier et al.",
        "antigen_type": "TAA", 
        "cancer": "Liver cancer",
        "gene": parts[0].strip(),
        "HLA": parts[1].strip(),
        "peptide": parts[3].strip()
     }
    return entry   
    
def create_mizukoshi_entry(parts):      # Done
    entry = {
        "source": "Paper - Mizukoshi et .",
        "antigen_type": "TAA", 
        "cancer": "Liver cancer",
        "gene": parts[2].strip(),
        "peptide": parts[4].strip()
     }
    return entry   
    
def create_reparaz_entry(parts):        # Done
    entry = {
        "source": "Paper - Reparaz et al.",
        "antigen_type": "TAA", 
        "cancer": "Liver cancer",
        "HLA": parts[1].strip(),
        "peptide": parts[3].strip(),
        "gene": parts[4].strip(),
        "binding_level": parts[5].strip()
   }
    return entry   


def load_taxids(hmrgd_file):
    """
    Read HMRGD file once into memory as:
        taxid -> body site
    """
#    lookup = {}

    with open(hmrgd_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter=",")
        header = next(reader)
        taxid_idx = [i for i, col in enumerate(header) if "taxid" in col.lower()]
        lookup_list = [{} for _ in range(len(taxid_idx))]
#       next(reader, None)  # Skip header

        for row in reader:
            if len(row) < 3:
                continue

            organism = row[0].strip()
            body_site = row[1].strip()
 #          taxid = row[2].strip()
            for tid in taxid_idx:
                taxid = row[tid].strip()
#               lookup[taxid] = {
                lookup_list[tid - 2][taxid] = {
                    "organism": organism,
                    "body_site": body_site,
                }

    return lookup_list

def extract_full_matches(inputfile, pep_lookup):
    """
    Parse BLAST output once, enrich each row from in-memory peptide lookup,
    and collect only full matches.
    """

    full_matches = []
    expected_extra_fields = (
        "antigen_type", "variant_type", "variant_tissue", "cancer", "gene",
        "wt_peptide", "aa_change", "mutation", "HLA",
        "MHCf_rank", "net4_aff", "net4_rank",
        "location", "source", "title"
    )

    with open(inputfile, "r", newline="") as file:
        for line in file:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < len(BLAST_COLUMNS):
                print("Not enough columns")
                continue

            # BLAST fields
            base_row = {
                BLAST_COLUMNS[i]: parts[i]
                for i in range(len(BLAST_COLUMNS))
            }
            base_row["organism_in_human"] = ""
            base_row["microbe_location"] = ""

            # Find peptide matches
            peptide = base_row["sseq"].strip().upper()
            matches = pep_lookup.get(peptide, [])

            # Normalise single dict -> list of dicts
            if isinstance(matches, dict):
                matches = [matches]

            # No match
            if not matches:
                matches = [{}]

            # Process each match
            for match in matches:
                if match:
                    row = deepcopy(base_row)
                    row.update(match)

                    # Ensure expected fields exist
                    for key in expected_extra_fields:
                        if key not in row:
                            row[key] = ""

                    if is_full_match(row):
                        full_matches.append(row)

    return full_matches


def is_full_match(row):
    return (
        float(row["pident"]) == 100.0
        and int(row["mismatch"]) == 0
        and int(row["nident"]) == int(row["length"])
        and int(row["qcovs"]) == 100
    )


def is_one_mismatch(row):
    return (
        int(row["mismatch"]) == 1
        and int(row["nident"]) == int(row["length"]) - 1
    )


def is_two_mismatch(row):
    return (
        int(row["mismatch"]) == 2
        and int(row["nident"]) == int(row["length"]) - 2
    )


def build_match_tables(full_matches, taxid_lookups):
    """
    Build final match table rows entirely in memory.
    """
    match_tables = [[] for _ in range(len(taxid_lookups))]
    seen = set()
    taxonomy_cache = {}

    for match in full_matches:
        sci_names = [x.strip() for x in match["sscinames"].split(";")]
        taxids = [x.strip() for x in match["staxids"].split(";")]

        for idx, taxid in enumerate(taxids):
            if not taxid:
                continue

            for lookup_idx in range(len(taxid_lookups)):
                hmrgd_entry = taxid_lookups[lookup_idx].get(taxid)
                if not hmrgd_entry:
                    continue

                match["organism_in_human"] = hmrgd_entry["organism"]
                match["microbe_location"] = hmrgd_entry["body_site"]
                
                if taxid not in taxonomy_cache:
                    taxonomy_cache[taxid] = get_taxonomy(taxid) if taxid else empty_taxonomy()

                lineage = taxonomy_cache[taxid]
                table_row = create_table_dict(match, lineage, taxid)

                row_key = tuple((col, table_row.get(col, "")) for col in MATCH_TABLE_COLUMNS)
                if row_key in seen:
                    continue

                seen.add(row_key)
                match_tables[lookup_idx].append(table_row)

    return match_tables


def create_table_dict(row, lineage_dict, taxid):
    return {
        "Peptide Sequence": row.get("sseq", ""),
        "Peptide length": row.get("qlen", ""),
        "% Identity": row.get("pident", ""),
        "Query coverage": row.get("qcovs", ""),
        "Subject": row.get("saccver", ""),
        "Peptide (BLAST results)": row.get("qaccver", ""),
        "MHC Class": get_mhc_class(row.get("HLA", "")),
        "HLA": get_hla(row.get("HLA", "")),
        "Organism name": row.get("organism_in_human", ""),
        "Taxid": taxid,
        "Kingdom": lineage_dict.get("kingdom", ""),
        "Phylum": lineage_dict.get("phylum", ""),
        "Class": lineage_dict.get("class", ""),
        "Order": lineage_dict.get("order", ""),
        "Family": lineage_dict.get("family", ""),
        "Genus": lineage_dict.get("genus", ""),
        "Species": lineage_dict.get("species", ""),
        "Mutation": lineage_dict.get("mutation", ""),
        "MHCf_rank": row.get("MHCf_rank", ""),
        "net4_aff": row.get("net4_aff", ""),
        "net4_rank": row.get("net4_rank", ""),
        "Antigen type": row.get("antigen_type", ""),
        "Variant type": row.get("variant_type", ""),
        "Cancer": row.get("cancer", ""),
        "Sample Tissue": row.get("variant_tissue", ""),
        "Gene": row.get("gene", ""),
        "WT Peptide": row.get("wt_peptide", ""),   
        "AA Change": row.get("aa_change", ""),   
        "Source": row.get("source", ""),
        "Microbe location": row.get("microbe_location", ""),
        "Title": row.get("stitle", "")
    }


def empty_taxonomy():
    return {
        "kingdom": "",
        "phylum": "",
        "class": "",
        "order": "",
        "family": "",
        "genus": "",
        "species": "",
    }


def get_taxonomy(taxid):
    """
    Call taxonkit once per unique taxid, with caching handled by caller.
    """
    result = subprocess.run(
        ["/srv/scratch/z8039617/neoantigens/bin/taxonkit", "lineage", "-R"],
        input=f"{taxid}\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=False,
    )

    if result.returncode != 0 or not result.stdout.strip():
        return empty_taxonomy()

    parts = result.stdout.rstrip("\n").split("\t")
    if len(parts) < 3:
        return empty_taxonomy()

    lineage_names = parts[1].split(";")
    lineage_ranks = parts[2].split(";")

    taxonomy = empty_taxonomy()
    for i, rank in enumerate(lineage_ranks):
        name = lineage_names[i].strip() if i < len(lineage_names) else ""
        rank = rank.strip()
 
        if rank == "kingdom":
            taxonomy["kingdom"] = name
        elif rank == "phylum":
            taxonomy["phylum"] = name
        elif rank == "class":
            taxonomy["class"] = name
        elif rank == "order":
            taxonomy["order"] = name
        elif rank == "family":
            taxonomy["family"] = name
        elif rank == "genus":
            taxonomy["genus"] = name
        elif rank == "species":
            taxonomy["species"] = name

    return taxonomy


def get_hla(hla):
    hla_str = hla
    if len(hla) > 0 and not hla_str.startswith("HLA-"):
        hla_str = "HLA-" + hla_str
    return hla_str
 
 
def get_mhc_class(hla):
#    if not hla or "-" not in hla:
    if not hla:
        return ""

 #   parts = hla.split("-")
 #   if len(parts) < 2:
 #       return ""

#    locus = parts[1]
    locus = hla
    if any(x in locus for x in ("DP", "DM", "DO", "DQ", "DR")):
        return "II"
    if any(x in locus for x in ("A", "B", "C")):
        return "I"
    return ""


def write_dict_csv(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
	