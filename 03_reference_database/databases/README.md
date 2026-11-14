# herbivore_honours_metabarcoding_check

Checking and updating the reference databases used for herbivore diet metabarcoding (trnL and rbcL plant markers, Kruger vegetation).

## Why this dir exists

The trnL and rbcL fasta databases got updated with new accessions (KVW versions), but the master CSV never got updated to match. This dir holds the check that found the gap, and the notebook that fixed it.

## Key files

- `trnL_DB_final.fasta` / `rbcL_DB_final.fasta` - original (pre-update) fasta databases
- `trnL_final_updated_KVW_DB.fasta` / `rbcL_final_updated_KVW_DB.fasta` - updated fasta databases
- `OG_Complete_reference_database_kyle.csv` - the old master CSV (Accession, Barcode, Family, Genus, Species). Confirmed to match the original fastas exactly, but missing all the new accessions.
- `Kruger_Veg_Species_list_kyle.csv` - target species list (Phylum, Class, Order, Family, Genus, Species) used for the coverage check.
- `Complete_reference_database_kyle_FINAL.csv` - the fixed master CSV, produced by the notebook (OG CSV + all newly added accessions).
- `database_checks.ipynb` - the R notebook with all the checks. Run top to bottom.
- `outputs/` - CSV results from each check block.

## What the notebook does (in order)

1. **Compare updated vs original fastas** (trnL and rbcL). Sequence/species/genus counts, full list of new accessions added, any removed.
2. **Check the OG CSV against the updated fastas.** Confirmed all 58 new accessions (14 trnL, 44 rbcL) were missing from the CSV.
3. **Build the final master database.** Appends the 58 new accessions onto the OG CSV and saves `Complete_reference_database_kyle_FINAL.csv`. Note: the new rows have a blank `Family` column, fasta headers don't carry family-level taxonomy, so that needs filling in manually.
4. **Coverage check.** Compares `Complete_reference_database_kyle_FINAL.csv` against `Kruger_Veg_Species_list_kyle.csv` to see how much of the target species list actually has a reference sequence, per marker and combined.
5. **Final sanity check.** Re-runs the accession check from step 2, but against the final database instead of the OG CSV. Should come back with zero mismatches, confirming the fix worked.

## Results summary (for reviewer comment, L200: "How many new accessions were added? Would be interesting to know how (in)complete the libraries previously were.")

**58 new accessions added, 0 removed:**

- trnL: 920 → 934 sequences (+14). 449 → 463 species (+14 new). 270 → 275 genera (+5 new).
- rbcL: 1238 → 1282 sequences (+44). 557 → 584 species (+27 new). 323 → 335 genera (+12 new).

The old master CSV (`OG_Complete_reference_database_kyle.csv`) reflected only the pre-update fastas. All 58 new accessions were missing from it, confirmed in Block 2. Fixed in `Complete_reference_database_kyle_FINAL.csv` (2158 → 2216 rows).

**Library completeness vs the Kruger target species list:** TODO, fill in from Block 4's `coverage_table` once run (`outputs/veg_species_coverage_summary.csv`). Gives species-level coverage, genus-level-only coverage, and no-match count/percent, per marker and combined.

## Other notes

- Found and fixed a header-parsing bug along the way: RefSeq-style accessions (`NC_065250.1`) have an internal underscore that broke the original split-on-underscore parser. Fixed by parsing from both ends of the header instead.
- The old `Reference_database/` working folder (Python scripts, old comparisons) was moved to Trash before this notebook was built. `OG_Complete_reference_database_kyle.csv` at the root is confirmed byte-identical to the copy that used to live in there.
