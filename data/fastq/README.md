# FASTQ input

Raw reads go here. Exact filenames pending.

Expected: one paired-end WTA library (BD Rhapsody Enhanced V2/V3 beads).
R1 = cell label + UMI, R2 = cDNA insert.

These were processed with the BD Rhapsody Sequence Analysis Pipeline v2.2.1
against `RhapRef_Mouse_WTA_2023-02` (Sample Tag version `mm`); that step is not
included in this repo, and its outputs are pre-staged in `data/bd_pipeline_out/`.

Not tracked in version control — these are the only files in this repo that
cannot be regenerated.
