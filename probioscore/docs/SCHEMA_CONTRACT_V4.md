# V4 schema contract

The raw V4 scoring path expects a directory with V4 TSV subfolders used by the embedded base engine. Required components include final matrices, ProbioSML category summaries, safety matrices, toxin matrices, mobile safety summaries, BGC outputs, CRISPR summaries and toxin antitoxin summaries.

For audit replay, the `--input-tsv` path accepts a prebuilt table with `genome_id` and `ProbioScore_Status`. Additional V4 derived columns are preserved.

Any external panel (e.g. the 634-genome validation panel) is scored with --mode prospective_frozen.
