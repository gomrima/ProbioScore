# AHP weight provenance (ProbioScore 1.0)

## Source

The AHP weights embedded in this release come from the sheet
`CURATED_AGG_WEIGHTS_CR0.10` of `AHP_ProScore_CALCULATOR_Results.xlsx`,
produced by the AHP calculator generator version 8 from a panel of 51
experts. For every section, only the expert pairwise matrices with a
consistency ratio at or below 0.10 were retained, then aggregated by the
geometric mean of judgments, and the section weights were taken as the
exact principal eigenvector of the aggregated matrix. Every retained
section matrix has a consistency ratio well below 0.10.

Retained matrices per section: pillars 15, safety 21, survival 17,
benefits 23 (out of 51 experts).

## Weights used by the engine

The positive utility axis uses the survival sub-pillar weights (B1 to B5)
and the benefit sub-pillar weights (C1 to C5) directly, and combines the
two positive pillars with a renormalized split.

Survival sub-pillars (sum to 1):

| code | criterion | weight |
|---|---|---|
| B1 | Acid and bile tolerance | 0.396756 |
| B2 | General stress response | 0.149943 |
| B3 | Mucus adhesion capacity | 0.292990 |
| B4 | Exopolysaccharide production | 0.089570 |
| B5 | Defense systems | 0.070741 |

Benefit sub-pillars (sum to 1):

| code | criterion | weight |
|---|---|---|
| C1 | Antimicrobial activity | 0.239838 |
| C2 | Beneficial metabolites | 0.284099 |
| C3 | Carbohydrate metabolism | 0.087796 |
| C4 | Antioxidant enzymes | 0.107218 |
| C5 | Immunomodulatory potential | 0.281048 |

Positive pillar split. The expert panel ranked the three pillars as
Safety 0.690418, Survival 0.179701, Benefits 0.129881. ProbioScore treats
safety as a non-compensatory veto handled by the fuzzy risk and pathogen
layer, so safety is not part of the compensatory positive utility. The two
positive pillars are therefore renormalized to sum to 1 over themselves:

| pillar | panel weight | renormalized |
|---|---|---|
| Survival and Colonization | 0.179701 | 0.580462 |
| Functional and Metabolic Benefits | 0.129881 | 0.419538 |

## Safety sub-pillar weights

The safety sub-pillar weights (A1 to A4) are recorded in the configuration
for completeness and traceability, but they are not consumed by the engine:
the safety axis is evaluated as a non-compensatory tiered burden, not as an
AHP weighted sum.

| code | criterion | weight |
|---|---|---|
| A1 | Antimicrobial resistance | 0.192455 |
| A2 | Virulence | 0.305226 |
| A3 | Toxin production | 0.395144 |
| A4 | Genomic instability | 0.107176 |
