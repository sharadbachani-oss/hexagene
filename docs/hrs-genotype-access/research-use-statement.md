# Draft Research Use Statement (RUS)

For the dbGaP Authorized Access application for **HRS study `phs000428`** (raw SNP
genotypes), and reusable as the long-form "intended use" for the Sensitive Health
Data Order Form.

> This is a starting draft built from the HexaGene repository description
> (sequence-intrinsic biophysical variant-pathogenicity scoring, validated against
> ClinVar/REVEL). **Edit the bracketed `[...]` placeholders** and trim to the word
> limit shown on the dbGaP form before submitting. Keep the named PI consistent
> with the eRA Commons account.

---

## Project title
Population-scale validation of a sequence-intrinsic biophysical variant-pathogenicity
score (HexaGene) against genotyped cohort data and longitudinal health outcomes.

## Principal Investigator
[Named PI — must match eRA Commons account], Merlin Digital, Dubai, UAE.

## Research Use Statement (long form)

We have developed HexaGene, a sequence-intrinsic biophysical scoring method that
predicts the pathogenicity of missense variants from local thermodynamic and
mechanical properties of the surrounding nucleotide context, independent of
evolutionary conservation. In prior work the score showed discrimination of
clinically classified variants (AUC ≈ 0.72 overall) and, importantly, orthogonality
to conservation-based predictors such as REVEL (r ≈ 0.27), giving it added value in
the classification "grey zone" where conventional predictors are uninformative.

The proposed work uses HRS genotype data to evaluate whether HexaGene scores carry
information beyond curated clinical labels when applied at population scale. Using
the genotyped HRS sample (`phs000428`), we will: (1) annotate each participant's
coding variants with HexaGene scores; (2) construct per-individual aggregate burden
measures from HexaGene scores; and (3) test the association of these measures with
HRS health phenotypes and longitudinal outcomes, benchmarking against established
references where available (e.g., HRS polygenic scores and clinical biomarkers). The
goal is methodological: to quantify the added predictive value, calibration, and
failure modes of a conservation-independent variant score in a real, deeply
phenotyped cohort.

We will not attempt to re-identify participants, will not contact participants, and
will not link these data to any other individually identifiable dataset beyond what
HRS authorizes (e.g., the HRS↔dbGaP cross-reference file under its own agreement).
Analyses are at the aggregate/statistical level and serve to improve a variant
interpretation method, not to make individual diagnoses.

## Non-technical summary
We are testing whether a new way of scoring DNA changes — based on the physics of
the DNA sequence rather than evolutionary conservation — predicts health outcomes in
a large, well-studied group of older adults. This helps judge whether the method
adds value for interpreting genetic variants of uncertain significance.

## Data security & cloud use
Data will be stored on [encrypted storage / institutional or approved cloud
environment], access-restricted to the named project personnel below, not
redistributed, and handled in accordance with the NIH Genomic Data Sharing Policy
and the study's Data Use Certification. [If using a cloud platform, name it and
confirm it meets dbGaP security expectations.]

## Personnel with access
- [PI name, role]
- [Analyst / collaborator names and roles]

## Data Use Limitations check
Before submitting, confirm on the dbGaP study page that the consent group(s) you
request permit your intended use. As a **for-profit, non-US** requester, explicitly
verify there is **no non-profit-only or US-only restriction** on the consent group
you select:
https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000428

---

### Reuse for the Sensitive Health Data Order Form
For PGS / APOE / telomere (the easy tier), you do not need the full RUS — a tightened
2–4 sentence version of the "Research Use Statement (long form)" above is enough.
That short version is already embedded in
[`draft-sensitive-data-order-email.md`](./draft-sensitive-data-order-email.md).
