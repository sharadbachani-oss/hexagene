# Getting access to HRS genetic / genotype data

A practical, click-by-click guide for obtaining Health and Retirement Study (HRS)
genetic data for use in the HexaGene research program (variant-pathogenicity
validation against a real cohort's genotypes + health phenotypes).

> **Important — what an assistant can and cannot do here.** The submission steps
> below require *your* personal credentials and, for raw genotypes, a
> *legally-binding agreement signed by your institution's authorized official*.
> Those steps must be performed by you, in your own browser, under your own
> identity — they cannot be automated or submitted on your behalf. This guide
> does the heavy lifting around them (the roadmap, the drafted research-use
> statement, the drafted order request) so that the parts only you can do are
> fast and correct.

---

## TL;DR — there are two access tiers, and you need both

| What you want | Tier | Mechanism | Difficulty | Prereqs |
|---|---|---|---|---|
| **Polygenic Scores (PGS)** | Sensitive Health Data | Order Form on the HRS data portal | **Easy** | The account you already have |
| **APOE / serotonin transporter / 2008 telomere** | Sensitive Health Data | Order Form on the HRS data portal | **Easy** | The account you already have |
| **Raw SNP genotypes V1 (2006–08) / V2 (2006–10)** | dbGaP controlled access (`phs000428`) | dbGaP Authorized Access application | **Hard** | eRA Commons account + institutional Signing Official + Data Access Committee approval |
| **Raw SNP genotypes V3 (2006–12)** | NIAGADS controlled access (`NG00119`) | NIAGADS Data Access Request (DAR) | **Hard** | Institutional verification / DAR approval |

You (Merlin Digital / Hexagene) **already hold the Sensitive Health Data tier** —
the April 2026 grant for Biomarker Data is exactly this tier. So PGS / APOE /
telomere are an *order-form addition*, not a new application. Raw genotypes are a
separate, much heavier process.

**Recommended order of operations:**
1. **Do the easy tier first** (PGS + APOE + telomere) — see
   [`draft-sensitive-data-order-email.md`](./draft-sensitive-data-order-email.md).
   Often enough on its own for a genetics product.
2. **Start the dbGaP track in parallel** if you genuinely need variant-level
   genotypes — it is the long pole (institutional registration + DAC review).

---

## Your current standing (from your records)

- **HRS data account holder:** Suhail Bachani
- **Login:** sharad.bachani@merlin-me.com (portal: https://hrsdata.isr.umich.edu)
- **Affiliation on file:** Merlin Digital, Dubai, UAE
- **Granted:** HRS Sensitive Health Data → Biomarker Data (April 2026)

> Note: the account name on file ("Suhail Bachani") differs from the name you sign
> with ("Sharad Bachani, CEO"). For the **easy tier** this is fine. For the
> **dbGaP tier** the Principal Investigator named on the application must be a real,
> identifiable researcher whose name matches their eRA Commons account — so decide
> up front *who the named PI will be* and keep it consistent everywhere.

---

## Tier 1 — Sensitive Health Data (PGS, APOE, telomere) — EASY

You already have an account with this tier, so the genetic Sensitive Health Data
products are an incremental request, not a fresh registration.

**Steps:**
1. Log in at **https://hrsdata.isr.umich.edu** (reset password at
   https://hrsdata.isr.umich.edu/user/password if needed).
2. Open the **Special Access Downloads** page
   (https://hrsdata.isr.umich.edu/data-products/special-access-downloads) and check
   whether **PGS**, **APOE & serotonin transporter**, and **2008 telomere** are
   already downloadable under your account. (Sometimes a grant for one Sensitive
   Health Data product covers others; sometimes each must be requested.)
3. For any that are **not** already available, submit the **Sensitive Health Data
   Order Form**: https://hrs.isr.umich.edu/data-products/sensitive-health-data-order-form
   — listing PGS, APOE/serotonin transporter, and 2008 telomere.
   A ready-to-send version of this request is in
   [`draft-sensitive-data-order-email.md`](./draft-sensitive-data-order-email.md).
4. Questions / follow-ups: `hrsdatareq@umich.edu` (requests) and
   `hrsquestions@umich.edu` (general).

**Product pages (what each contains):**
- PGS: https://hrsdata.isr.umich.edu/data-products/polygenic-score-data-pgs
- APOE & serotonin transporter: https://hrsdata.isr.umich.edu/data-products/apoe-and-serotonin-transporter-alleles
- 2008 telomere: https://hrsdata.isr.umich.edu/data-products/2008-telomere-data

**Typical turnaround:** days to ~2 weeks.

---

## Tier 2 — Raw SNP genotypes via dbGaP (`phs000428`) — HARD

This is real controlled-access human genetics. There is no shortcut: NIH gates it
behind verified institutional identity and a legal data-use certification.

### Prerequisites (the long pole)
1. **An eRA Commons account for the named PI.** Before that can exist, the
   **institution** (Merlin Digital / Hexagene) must be **registered in eRA Commons**,
   and must designate at least one **Signing Official (SO)** and an account
   administrator. Foreign and commercial organizations *can* register — it is just
   a multi-step, multi-week process.
   - eRA Commons institutional registration:
     https://www.era.nih.gov/register-accounts/register-your-organization-in-era-commons.htm
2. Confirm the study's **Data Use Limitations (DUL)** allow your use. Open the HRS
   dbGaP study page and read the consent groups / DUL **before** applying — some
   consent groups restrict to non-profit or US use, and as a for-profit, non-US
   entity you must confirm your intended use is permitted:
   - dbGaP study: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000428

### Application steps (performed by you, as PI, in your browser)
1. Sign in to the **dbGaP Authorized Access** system with the PI's eRA Commons
   credentials: https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=login
2. **Create a Data Access Request project** and select study **`phs000428`** and the
   consent group(s) you need.
3. Paste your **Research Use Statement** (drafted for you in
   [`research-use-statement.md`](./research-use-statement.md)) and a non-technical
   summary.
4. **Submit** → your institution's **Signing Official electronically co-signs the
   Data Use Certification (DUC)** → the **HRS Data Access Committee** reviews.
5. On approval, download the genotype + imputed dosage files from dbGaP (see HRS's
   walkthrough video: https://hrs.isr.umich.edu/documentation/video-tutorials/9110).

### Then: link genotypes back to HRS survey data
The dbGaP files use de-identified IDs. To merge with HRS survey/biomarker variables
you must separately obtain the **HRS↔dbGaP Cross-Reference file**:
- Submit the **Genetic Data Cross-Reference Request Form**, **and**
- Submit a signed **Genetic Data Access Use Agreement**.
- Both are linked from: https://hrs.isr.umich.edu/data-products/genetic-data/access-summary

**Typical turnaround:** weeks-to-months — dominated by eRA Commons institutional
registration (if not already done) and DAC review.

### Alternative for V3 (2006–2012): NIAGADS
Genotype Data Version 3 is distributed by NIAGADS, not dbGaP, via a **Data Access
Request (DAR)**:
- Dataset NG00119: https://dss.niagads.org/datasets/ng00119/

---

## Honest assessment for a Dubai-based commercial entity

- **Tier 1 (PGS/APOE/telomere):** straightforward — you already qualify. **Start here.**
- **Tier 2 (raw genotypes):** achievable but non-trivial. The two real obstacles are
  (a) standing up an **eRA Commons institutional registration + Signing Official**
  for a for-profit, non-US company, and (b) satisfying the study's **Data Use
  Limitations** as a commercial/foreign requester. Budget months, and confirm the
  DUL permits for-profit use *before* investing in the registration.
- If PGS turns out to cover the HexaGene validation need, you may not need Tier 2 at
  all.

---

## Sources
- HRS Genetic Data Products — https://hrs.isr.umich.edu/data-products/genetic-data/products
- HRS Genetic Data Access Summary — https://hrs.isr.umich.edu/data-products/genetic-data/access-summary
- HRS Sensitive Health Data Order Form — https://hrs.isr.umich.edu/data-products/sensitive-health-data-order-form
- HRS Genotype Data V2 — https://hrs.isr.umich.edu/data-products/genetic-data/genotype-data-v2
- dbGaP study phs000428 — https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000428
- NIAGADS NG00119 (V3) — https://dss.niagads.org/datasets/ng00119/
- Downloading HRS genetic data from dbGaP (video) — https://hrs.isr.umich.edu/documentation/video-tutorials/9110
- eRA Commons — register your organization — https://www.era.nih.gov/register-accounts/register-your-organization-in-era-commons.htm
