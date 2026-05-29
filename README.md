# ortholog-pipeline

Loads gene orthologs into RGD.

Two flows live in this pipeline:

1. **Alliance (AGR) flow** (`--agrOrthologs`) &mdash; downloads the Alliance combined ortholog TSV and populates the `AGR_ORTHOLOGS` table.
2. **HCOP / NCBI flow** (`--species ...`) &mdash; populates `GENETOGENE_RGD_ID_RLT` (orthologs) and `RGD_ASSOCIATIONS` (weak orthologs) per species using a four-tier priority cascade: **manual &rarr; Alliance &rarr; HCOP &rarr; NCBI**. Alliance candidates come from `AGR_ORTHOLOGS` (mutual-best rows), so the AGR flow must be run first; `runAll.sh` enforces this.

## Entry point

Main class: `edu.mcw.rgd.dataload.OrthologRelationLoadingManager` (wired through Spring from `properties/AppConfigure.xml`).

Command-line flags:

| Flag | Effect |
|---|---|
| `--species <name>` | Process a single species (`rat`, `mouse`, `dog`, `pig`, &hellip;). `--species human` is rejected. |
| `--species all` | Iterate over all searchable species (excluding `HUMAN`) and load each one. |
| `--agrOrthologs` | Run the Alliance TSV loader (`AgrTsvLoader`) instead of the HCOP/NCBI flow. |
| `--fixXRefDataSet` | Maintenance: dedupe `XREF_DATA_SET` values in orthologs and `weak_ortholog` associations (fixes duplicate fields emitted by HCOP). |

## HCOP / NCBI flow

The HCOP/NCBI flow always runs **per species**, processing only the human&harr;species pair on a given invocation. `--species all` loops over every searchable RGD species (excluding `HUMAN`). For each run the manager logs old / new ortholog and `weak_ortholog` association counts, the difference, and elapsed time.

Before any per-species work, the manager calls `checkAllianceFreshness()`: the run aborts immediately if `AGR_ORTHOLOGS` is empty or its newest `last_update_date` is older than `agrMaxAgeDays` (default 60, configurable on the `orthologRelationLoadingManager` bean). Run `--agrOrthologs` (or `loadAgrOrthologs.sh`) first to refresh.

### 1. Download and parse (`OrthologRelationFile`, `OrthologRelationParser`)

- If the species is in the HCOP set (Human, Rat, Mouse, Dog, Pig), the HCOP file is downloaded **first** (`https://storage.googleapis.com/.../human_all_hcop_sixteen_column.txt.gz` &rarr; `data/hgnc.data`), parsed and filtered down to relations between human and the target species, then the NCBI file is downloaded (`https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz` &rarr; `data/ncbi.data`) and its relations are appended. For non-HCOP species only the NCBI file is read.
- Files are kept gzip-compressed; previous downloads are kept under date-stamped backup names.
- A sanity check throws if fewer than 5,000 human/species relations come out of the parse.
- Each surviving line becomes an `OrthologRelation` (src/dest EntrezGene IDs, species keys, source name &mdash; `HGNC` or `NCBI` &mdash; and a dataset/evidence string).

### 2. Resolve EntrezGene IDs to RGD IDs (`OrthologRelationLoader.matchRgdId`)

Each side of each relation is resolved (in parallel, with a memoization map):

- Single ACTIVE RGD gene &rarr; use that RGD ID.
- Multiple ACTIVE genes &rarr; counted as `multipleMatchC`, logged to `multipleMatch.log`, relation dropped.
- No ACTIVE gene but one WITHDRAWN gene with a single active replacement &rarr; use the replacement.
- WITHDRAWN with no active replacement &rarr; counted as `withdrawnC`, logged to `withdrawn.log`, relation dropped.
- No match at all &rarr; counted as `unMatchedC`, logged to `unmatched.log`, relation dropped.

Successful matches are written to `matched.log`. Relations that did not resolve on both sides are dropped at the start of `loadData` (`dropUnmappedRelations`).

### 3. Group by human RGD ID (`buildGroups`, `OrthologGroup.add`)

All relations are grouped on the **human RGD ID**. Inside one group:

- Duplicate (src,dest) pairs from the same data source are merged; the dataset/evidence string is union-ed.
- A duplicate (src,dest) pair where one source is HGNC and the other is NCBI is merged by appending `NCBI` to the HGNC dataset.
- For each group, reverse (species&rarr;human) relations are synthesized via `buildComplementaryRelations` so both directions of the ortholog are loaded.

### 4. Pick the strong ortholog &mdash; four-tier priority cascade (`qcGroups`, `generateOrtholog`)

Within a group, relations are partitioned by species pair (e.g. `human|rat`, `rat|human`). For each species pair the loader picks **one** strong ortholog and lines up everything else as weak orthologs.

The cascade for each `(srcRgdId, destSpeciesTypeKey)` is **manual &rarr; Alliance &rarr; HCOP &rarr; NCBI**; the first non-empty tier wins and lower tiers are ignored for selection (their relations may still survive as `weak_ortholog` associations &mdash; see step 6). Priority numbers used elsewhere in the code: `RGD=4, Alliance=3, HGNC=2, NCBI=1` (see `OrthologRelationDao.sourcePriority`).

1. **Manual** &mdash; `getManualOrthologs(srcRgdId, destSpeciesTypeKey)` returns orthologs with `XREF_DATA_SRC='RGD'`. If exactly one exists it is used as the strong ortholog. If multiple exist, the conflict is logged and no ortholog is picked for that pair this round.
2. **Alliance** &mdash; `getAllianceOrthologs(srcRgdId, destSpeciesTypeKey)` queries `AGR_ORTHOLOGS` for mutual-best rows (`is_best_score='Y' AND is_best_rev_score='Y'`) in either direction. Synthetic `Ortholog` objects are produced with `XREF_DATA_SRC='Alliance'` and `XREF_DATA_SET=methods_matched`. Multiple matches log a conflict and no ortholog is picked.
3. **HCOP** &mdash; best-fit over the incoming relations filtered by `dataSource='HGNC'`.
4. **NCBI** &mdash; best-fit over the incoming relations filtered by `dataSource='NCBI'`.

The best-fit picker (tiers 3 and 4) breaks ties in this order:
  1. Only one relation &rarr; pick it.
  2. Largest number of evidence entries (commas in the dataset string) &rarr; pick the winner.
  3. Destination gene symbol matches the source gene symbol (case-insensitive) &rarr; pick that one.
  4. Otherwise sort destinations alphabetically and pick the first.

Tier choices are counted under `ORTHOLOG SOURCE STATS` and tie-break methods under `ORTHOLOG BEST FIT STATS` in `process.log`.

### 5. Insert / match / downgrade (`loadGroups`, `getKeyForMatchingOrtholog`)

For each chosen strong ortholog, the DAO inspects what is already in `GENETOGENE_RGD_ID_RLT` for the same (srcRgdId, destSpeciesTypeKey):

- **Match**: an existing row already points to the same destination gene &rarr; bump `last_modified_date`.
- **Insert**: no existing row for that pair &rarr; insert the new ortholog.
- **Conflict** (existing row points to a different destination):
  - `compareOrthologs` ranks by **source priority first** (RGD > Alliance > HGNC > NCBI), then by evidence count, matching symbol, alphabetic. The higher-priority candidate wins outright regardless of evidence count.
  - If the incoming candidate wins, the existing row is deleted and the incoming one is inserted.
  - If the existing one wins, the incoming relation is downgraded to a `weak_ortholog` association.
- **Multiple existing rows** for the same pair (data drift): keep the strongest, queue the rest for deletion.

All deletes go through `deleteStaleOrtholog`, which refuses to remove a manual (`RGD`) ortholog and refuses to remove the sole remaining row for a (srcRgdId, destSpeciesTypeKey) pair.

### 6. Weak orthologs in `RGD_ASSOCIATIONS` (`processOrthologAssociations`)

The non-winning relations within each group, plus any orthologs downgraded in step 5, are collected as `weak_ortholog` associations. The loader then:

- Drops any association whose two genes already have a strong ortholog between them (`areGenesOrthologous`).
- Compares the incoming set to existing `weak_ortholog` rows for the species pair via `OrthoAssociationSyncer`, then inserts / updates / deletes accordingly.
- Before issuing deletes, it checks complement rows (`checkComplementOrthologs`, `checkComplementAssociations`) so that an asymmetry that is being repaired in the same run is not flapped delete-then-insert.

### 7. Stale ortholog and duplicate cleanup

- `handleStaleOrthologs` looks at orthologs for the species pair whose `last_modified_date` was not bumped this run and tries to delete them, again skipping manual orthologs and never removing the last row for a pair. Alliance and HCOP/NCBI rows are re-touched by the cascade each run, so they don't appear stale unless their underlying source actually stopped emitting them.
- `deleteDuplicateNonManualOrthologs` runs at the end of `loadData`. For every `(src_rgd_id, dest_rgd_id)` pair with multiple rows, an Oracle `ROW_NUMBER() OVER (PARTITION BY src_rgd_id, dest_rgd_id ORDER BY CASE xref_data_src ... DESC)` keeps the highest-priority row and queues the rest for deletion &mdash; but only rows with `created_by = orthologPipelineId` are actually deleted, so manual rows and rows owned by other pipelines are never touched.

The net effect for curators: **manual orthologs in RGD are never overwritten, never deleted, and always preferred over any incoming HCOP / NCBI / Alliance assignment**. Below manual, Alliance mutual-best wins over HCOP, and HCOP wins over NCBI. Losers at any tier land in `RGD_ASSOCIATIONS` as weak orthologs instead.

## Source files

| File | Role |
|---|---|
| `OrthologRelationLoadingManager` | Main class; CLI parsing, run orchestration, metrics summary, RGD logging. |
| `OrthologRelationFile` | Downloads HCOP or NCBI source file per species (uses `FileDownloader2`, gzip + date-stamped backups). |
| `OrthologRelationParser` | Parses downloaded file into `OrthologRelation` records for the target species (with sanity-check minimum). |
| `OrthologRelationLoader` | Matches relations to RGD genes, builds groups, picks strong vs weak orthologs, writes to DB. |
| `OrthologRelationDao` | DAO for orthologs and associations; defines the "manual orthologs always win" rules. |
| `OrthologRelation` / `OrthologGroup` | Data classes &mdash; one relation and a group of relations keyed by human RGD ID. |
| `OrthoAssociationSyncer` | Extends `AssociationSyncer` (from rgdcore); disables indel optimization for ortholog associations. |
| `AgrTsvLoader` | Standalone loader for the Alliance TSV &rarr; `AGR_ORTHOLOGS` table. Refuses to run if obsolete-row deletion would exceed `obsoleteOrthologsDeleteThreshold` (default `10%`). |
| `HomologeneLoader` | Standalone tool (own `main`) that loads from NCBI HomoloGene. Not wired into the production flow. |
| `OrthoTool` | Standalone reporting tool (own `main`); emits a tab-delimited ortholog report for an input gene list. |

## Configuration &mdash; `properties/AppConfigure.xml`

Spring beans:

- `orthologRelationLoadingManager` &mdash; `xrefDataSrc=HGNC`, `agrMaxAgeDays=60` (per-species run aborts if `AGR_ORTHOLOGS` is empty or newest row is older than this many days)
- `orthologRelationLoader` &mdash; `defaultDataSetName=HCOP`, `createdBy=70`
- `orthologRelationDao` &mdash; `directOrthologTypeKey=11`, `transitiveOrthologTypeKey=13`
- `orthologRelationFile` &mdash;
  - HCOP source: `https://storage.googleapis.com/public-download-files/hcop/human_all_hcop_sixteen_column.txt.gz` &rarr; `data/hgnc.data`
  - NCBI source: `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz` &rarr; `data/ncbi.data`
  - HCOP species: Human, Rat, Mouse, Dog, Pig (everything else goes through NCBI)
- `agrTsvLoader` &mdash;
  - Source: `https://fms.alliancegenome.org/download/ORTHOLOGY-ALLIANCE_COMBINED.tsv.gz`
  - Processed species: rat, human, mouse, zebrafish, fly, roundworm, yeast, African clawed frog, tropical clawed frog
- `orthologRelationParser`, `rgdLogger`

The Oracle DataSource bean is loaded externally via `-Dspring.config=...` (see `run.sh`).

## Logging &mdash; `properties/log4j2.xml`

Per-flow loggers (each routes to its own monthly-rolled file under `logs/`):

- **HCOP/NCBI flow**: `status.log`, `detail.log`, `process.log`, `matched.log`, `unmatched.log`, `multipleMatch.log`, `withdrawn.log`, `assoc.log`, `rat_genes_without_ortholog.log`, `crosslinked_orthologs.log`, `inserted_orthologs.log`, `deleted_orthologs.log`
- **AGR flow**: `agr_status.log`, `agr_summary.log`, `agr_debug.log`, `inserted_agr_orthologs.log`, `inserted_agr_genes.log`, `deleted_agr_orthologs.log`

Status loggers also mirror INFO-level output to the console.

## Build

- Java 17, Gradle (`./gradlew build`)
- Dependencies: `commons-dbcp2`, `log4j-core`, `ojdbc11`, `spring-jdbc`, plus `rgdcore` (`lib/rgdcore_*.jar`)
- Distribution artifact: `build/distributions/ortholog-pipeline.zip`

## Run scripts &mdash; `src/main/dist/`

| Script | Purpose |
|---|---|
| `run.sh` | Base launcher; sets `-Dspring.config` and `-Dlog4j.configurationFile`, forwards args to the jar, tees to `run.log`. |
| `loadSpecies.sh <species>` | Load one species via `run.sh --species <species>`. |
| `runAll.sh` | Refreshes Alliance (`--agrOrthologs`) first, aborts if that fails, then runs for all searchable species (`loadSpecies.sh all`); emails `logs/status.log` on completion. |
| `loadAgrOrthologs.sh` | Run AGR flow; emails `logs/agr_summary.log` on completion. |
| `fixXRefDataSrc.sh` | Run maintenance `--fixXRefDataSet`; emails the result. |
