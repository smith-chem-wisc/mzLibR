# Changelog

## mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb’s `since.mzlibr` in
the bridge’s verb specs, which every help page’s **Since** section
renders.

This release projects mzLib 1.0.595. The verbs marked *needs the pyMzLib
0.2.0 bridge*, *needs the pyMzLib 0.3.0 bridge*, *needs the pyMzLib
0.4.0 bridge* or *needs the pyMzLib 0.6.0 bridge* below are new in the
bridge that release publishes; with an older bridge they are refused,
naming the release they need.
[`mzlibr_install_bridge()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md)
installs the bridge pyMzLib 0.6.0 publishes (mzLib 1.0.595, `d10017ce`),
verified against that release’s SHA256SUMS.

### Statistics, designs and PTM-QTL (the pyMzLib 0.6.0 bridge)

Ported from pyMzLib 0.6.0, verb for verb; each *needs the pyMzLib 0.6.0
bridge*. Every number is mzLib’s, and every recording the help pages and
articles replay is pyMzLib’s, byte for byte. New articles:
[`vignette("diff")`](https://smith-chem-wisc.github.io/mzLibR/articles/diff.md)
and
[`vignette("ptmqtl")`](https://smith-chem-wisc.github.io/mzLibR/articles/ptmqtl.md);
the differential-abundance article gains a section for each new
statistic.

- **`stats_fit(contrasts = )`: any comparison, with a confidence
  interval** (mzLib
  [\#1452](https://github.com/smith-chem-wisc/mzLibR/issues/1452),
  `EmpiricalBayes.ModerateContrasts`). Give the contrasts file
  `diff_design(out = )` writes beside the design (`contrasts_<i>.tsv`),
  or your own: a contrast id column, then exactly the design’s
  coefficient columns, one weight each. Each contrast is moderated with
  the same prior as the coefficients, has its own Benjamini-Hochberg
  family, and carries `confidence_low` and `confidence_high`, as limma’s
  `contrasts.fit()`, `eBayes()` and `topTable(confint = )`. The answers
  are a second table, `fit$contrasts` (an `mzlibr_moderated_contrasts`
  with `confidence_level`, `tested_contrasts`, `weights` and `records`),
  so the coefficient table is unchanged; `coefficients` may now be empty
  and defaults to `character(0)`. `confidence` sets the level, as a
  fraction; left `NULL`, it is mzLib’s default (0.95).
  `fit$contrasts_file` names the file.
- **`stats_fit(estimator = )`: mzLib’s marginal-likelihood variance
  prior** (mzLib
  [\#1368](https://github.com/smith-chem-wisc/mzLibR/issues/1368)).
  `"marginal_likelihood"` fits the prior with each feature’s own
  residual df, so the unequal df that missing values cause are handled
  exactly; `"moments_legacy"` (limma’s `eBayes(legacy = TRUE)`) stays
  mzLib’s default, and the new `fit$prior$estimator` reports which was
  used. On the MALAT1 dilution series it lowers the prior df from 19.6
  to 11.0 and drops the false FLuc call at half load.
- **An older bridge is named, not trusted.** A bridge before pyMzLib
  0.6.0 lists `stats fit` but ignores `--contrasts` and `--estimator`;
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  then raises an `mzlib_usage_error` naming the 0.6.0 bridge instead of
  returning a fit that did not do what was asked.
- **[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)’s
  caveats name the variance-prior estimator.** The first caveat now
  begins `moments_legacy:` (or `marginal_likelihood:`), the
  `residual_df_differ` caveat names `estimator = "marginal_likelihood"`,
  and the `bh_adjusted` caveat says “per coefficient and per contrast”.
  Every number is unchanged; only code that matches caveat text needs a
  look. On fewer than 3 features `prior$scale` is `NA` and `prior$df` 0
  (mzLib
  [\#1398](https://github.com/smith-chem-wisc/mzLibR/issues/1398));
  [`?stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  now says bridges before pyMzLib 0.6.0, not 0.5.1, reported the first
  feature’s variance there.
- **[`stats_combine_p()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_combine_p.md)**:
  one p-value per feature from several independent studies, by Fisher’s
  method or Stouffer’s, weighted or not (mzLib
  [\#1367](https://github.com/smith-chem-wisc/mzLibR/issues/1367),
  `PValueCombination`). `method` is required. Give each test’s `effect`
  and every p-value is read as two-sided and converted to the one-sided
  p-value for effect \> 0 first, so opposite effects cancel. The tests
  hold the recordings to R’s
  [`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) and
  [`qnorm()`](https://rdrr.io/r/stats/Normal.html), recorded and
  computed here.
- **[`stats_correlate()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_correlate.md)**:
  Spearman’s rho of every feature with a per-sample trait, optionally
  partial given covariates (`given`; mzLib
  [\#1367](https://github.com/smith-chem-wisc/mzLibR/issues/1367),
  [\#1446](https://github.com/smith-chem-wisc/mzLibR/issues/1446),
  `SpearmanCorrelation`). Each row carries rho, a two-sided p-value, the
  samples used and how the p-value was computed (`"exact"`,
  `"asymptotic"` or `"not_estimable"`, with rho and p `NA`). The tests
  hold every recorded row to
  [`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) run here.
- **[`stats_fit_logistic()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_logistic.md)**:
  a logistic regression per feature of a 0/1 table, as R’s
  `glm(family = binomial)` (mzLib
  [\#1367](https://github.com/smith-chem-wisc/mzLibR/issues/1367)): is a
  feature’s detection related to the condition? A feature the design
  separates is `"separated"`, with `NA` statistics, instead of R’s huge
  coefficient. `max_iterations` and `tolerance` are mzLib’s when left
  `NULL`. The tests hold the recording to R 4.6.1’s
  [`glm()`](https://rdrr.io/r/stats/glm.html) output and to
  [`glm.fit()`](https://rdrr.io/r/stats/glm.html) run here.
- **[`stats_fit_mixed()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit_mixed.md)**:
  a linear mixed model per feature with a random intercept per group of
  samples, as `nlme::lme(random = ~ 1 | group)` (mzLib
  [\#1367](https://github.com/smith-chem-wisc/mzLibR/issues/1367)),
  tested by Wald t on nlme’s containment df. `estimator = "reml"` or
  `"ml"`; left `NULL`, mzLib’s default (REML), which the result reports.
  The tests hold the recordings to nlme 3.1.169’s output and to
  [`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html) run here.
- **[`stats_entrapment_fdp()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_entrapment_fdp.md)**:
  the lower bound, combined and paired estimates of a search’s false
  discovery proportion at each q-value threshold (mzLib
  [\#1385](https://github.com/smith-chem-wisc/mzLibR/issues/1385), Wen
  et al. 2025, as FDRBench). `ratio` is required; `paired` is reported
  only at a ratio of exactly 1; `thresholds` left `NULL` is mzLib’s
  default sweep, 0.001 to 0.100.
- **[`diff_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_design.md),
  with
  [`factor_spec()`](https://smith-chem-wisc.github.io/mzLibR/reference/factor_spec.md)
  and
  [`covariate_spec()`](https://smith-chem-wisc.github.io/mzLibR/reference/covariate_spec.md)**:
  a sample sheet to mzLib’s design matrix per stratum (mzLib
  [\#1450](https://github.com/smith-chem-wisc/mzLibR/issues/1450),
  `AnalysisDesign`), equal to R’s
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) (the
  tests run it here), with its rank, aliased terms, warnings and the
  contrasts it cannot run; one row per default contrast with **mzLib’s
  verdict on whether these samples can estimate it at all**. Each
  stratum’s `matrix` is an R matrix named by sample id, and `strata` is
  a list named by stratum key. With `out` it writes `design_<i>.tsv` and
  `contrasts_<i>.tsv`, which
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  reads.
- **[`diff_observations()`](https://smith-chem-wisc.github.io/mzLibR/reference/diff_observations.md)**:
  every peptide in every biological sample from a FlashLFQ or
  MetaMorpheus peptide table (mzLib
  [\#1455](https://github.com/smith-chem-wisc/mzLibR/issues/1455)): log2
  intensity, fractions summed and technical replicates averaged, or
  **`NA`, never 0, with the reason in `state`**. `runs` is a data.frame
  of `file` and `sample`, and optionally `fraction` and
  `technical_replicate`; `basis` has no default. mzLib’s shared-peptide
  normalization is always measured and reported, and subtracted with
  `normalize = TRUE`. `responses_out` writes the table
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  reads.
- **[`ptmqtl_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_occupancy.md),
  [`ptmqtl_pairs()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_pairs.md)
  and
  [`ptmqtl_trait_effects()`](https://smith-chem-wisc.github.io/mzLibR/reference/ptmqtl_trait_effects.md)**:
  PTM site occupancy per run with each cell’s state (a censored cell has
  no `fraction`, never a 0), physical and co-varying PTM pairs (`type`
  is required; `remove_run_level` takes out a batch shift), and a
  per-run trait’s effect on logit occupancy with a random intercept per
  biological replicate (mzLib
  [\#1369](https://github.com/smith-chem-wisc/mzLibR/issues/1369),
  [\#1430](https://github.com/smith-chem-wisc/mzLibR/issues/1430),
  [\#1446](https://github.com/smith-chem-wisc/mzLibR/issues/1446)).
  **Positions are 1-based, unlike
  [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)**.
- The vendored specs are re-synced from bridge `60b1d259`.
  `scripts/sync-specs.R` now reads a bare `y` or `n` in a spec as the
  string the bridge’s own YAML reader takes it for, where R’s `yaml`
  package read it as a boolean (`stats correlate`‘s and `ptmqtl pairs`’
  `n` column).

### Data: checksums, alignment, orthology, spectral libraries and SDRF scoping (the pyMzLib 0.6.0 bridge)

Ported from pyMzLib 0.6.0, verb for verb; each *needs the pyMzLib 0.6.0
bridge* and is refused, naming that release, by an older one. Every
option left out is not sent, so mzLib’s default applies and the result
says what was used. Every recording the help pages and articles replay
is pyMzLib’s, byte for byte. A new article,
[`vignette("library")`](https://smith-chem-wisc.github.io/mzLibR/articles/library.md);
the PRIDE, protein-database and SDRF articles gain a section for each
new verb.

- **[`pride_file_checksums()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_file_checksums.md):
  the MD5 and exact size PRIDE recorded for each deposited file** (mzLib
  [\#1427](https://github.com/smith-chem-wisc/mzLibR/issues/1427), wire
  `pride checksums`). A data.frame of `file_name`, `md5`, `size_bytes`
  and `project_accession`, its row names the file names, so
  `checksums["README.txt", ]` looks a file up by the name it was
  downloaded under. `size_bytes` is the exact download length, which
  neither other listing gives (PXD000001’s `.xml.gz`: 10,668,000 bytes,
  where
  [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  says 10,677,205), and `md5` is not the manifest’s `checksum` (a SHA-1
  where present). It is a third list: PRIDE’s `generated/` files have no
  row. Zero rows is an answer; an unknown accession, or a project not
  yet published, raises `mzlib_project_not_found`; a list mzLib cannot
  read stays an `mzlib_bridge_error` (`MzLibException`). No timeout of
  mzLibR’s own.
  [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and
  [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  still check no bytes, and now say so; the PRIDE article shows the
  check with [`file.size()`](https://rdrr.io/r/base/file.info.html) and
  [`tools::md5sum()`](https://rdrr.io/r/tools/md5sum.html).
- **[`proteins_align()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_align.md):
  map a residue onto an isoform, an ortholog or another release** (mzLib
  [\#1428](https://github.com/smith-chem-wisc/mzLibR/issues/1428), wire
  `proteins align`). A global affine alignment scored with BLOSUM62, for
  a list of `c(source, target)` pairs in one bridge call: `records` (one
  row per pair: score in half-bits, identical and substituted residues,
  residues facing a gap, the aligned rows) and `residue_map` (one row
  per residue, its counterpart’s position and letter, `NA` facing a gap,
  never the nearest residue). `pair_index` is 1-based.
  [`proteins_map_source_to_target()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_source_to_target.md)
  and
  [`proteins_map_target_to_source()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_map_target_to_source.md)
  look positions up, and
  [`proteins_alignments()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_alignments.md)
  splits a result into one alignment per pair. `gap_open`, `gap_extend`,
  `end_gaps` and `max_cells` default to `NULL`, mzLib’s (11, 1,
  `"free"`, 200,000,000); `residue_map = FALSE` sends
  `--no-residue-map`. One pair is
  `proteins_align(list(c(source, target)))`, as the spec spells it:
  there is no separate one-pair function, which is pyMzLib’s
  `proteins.align()`.
- **[`proteins_parse_accessions()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_parse_accessions.md):
  entry, isoform or version, and applied variants** (mzLib
  [\#1382](https://github.com/smith-chem-wisc/mzLibR/issues/1382), wire
  `proteins parse-accession`). `P04637-2_S70N` is entry `P04637`,
  isoform 2, with `S70N`; `NP_000537.3_R72P` is RefSeq `NP_000537`,
  version 3. It parses and never repairs: a FASTA header, a `DECOY_`
  prefix, a `|`-joined group or mzLib’s `P04637_2` collision counter is
  `unrecognized` and kept verbatim, and every accession is sent exactly
  as given, not trimmed.
- **[`proteins_build_orthology()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_build_orthology.md):
  an Ensembl Compara orthology snapshot** (mzLib
  [\#1381](https://github.com/smith-chem-wisc/mzLibR/issues/1381), wire
  `genes orthology`). From each species’ GTF (or gene table), every
  species’ homology dump - mzLib refuses a call without one - and the
  gene-tree file, all local; `gtfs` and `gene_sets` are named character
  vectors, production name to path. It reports what went in
  (`gene_sets`, `gene_trees`, `dumps`), the relationships per species
  pair (`pairs`;
  [`proteins_orthology_pair()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_orthology_pair.md)
  finds one in either order), a `snapshot_id` to record with any result,
  with `pair_status` every gene’s status in one other species
  (`has_ortholog`, `no_edge_in_shared_tree`,
  `tree_lacks_target_species`, `not_in_any_tree`), and with `out` the
  Parquet snapshot it wrote (`written`). mzLibR reads no Parquet; the
  proteins article queries it with `arrow`. mzLib’s three-species test
  world is copied into `tests/orthology/`, kept byte-exact by
  `.gitattributes` because the recordings’ sha256s are over its bytes.
- **[`library_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read.md)
  and
  [`library_read_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/library_read_many.md):
  mzLib’s `.msl` spectral libraries** (wire `library read`). One row per
  precursor, or per fragment ion with `fragments = TRUE`, in the
  library’s stored order; `retention_time_unit` is per row (iRT or
  minutes, never converted), and mzLib’s “not available” values (ion
  mobility 0.0, a NaN q-value, NCE 0, an empty protein) are `NA`.
  `precursor_index` is 1-based. `index_only_threshold` defaults to
  `NULL`, mzLib’s 1 GiB; `limit` and `offset` count precursors. The
  many-libraries form returns one long table and a `files` data.frame,
  in one bridge call.
- **[`sdrf_isobaric_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_isobaric_design.md):
  a TMT, TMTpro, iTRAQ or DiLeu design from an SDRF** (wire
  `sdrf isobaric-design`). One row per file and annotated channel, or
  every reason there is none. It guesses neither the kit nor the plex:
  declare exactly one of `plex_column`, `plex_file_name_pattern` and
  `single_plex`. **Its numbers are 1-based, as the SDRF writes them**,
  where
  [`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)’s
  are 0-based; the split is mzLib’s, and mzLibR converts neither. `out`
  names a folder for MetaMorpheus’s `TmtDesign.txt`.
- **[`sdrf_restrict()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_restrict.md):
  an SDRF cut to the files one search read** (wire `sdrf restrict`).
  mzLib joins `X-calib.mzML` to the row naming `X.raw` by stem, adds
  `comment[searched data file]`, and reports what it could not place
  (`searched_without_row`, `dropped_data_files`, `ambiguous`,
  `shared_stem`). Name an element of `searched_files` by its acquired
  file to override the suffix rule. The result is an `mzlibr_sdrf`;
  [`sdrf_searched_data_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_searched_data_files.md)
  reads the new column. Restrict first, then read a design: both design
  readers key on that column.

*Same file, different value.*

- **`sdrf_pool(out = )`’s `written$path` names the file mzLib wrote.**
  mzLib adds `.sdrf.tsv` to a name that does not already end in it, so
  `out = "x.tsv"` writes `x.tsv.sdrf.tsv`. With the pyMzLib 0.6.0 bridge
  `written$path` says `x.tsv.sdrf.tsv`; older bridges reported `x.tsv`,
  a file that does not exist. The file written is the same, and an `out`
  already ending in `.sdrf.tsv` is unchanged. `sdrf_restrict(out = )`
  reports its file the same way.

*For contributors.* The replay bridge (`R/replay.R`) gains pyMzLib
0.6.0’s rules for these verbs: a recording answers only a call with the
same pairs, accessions, searched files or species on stdin, the same
`--sequences`, `--no-residue-map`, `--pair-status` and
`--restrict-gene-trees` choice, and the options it was recorded with;
and, as in pyMzLib, a recording whose echo of an option is null never
answers a call that sends it. The test harness gains
`recorded_bridge_call()`, which runs a public function through the real
transport and verb check against a recording and returns what it sent.

### Defaults come from mzLib (the pyMzLib 0.5.0 bridge)

- **[`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  takes its defaults from mzLib, so a call that names no protease or
  dissociation type gives different peptides and ions.** Every digestion
  argument (`protease`, `dissociation`, `missed_cleavages`,
  `min_length`, `max_length`, `max_modifications`, `max_isoforms`,
  `terminus`) now defaults to `NULL`, mzLibR sends nothing for it, and
  mzLib’s `DigestionParams` default applies. Two of the old defaults
  were mzLibR’s own and change:
  - protease: mzLib’s default is `"trypsin"`, which also cuts before
    proline. The old default was `"trypsin|P"`, which does not.
    Albumin’s default digest goes from 303 to 310 peptidoforms, and from
    195 to 202 distinct sequences.
  - dissociation: mzLib names no default, so the bridge uses HCD (b and
    y ions) until it does (smith-chem-wisc/mzLib#1488). The old default
    was ETD (c and z-dot ions).

  The other defaults (missed cleavages, the length window, the
  modification and isoform caps, terminus) have the same values as
  before, now set by mzLib. **To get the old answer, pass
  `protease = "trypsin|P", dissociation = "ETD"`.** A digest records
  `min_length`, `max_length` and `missed_cleavages` as `NA` when they
  were left to mzLib, since the bridge does not report those three back,
  and prints them as mzLib’s default. A `max_length` below `min_length`
  is refused only when both are given. The help page’s counts are
  re-measured at mzLib `d10017ce`.
- **[`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  and
  [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
  leave the page size to mzLib.** `page_size` now defaults to `NULL` and
  is sent only when given, so mzLib’s own paging default applies (100
  today, as before). The results never depended on it. A given page size
  is checked as before.
- **[`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
  [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
  and
  [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
  add no timeout of their own.** They defaulted to `timeout = 300`,
  which turned a PRIDE outage that hangs into an `mzlib_timeout`. mzLib
  already bounds each request
  ([\#1412](https://github.com/smith-chem-wisc/mzLibR/issues/1412):
  three retries after 5, 20 and 60 s, each attempt up to 100 s), so they
  now default to `NULL` and a hang is reported as
  `mzlib_service_unavailable`, after about 8 minutes per request. Pass
  `timeout = 300` for the old limit; a number still cuts the call off
  with `mzlib_timeout`. The PRIDE article’s outage handling now waits
  minutes and tries once more, since mzLib has already retried, and
  [`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  says so.
- [`?stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  no longer says the bridge’s spline default is 4. With `trend = TRUE`
  and no `spline_basis`, mzLibR has always sent none; the 0.5.0 bridge
  now asks mzLib for none either, so mzLib sizes the trend as limma does
  and a trend fit on fewer than 30 features gets limma’s numbers.
  `prior$spline_basis_count` says what was used.
- The help pages are re-rendered from the bridge’s specs, re-checked
  against mzLib `d10017ce`. PXD000001’s REST manifest listed 9 of the 14
  files in its FTP tree on 2026-10-21, and
  [`?pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
  says so rather than 8 of 13.

### mzLib 1.0.595

Built from mzLib 1.0.595 (tag commit `d10017ce`), which carries mzLib
1.0.594 as well.

- **DIA-NN 2.x’s `report.parquet` reads** (mzLib
  [\#1384](https://github.com/smith-chem-wisc/mzLibR/issues/1384), file
  type `DiaNnReportParquet`). DIA-NN 2.x writes its main report only as
  parquet. mzLib reads it into the same record as the 1.8/1.9
  `report.tsv`, so it is a fifth quantifiable type:
  [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md),
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  and
  [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  take it.
  **[`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  lists 39 types, up from 38, and 5 offer the `quantifiable` view, up
  from 4**; 19 still offer none. Like the TSV report it is recognised by
  its columns (`Precursor.Id`, `Stripped.Sequence`, `Run`), not its
  name, and any other `.parquet` is an `mzlib_usage_error`. **Its
  `retention_time_unit` is `"minutes"`**, as for the TSV report.
  `absent_fields` is always empty for it, because there is no text
  header to judge from, so the columns 2.x stopped writing (`File.Name`,
  `PG.Quantity`, `MS2.Scan` and others) are `NA` in
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  without being named there, and `ms2_scan_number` reads 0.
- **`peptidoform_convert(target = "ProForma")` names UniProt
  modifications by UNIMOD accession** (mzLib 1.0.594,
  [\#1403](https://github.com/smith-chem-wisc/mzLibR/issues/1403)).
  `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]-SEQK`, where
  it was written back by name; the caveat that said so is retired, so a
  result carries five caveats, not six. The `pro_forma` column
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  computes for a `.psmtsv` written without one gains the same
  accessions.

Every other result mzLib 1.0.594 and 1.0.595 changed, listed after the
fact: pyMzLib measured each with the bridge at `d10017ce`, and the
conversions, the
[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
prior and the `entrapment_fraction` column were re-checked from R
against the pyMzLib 0.5.0 bridge. The help pages carry them too,
rendered from the bridge’s specs.

*Calls that used to succeed are now refused.*

- **ProForma input that lost part of itself** (mzLib
  [\#1462](https://github.com/smith-chem-wisc/mzLibR/issues/1462)).
  `peptidoform_convert(..., source = "ProForma")` converted
  `PEPTIDE-[UNIMOD:2][UNIMOD:35]` with one C-terminal modification
  dropped, and read `PEPTIDE-[UNIMOD:2]K` as `PEPTIDEK`. Stacked
  terminal modifications now fail with `failure_reason`
  `"IncompatibleModifications"` (`incompatible_items`: “Multiple
  C-terminal modifications (2 found).”), and a residue after a
  C-terminal modification, a leading `-` or a modification with no
  sequence fail with `"InvalidSequence"`.
- **SUMO remnants named `... on D`** (mzLib
  [\#1432](https://github.com/smith-chem-wisc/mzLibR/issues/1432)).
  mzLib’s human SUMO-1 and SUMO-2/3 remnants
  (`DVIEVYQEQTGG (SUMO-1 Site human)`,
  `DVFQQQTGG (SUMO-2/3 Site human)`) now target K, where they
  targeted D. A full sequence that names the `on D` form no longer
  converts; the `on K` form does.
- **TMTpro on S, T, Y or H** (mzLib
  [\#1419](https://github.com/smith-chem-wisc/mzLibR/issues/1419)).
  mzLib used to turn these into TMT, 224 Da, so
  `PEPS[Common Fixed:TMTpro on S]K` converted at the wrong mass. They
  now fail to convert. TMTpro on K is unchanged (+304.2071 Da).

*Calls that used to fail now succeed.*

- **[`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  with gapped replicate numbers** (mzLib
  [\#1422](https://github.com/smith-chem-wisc/mzLibR/issues/1422)). A
  design whose biological replicates skip a number, with fractions,
  raised an error; it now quantifies. So does the numbering the bridge
  uses when `spectra` gives no `biological_replicate`, which skips at
  every condition after the first. `bayesian_protein_quant = TRUE` on
  such a design now writes `BayesianProteinQuant.tsv` to
  `output_directory`.
- **A condition whose name ends in digits** (mzLib
  [\#1424](https://github.com/smith-chem-wisc/mzLibR/issues/1424)).
  Condition `A1` biorep 0 beside condition `A` biorep 10 raised an
  `mzlib_bridge_error` in
  [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  and
  [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md),
  because the two spell the same text. It now quantifies.

*Same file, different value.*

- **Entrapment and contaminant flags read every joined label** (mzLib
  [\#1371](https://github.com/smith-chem-wisc/mzLibR/issues/1371)). In
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  and
  [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
  `is_decoy`, `is_contaminant` and `is_entrapment` are `TRUE` when any
  `|`-joined parent is D, C or E. A protein group labelled `T|C` is now
  a contaminant and `T|ET` entrapment; on a `.psmtsv` the “E” in
  MetaMorpheus’s “Output too long for Excel” no longer reads as
  entrapment. Labels MetaMorpheus writes itself read as before.
- **[`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  of a `.psmtsv` or `.osmtsv` gains an `entrapment_fraction` column**
  (mzLib
  [\#1371](https://github.com/smith-chem-wisc/mzLibR/issues/1371)): the
  share of each PSM that is an entrapment discovery, from 0 to 1. A
  shared `T|ET` counts half when the sequences differ and nothing when
  both proteins carry the same peptide. Use it, not `is_entrapment`, to
  count entrapment discoveries. It is `NA` where mzLib cannot line the
  label up with the full sequence (mzLib’s `NaN`).
- **Lactylation on K weighs one hydrogen more** (mzLib
  [\#1417](https://github.com/smith-chem-wisc/mzLibR/issues/1417)).
  mzLib’s formula is now C3H4O2, 72.021129 Da, where it was C3H3O2,
  71.013304 Da. `peptidoform_convert(..., target = "MassShift")` of
  `PEPK[Less Common:Lactylation on K]R` gives `PEPK[+72.021129]R`, and
  the `"ProForma"` and `"Unimod"` targets write `UNIMOD:2114`. A
  `.psmtsv`’s `pro_forma` column changes the same way.
- **`flashlfq_quantify(normalize = TRUE)`** (mzLib
  [\#1422](https://github.com/smith-chem-wisc/mzLibR/issues/1422)). The
  reference is now the first condition’s lowest biological replicate,
  and each fraction’s lowest technical replicate, taken in numeric
  order. A design numbered from 1 is now normalized; before, biorep
  normalization did nothing for it. A gapless design numbered from 0 and
  listed in order is unchanged.
- **[`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  on fewer than 3 features** (mzLib 1.0.594,
  [\#1398](https://github.com/smith-chem-wisc/mzLibR/issues/1398)).
  mzLib fits no prior, as limma does: `prior$df` is 0, each feature
  keeps its own variance, and t and p are ordinary per-feature t-tests.
  They used to borrow a prior fitted to the two features. The 0.5.0
  bridge still reports the first feature’s own variance as `prior$scale`
  there; bridges from pyMzLib 0.5.1 report `NA`, and each feature’s
  value is in the `prior_variance` column.
  [`?stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  says both. Tables of 3 or more features are unchanged.
- **[`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  with a `|P` protease** (mzLib
  [\#1320](https://github.com/smith-chem-wisc/mzLibR/issues/1320)).
  `protease = "trypsin|P"` cut a protein that begins `RP` or `KP` after
  its first residue; it no longer does, and the other `|P` proteases
  change the same way. mzLib’s default `"trypsin"`, and any protein that
  begins with M, are unaffected.
- **PRIDE: a failure that can pass is retried** (mzLib 1.0.594,
  [\#1412](https://github.com/smith-chem-wisc/mzLibR/issues/1412)).
  Every `pride_` function retries such a failure three times, after 5,
  20 and 60 seconds, so an outage takes longer to report. An FTP listing
  that comes back empty, which
  [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
  used to return as no files (and so `mzlib_project_not_found`), is
  retried and then raised as an `mzlib_bridge_error` with `error_type`
  `"HttpRequestException"`. A failed download leaves `<name>.partial`
  and `<name>.partial.validator` behind, so the next call resumes.
- **PRIDE downloads are not checked yet.** mzLib can check a download
  against PRIDE’s checksums and catch an EBI error page served at the
  right size
  ([\#1427](https://github.com/smith-chem-wisc/mzLibR/issues/1427),
  [\#1449](https://github.com/smith-chem-wisc/mzLibR/issues/1449)), but
  not on the path
  [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and
  [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  use, so a file can be reported as downloaded when its bytes are wrong.
  [`?pride_download`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and the PRIDE article say to check large files yourself.
- [`?stats_fit`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  no longer says contrasts are not implemented: mzLib has them
  (`EmpiricalBayes.ModerateContrasts`,
  [\#1452](https://github.com/smith-chem-wisc/mzLibR/issues/1452)), but
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  does not expose them yet.
- The psmtsv caveat on
  [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md)
  cites `SpectrumMatchFromTsv.cs:122` and `:197`; mzLib
  [\#1371](https://github.com/smith-chem-wisc/mzLibR/issues/1371) moved
  the lines it cited. The claim is unchanged.

### mzLib 1.0.593

- [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
  [`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  and
  [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
  read an RNA search’s `AllQuantifiedTranscriptGroups.tsv` and
  `AllQuantifiedOligos.tsv` (mzLib
  [\#1388](https://github.com/smith-chem-wisc/mzLibR/issues/1388)).
  mzLib reads them with subclasses of the protein-group and peptide
  readers, so the columns keep their protein names: `protein_group_name`
  is the transcript group, `sequence` the oligonucleotide, and occupancy
  names RNA modifications.
- MetaMorpheus protein-group tables written without quantification now
  read (mzLib
  [\#1365](https://github.com/smith-chem-wisc/mzLibR/issues/1365)):
  `AllProteinGroups.tsv` and each file’s `<file>_ProteinGroups.tsv`
  failed with `Tsv file type not supported`, and now read as
  `MetaMorpheusQuantifiedProteinGroups`, with `intensity` in
  `absent_fields`.
- Four shipped modifications now write their Unimod accession in
  `pro_forma` (mzLib
  [\#1328](https://github.com/smith-chem-wisc/mzLibR/issues/1328)):
  `GG (Ubiquitination Site)` as `[UNIMOD:121]`, both `Myristoylation`
  entries as `[UNIMOD:45]` and `EQIGG` as `[UNIMOD:846]`, instead of by
  name. Masses are unchanged.
- A read fault on an existing `.mzid` is an `mzlib_bridge_error` naming
  the file (mzLib
  [\#1362](https://github.com/smith-chem-wisc/mzLibR/issues/1362)),
  where it was an `IOException`. A missing file is still an
  `mzlib_usage_error`.
- PRIDE download errors name the file and the host, never the URL (mzLib
  [\#1350](https://github.com/smith-chem-wisc/mzLibR/issues/1350)), so a
  reviewer token in a query string cannot reach a log. They are still
  `mzlib_service_unavailable`.
- The help pages no longer state how many file types mzLib recognises:
  [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  asks mzLib, and a number written down goes stale at the next release.

### Reading files

- [`readers_formats()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_formats.md)
  and
  [`readers_identify()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify.md)
  say what a file is and which cross-format views it offers.
- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  reads scans from mzML, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF
  and msalign: headers always, peaks on request.
- [`readers_read_results()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_results.md),
  [`readers_read_features()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_features.md)
  and
  [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  read the three cross-format views of search results;
  [`readers_read_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_records.md)
  reads any recognised file into that format’s own fields.
- [`readers_retention_time_in_minutes()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_retention_time_in_minutes.md)
  converts using each file’s declared unit, and refuses when there is
  none - for a `_many` batch, each file’s rows by that file’s unit.
- Many files in one call:
  [`readers_identify_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_identify_many.md)
  and a `_many` form of every reader hand the whole list to one bridge
  process with `threads` and `on_error` stated, and return one long
  `records` table plus a `files` data.frame of per-file facts, one row
  per input (*needs the pyMzLib 0.2.0 bridge*).
- Every read reports the per-file block: `reader`, `rows_not_read`,
  `retention_time_unit`, `caveats`, and the three ways a column can be
  empty kept apart - `absent_fields` (the format has no source for it,
  so it is `NA` in every row), `failed_fields` and `excluded_fields`.
- [`readers_read_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_spectra.md)
  reports the run’s `source`: instrument model and its PSI-MS accession,
  serial number, and acquisition start time with whether it is UTC
  (mzLib
  [\#1349](https://github.com/smith-chem-wisc/mzLibR/issues/1349)).
- [`readers_read_matches()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_matches.md)
  reads mzIdentML (`.mzid`, `.mzid.gz`) with each match’s `q_value`,
  `rank` and `pass_threshold`, and with `scores = TRUE` one row per
  engine score; mzIdentML items mzLib does not represent are listed in
  `skipped` with the reason.
- New
  [`readers_read_protein_groups()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_protein_groups.md),
  [`readers_read_quantified_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_quantified_peptides.md)
  and
  [`readers_read_occupancy()`](https://smith-chem-wisc.github.io/mzLibR/reference/readers_read_occupancy.md)
  read MetaMorpheus and FlashLFQ quantification tables long, one row per
  record per sample, and per-site PTM occupancy (mzLib
  [\#1347](https://github.com/smith-chem-wisc/mzLibR/issues/1347);
  *needs the pyMzLib 0.2.0 bridge*).
- [`mzlibr_bridge_version()`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md)
  reports `verbs`, every command the bridge dispatches.
- mzLib now recognises Pytheas, mzIdentML and its gzipped form, and the
  MetaMorpheus and FlashLFQ quantification tables.
- `limit = 0` is accepted by every reader, as the bridge allows: the
  envelope alone.

### SDRF-Proteomics

- [`sdrf_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_read.md)
  reads one experimental-design file and
  [`sdrf_pool()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_pool.md)
  pools several into one analysis table with provenance;
  [`sdrf_value()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_value.md),
  [`sdrf_all()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_all.md),
  [`sdrf_records()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_records.md)
  and the other `sdrf_` accessors read them.
- [`sdrf_validate()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate.md),
  [`sdrf_lint()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_lint.md),
  [`sdrf_assess()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess.md),
  [`sdrf_samples()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples.md)
  and
  [`sdrf_parse_ages()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_parse_ages.md)
  ask mzLib whether an SDRF document is well-formed, whether several
  agree on their vocabulary, whether it describes its samples, and what
  its ages are in years;
  [`sdrf_validate_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_validate_many.md),
  [`sdrf_assess_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_assess_many.md)
  and
  [`sdrf_samples_many()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_samples_many.md)
  read a corpus in one bridge call (mzLib
  [\#1207](https://github.com/smith-chem-wisc/mzLibR/issues/1207),
  [\#1325](https://github.com/smith-chem-wisc/mzLibR/issues/1325),
  [\#1326](https://github.com/smith-chem-wisc/mzLibR/issues/1326),
  [\#1333](https://github.com/smith-chem-wisc/mzLibR/issues/1333),
  [\#1335](https://github.com/smith-chem-wisc/mzLibR/issues/1335)).
- New
  [`sdrf_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design.md)
  reads the label-free experimental design MetaMorpheus and FlashLFQ
  take out of an SDRF, with mzLib’s `SdrfLabelFreeDesign`, or every
  reason it refuses to - all at once, and as a result rather than an
  error.
  [`sdrf_design_spectra()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_spectra.md)
  and
  [`sdrf_design_run_design()`](https://smith-chem-wisc.github.io/mzLibR/reference/sdrf_design_run_design.md)
  hand it to
  [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  and
  [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md);
  `out` writes MetaMorpheus’s `ExperimentalDesign.tsv`. Its replicate
  and fraction coordinates stay 0-based, as the quant functions take
  them (mzLib
  [\#1363](https://github.com/smith-chem-wisc/mzLibR/issues/1363);
  *needs the pyMzLib 0.3.0 bridge*).

### Isobaric kits

- New
  [`isobaric_kits()`](https://smith-chem-wisc.github.io/mzLibR/reference/isobaric_kits.md)
  lists mzLib’s TMT, TMTpro, iTRAQ and DiLeu kits with every channel’s
  label, theoretical reporter-ion m/z and matching window, or one kit by
  MetaMorpheus’s name for it; `channel_index` is 1-based (mzLib
  [\#1375](https://github.com/smith-chem-wisc/mzLibR/issues/1375);
  *needs the pyMzLib 0.3.0 bridge*).

### PRIDE, peptidoforms and quantification

- [`pride_list_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
  [`pride_list_ftp_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md),
  [`pride_download()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download.md)
  and
  [`pride_download_files()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_download_files.md)
  list and fetch PRIDE Archive projects.
- [`pride_search()`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_search.md)
  finds PRIDE projects by keyword, one row per project, with dates as
  `Date`.
- [`peptidoform_fragments()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_fragments.md)
  digests a UniProt entry and fragments its peptidoforms.
- New
  [`peptidoform_convert()`](https://smith-chem-wisc.github.io/mzLibR/reference/peptidoform_convert.md)
  rewrites full sequences in another notation with mzLib’s
  `SequenceConversionService`. The main use is MetaMorpheus full
  sequences to Unimod accessions: `[UniProt:N-acetylserine on S]SEQK`
  becomes `[UNIMOD:1]SEQK`, and dimethyllysine becomes `UNIMOD:36`. One
  row per input, in order, in `records`, with mzLib’s own verdict
  (`converted`, `converted_with_warnings` or `failed`), the
  modifications it could not write, and its warnings. `source` and
  `target` take any notation mzLib has registered; every result lists
  them. `mode` is mzLib’s `SequenceConversionHandlingMode`. The ProForma
  target names UniProt modifications by UNIMOD accession too, since
  mzLib [\#1403](https://github.com/smith-chem-wisc/mzLibR/issues/1403)
  (see mzLib 1.0.595 above; at mzLib 1.0.593 it wrote them back by name,
  mzLib#1401). (pyMzLib
  [\#75](https://github.com/smith-chem-wisc/mzLibR/issues/75); *needs
  the pyMzLib 0.4.0 bridge*.)
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  runs FlashLFQ label-free quantification with match-between-runs.
- [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  re-runs FlashLFQ’s protein median polish from a
  `QuantifiedPeptides.tsv` under a new experimental design, without
  re-reading spectra.
- [`flashlfq_quantify()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_quantify.md)
  no longer warns about mzLib#1111 when `max_threads` is not 1:
  mzLib#1155 fixed it inside the pinned bridge. The default stays 1 here
  (pyMzLib and mzLibRust use -1).
- [`flashlfq_median_polish()`](https://smith-chem-wisc.github.io/mzLibR/reference/flashlfq_median_polish.md)
  with no design no longer hangs in terminal R. A call that sends the
  bridge nothing on stdin now gives it an empty stdin rather than R’s
  own, which a terminal never closes. Every verb goes through the same
  transport, so none can wait on it (pyMzLib
  [\#73](https://github.com/smith-chem-wisc/mzLibR/issues/73)).

### Differential abundance

- New
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md)
  fits one linear model per feature of a feature-by-sample table and
  tests each named coefficient with limma’s empirical-Bayes moderated t,
  Benjamini-Hochberg adjusted: `lmFit()` then `eBayes(legacy = TRUE)`,
  computed by mzLib’s `LinearModel` and `EmpiricalBayes`. A feature that
  cannot be fitted is reported with its reason, never dropped, and
  `residual_df_differ` says when default limma would use a different
  prior estimator.
- New
  [`stats_adjust()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_adjust.md)
  Benjamini-Hochberg adjusts p-values from anywhere, keeping every
  position; `NA` is untested and not counted in m.
- New
  [`stats_meta()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_meta.md)
  pools one estimate per study into a DerSimonian-Laird random-effects
  estimate per feature, with direction agreement and leave-one-out
  sensitivity.
- The tests hold all three to limma, metafor and
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) to 1e-8
  relative - against their recorded output always, and against limma and
  metafor run in the test when they are installed (both are now in
  `Suggests`). The reference tables ship in `extdata/stats/` so the
  examples and the new article read real files (mzLib
  [\#1341](https://github.com/smith-chem-wisc/mzLibR/issues/1341),
  [\#1357](https://github.com/smith-chem-wisc/mzLibR/issues/1357);
  *needs the pyMzLib 0.3.0 bridge*).

### Protein databases

- New
  [`proteins_read()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_read.md)
  reads UniProt XML or FASTA databases into one row per protein -
  organism, NCBI taxon, gene names, length, mass - with GO terms and
  Ensembl gene links on request;
  [`proteins_resolve_genes()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_resolve_genes.md)
  resolves proteins to stable Ensembl gene ids against a gene set you
  pin;
  [`proteins_classify_peptides()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_classify_peptides.md)
  calls each peptide `Unique`, `SharedWithinGene`, `SharedAcrossGenes`
  or `NotInDatabase`, with I and L one residue. One database or many, in
  one bridge call (mzLib
  [\#1336](https://github.com/smith-chem-wisc/mzLibR/issues/1336),
  [\#1338](https://github.com/smith-chem-wisc/mzLibR/issues/1338),
  [\#1348](https://github.com/smith-chem-wisc/mzLibR/issues/1348);
  *needs the pyMzLib 0.2.0 bridge*).
- New
  [`proteins_annotate_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_annotate_go.md)
  annotates a stored MetaMorpheus protein-group table with Gene Ontology
  terms: one row per (group, term) that any member holds, directly or
  through an ancestor, naming the members that carry it, so consensus
  and direct-only views are filters on the rows. Every non-decoy group
  gets a row, a term-less one saying why. It reads the go.obo you name
  and never downloads one;
  [`proteins_update_go()`](https://smith-chem-wisc.github.io/mzLibR/reference/proteins_update_go.md)
  is the one function that fetches the current release, keeping the
  previous file beside it (mzLib
  [\#1353](https://github.com/smith-chem-wisc/mzLibR/issues/1353),
  [\#1366](https://github.com/smith-chem-wisc/mzLibR/issues/1366);
  *needs the pyMzLib 0.3.0 bridge*).

### Documentation

- Every help page for a wire verb now carries that verb’s facts,
  generated from the bridge’s per-verb spec rather than retyped: what it
  wraps in mzLib, every parameter with its unit, range and default,
  every returned field and column with its unit and what `NA` means, the
  errors it raises as condition classes, caveats, the same verb in
  pyMzLib and mzLibRust, references and `since`. `scripts/name-parity.R`
  checks the hand-written `@param` and `@return` text against the same
  specs, and CI fails on any difference.
- Every exported function has an example, and `R CMD check` runs them
  all: a stand-in bridge answers each call from a recording of the real
  one, so the examples need no bridge, .NET or network.
- New
  [`?mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md)
  documents the condition classes and how to handle each.
- A pkgdown site, with an article for each module.
- The articles are now vignettes, so `R CMD check` builds and runs every
  one of them, on every platform, against the replay bridge. Each opens
  with a question -\> function -\> mzLib table, ends with what to cite
  (rendered from the specs’ DOIs, each checked to resolve), and states
  no count of formats or verbs and no mzLib version in its prose;
  `scripts/docs-lint.R` holds all of that in CI.
- The articles teach on real data recorded from the real bridge: a
  FlashLFQ run with match-between-runs on mzLib’s K562 pair (and median
  polish reproducing its protein intensities for every group), the whole
  albumin digest, PXD000001’s live FTP listing, an RNA search’s
  transcript groups, and the MALAT1 dilution series for
  [`stats_fit()`](https://smith-chem-wisc.github.io/mzLibR/reference/stats_fit.md).
- Corrected from those runs: the K562 peptide roll-up shows 21 of the
  140 MBR transfers, not 52; `use_pep_q_value` filters nothing, it
  changes the q-value FlashLFQ carries; PXD000001’s FTP tree now holds
  14 files.
- `knitr` and `rmarkdown` join `Suggests`, for the vignettes. The
  package still imports nothing but base R.
