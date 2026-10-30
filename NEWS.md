# mzLibR 0.1.0 (unreleased)

The first release. Versions follow each wire verb's `since.mzlibr` in the bridge's verb specs,
which every help page's **Since** section renders.

This release projects mzLib 1.0.596. The verbs marked *needs the pyMzLib 0.2.0 bridge*, *needs
the pyMzLib 0.3.0 bridge*, *needs the pyMzLib 0.4.0 bridge*, *needs the pyMzLib 0.6.0 bridge* or
*needs the pyMzLib 0.7.0 bridge* below are new in the bridge that release publishes; with an older bridge they are refused, naming
the release they need. `mzlibr_install_bridge()` installs the bridge pyMzLib 0.7.0 publishes
(mzLib 1.0.596, `cc67b7a2`), verified against that release's SHA256SUMS.

## mzLib 1.0.596 (the pyMzLib 0.7.0 bridge)

`mzlibr_install_bridge()` now installs the pyMzLib 0.7.0 bridge, built from mzLib 1.0.596 (tag
commit `cc67b7a2`), up from 1.0.595. Several existing functions give a different answer, most
because mzLib changed them; each is listed below, as pyMzLib 0.7.0's upgrading notes list them.
Every recording the help pages and articles replay was re-recorded from the 0.7.0 bridge (pyMzLib's,
byte for byte), and the live tests hold each to what that bridge returns.

* **`flashlfq_quantify()`'s `max_threads` now defaults to `-1` (every core), not 1**, as in mzLib,
  the bridge, pyMzLib and mzLibRust. The answer does not change. The default was 1 while
  mzLib#1111 made multithreaded match-between-runs nondeterministic; mzLib#1155 fixed that, and
  mzLib's K562 pair re-measured identical: five runs at `-1` and one at `1` gave the same results.
  Pass `max_threads = 1` to keep the old setting; on a two-run pair it is a little faster.

*New, each needing the pyMzLib 0.7.0 bridge.* Ported from pyMzLib 0.7.0, verb for verb; every
number is mzLib's, and every recording the help pages and articles replay is pyMzLib's, byte for
byte. Every option left out is not sent, so mzLib's default applies and the result reports it.

* **`proteins_map_residues()`: a residue across species, through Ensembl Compara's gene-tree
  alignment** (wire `genes map-residues`, mzLib #1465, `ComparaGeneTreeAlignment`). For each
  `(source protein, position, target protein)` query - a data.frame or a list of
  `list(source, position, target)` - the residue of the target in the same alignment column, or
  mzLib's outcome saying why there is none: `"GapInTarget"` (never the nearest residue),
  `"DifferentAlignments"`, `"SourceNotInAlignment"`, `"TargetNotInAlignment"`. It reads the
  release's alignment file, its gene-tree file and your species' GTFs (`gtfs`, `gene_sets`, as
  `proteins_build_orthology()` takes them), all local; with none, mzLib keeps every gene and a
  caveat says so. `proteins = TRUE` returns the kept Ensembl proteins with their sequences, to chain
  with `proteins_align()` from and back to UniProt positions; `alignment_index` is 1-based in R. A
  position past the end of its source protein fails the whole call with mzLib's message, and
  mzLibR does not check lengths itself. The proteins article gains the section, with the
  UniProt-to-UniProt chain for human p53's Ser15; the inputs are a verbatim cut of Ensembl release
  116 in `tests/map_residues/`.
* **`library_calibrate_rt()`: a run's minutes onto a spectral library's iRT scale, and back**
  (wire `library calibrate-rt`, mzLib #1387, `IrtCalibration`). From anchors you choose (a
  data.frame of `group`, `rt_minutes`, `irt`), each group (run) is fitted on its own with mzLib's
  robust LOWESS (or `kind = "linear"`), and the queries sent with them are converted either way:
  where a library precursor should elute in this run, or an observed minute on the library's scale,
  with `is_extrapolated` outside the anchors. `minimum_anchors`, `bandwidth`, `knots` and
  `robustness_iterations` are mzLib's `IrtCalibrationOptions`, each sent only when given (mzLib's
  defaults at the pin: 20 anchors, 0.3, 100 points, 3 iterations). `on_error = "skip"` reports a
  run mzLib refuses in its row of `groups` (`error_kind`, `error_type`, `error_message`) instead of
  failing. The fitted curve cannot be returned - mzLib keeps its points private - so queries travel
  in the same call as their anchors. The spectral-library article gains the section, on two DIA-NN
  reports from mzLib's test data, shipped in `extdata/calibration/`.
* **`stats_fit_latent_mixed()`: a random-intercept mixed model per feature on a latent scale, from
  samples measured on it, counted through it (k of n trials), or both** (wire
  `stats fit-latent-mixed`, mzLib #1469, `LogitNormalMixedModel.Fit`). A PTM site's occupancy is
  the usual case. Counts are read through the detection `link` you choose - `"identity"`,
  `c(intercept, slope[, depth_slope[, overdispersion]])`, or `stats_fit_link()`'s result -
  required, because mzLib has no default. Rows add `measured`, `counted` and
  `residual_variance_given`, and `"separated"` is a status, never a p-value. The tests hold the
  recordings to lme4's `glmer()` and to `glmmTMB`, whose fits of mzLib's reference input ship in
  `extdata/stats/`. A function of its own, so `stats_fit_mixed()`'s answer does not change.
* **`stats_fit_link()`: fit that detection link** (wire `stats fit-link`, mzLib #1469) from cells
  with a latent mean and SD and a count, on detection (any success or none) by default;
  `with_depth`, `with_overdispersion` and `detection_only` are sent only when given, and the result
  reports what mzLib used. The link is an `mzlibr_detection_link`, `NULL` when mzLib's search did
  not converge.
* **`ptmqtl_pairs(include_bounds = TRUE)`** lets a run where a site was seen in only one form into
  co-varying pairs as a tied bound, never an invented value: a covered-zero or floor ranks below
  every fraction, a ceiling above (mzLib #1468). A site often below detection can then be tested at
  all; a pair whose two sites are both mostly one bound is `"both_sites_mostly_bounds"`, not tested.
* **`ptmqtl_pairs(run_covariate = )`** correlates two sites given each run's depth (or any per-run
  value), a partial Spearman correlation (mzLib #1468). The covariate's name, file and the runs it
  leaves out are in the result.
* **`ptmqtl_trait_effects(scale = "log2_odds")`** reports effects in log2 odds, where +1 doubles the
  odds (mzLib #1464): mzLib divides the natural-log fit by ln 2, so `df`, `p_value` and `q` do not
  change. Left out, mzLib's default stays natural-log odds.

The articles teach each (`vignette("proteins")`, `vignette("library")`, `vignette("stats")`,
`vignette("ptmqtl")`), and their citation blocks now carry the new verbs' references. The parity
check lists no pyMzLib 0.7.0 function or option as not yet ported.

*Calls that used to fail now succeed.*

* **`sdrf_isobaric_design()` on a TMT10 SDRF that names its last channel `TMT131`** (mzLib #1467).
  PRIDE spells TMT10's last channel `TMT131`, and it was refused as not a channel of TMT10. It now
  reads as channel `131N`. `TMT131` is still refused for TMT11, where `131N` and `131C` are both
  channels.
* **`sdrf_design()` on an SDRF with a missing replicate or a numbering gap** (mzLib #1426). A
  biological replicate of `not available`, or no replicate column at all, was refused; mzLib now
  adds a number after the highest one the condition uses, and says so in `notes`. A gap in the
  fraction or technical replicate numbers is noted instead of refused.
* **A `peptidoform_convert()` batch with two modifications on one residue** (mzLib #1457). One such
  line, for example `PEPM[A][B]TIDE`, failed the whole call with an `mzlib_bridge_error`
  (`ArgumentException`). Now only that line fails (below), and the rest of the batch converts.

*Calls that used to succeed are now refused.*

* **Malformed modification text in `peptidoform_convert()`** (mzLib #1457), with `source =
  "mzLib"` (the default) or `"MassShift"`. An empty modification (`PEP[]TIDE`), a modification that
  does not follow an unmodified residue (two leading groups, or two on one residue) and anything
  after the C-terminal modification were converted with a modification lost or moved. Each is now a
  `"failed"` row with `failure_reason` `"UnknownFormat"`. An ambiguous MetaMorpheus sequence
  (`|`-joined candidates) still converts with a warning, unless a candidate ends in a C-terminal
  modification or the next one's N-terminal modification lands on an already modified residue.

*Same failure, different error.*

* **A `.parquet` another process holds, or that you may not read** (mzLib #1478). Every reader
  raised an `mzlib_usage_error` (`Parquet file type not supported`). It now raises an
  `mzlib_bridge_error` naming the I/O failure (`IOException` or `UnauthorizedAccessException`). A
  file that is not Parquet at all is still the usage error.

*Same file, different value.*

* **`sdrf_design()` keeps the SDRF's biological replicate numbers** (mzLib #1426). It used to rank
  them within each condition, so an SDRF that numbers a condition's samples 22, 23, 24 gave 0, 1, 2.
  It now gives 21, 22, 23 (0-based), gaps included, because replicate 4 of one condition can be the
  same subject as replicate 4 of another, and a gap can be a lost sample. `notes` flags each
  condition whose numbers are not 1..N, where it used to list the renumbering, and the result
  carries a seventh caveat saying so. A design written with `out` carries those numbers too;
  MetaMorpheus refuses one with gaps until a release carries smith-chem-wisc/MetaMorpheus#2928.
* **`diff_design()` numbers contrast ids once over the whole analysis** (mzLib #1476). With
  `strata`, `c1` used to be each stratum's own first comparison; it is now one comparison in every
  stratum, and a stratum that cannot run it skips the id. A stratum's `not_run` lists the skipped
  contrast under its own label (`condition=B vs condition=A`, where it said
  `condition (one level: A)`) and gains an `id` column, `NA` when no default contrast stands for
  the entry. Two caveats are reworded to say so. Designs without strata are unchanged.
* **`ptmqtl_occupancy()` reports more `"not_detected"` cells** (mzLib #1468). mzLib now defines
  `not_detected` as a covered-zero: peptidoforms covering the site were identified in the run and
  none carried the modification. The covering forms no longer have to be quantified, so a site
  covered only by unquantified forms gets a row it did not get before. `fraction` stays `NA`.
* **`ptmqtl_occupancy()` gains two columns, and `out` writes them**: `covering_peptidoforms` and
  `modified_peptidoforms`, after `unmodified_quantified` (mzLib #1468), the k of n of each cell,
  counting identifications rather than intensity. Read a written table by header name, not by
  position; `ptmqtl_trait_effects()` reads both the old and the new layout. A new caveat says what
  the counts are. The shipped `extdata/ptmqtl/trait_occupancy.tsv` is the new layout.
* **`ptmqtl_pairs(type = "covarying")` reports `exclude_ceiling` `NA`** when you do not pass
  `keep_ceiling = TRUE`, where it reported `TRUE`. mzLib made the option nullable and decides it
  inside `PtmPairEngine.CoVarying` (mzLib #1468); mzLibR does not restate mzLib's rule to name the
  value. The pairs themselves are unchanged: ceilings are still left out by default.
* **`ptmqtl_pairs()` says why a co-varying pair was not tested**, in a new `test_status` column
  (`"tested"`, `"not_estimable"`, `"too_few_shared_runs"`, `"both_sites_mostly_bounds"`; `NA` for
  physical pairs), beside `bounds_included` and `run_covariate_removed`. The result gains
  `include_bounds`, `min_untied_share`, `run_covariate_file`, `run_covariate_name` and
  `run_covariate_runs_missing`, mzLib's echo of the two new options (above). Its caveat now
  counts untested pairs by reason ("N co-varying pairs were not tested (...)"), where it said "N
  co-varying pairs are not_estimable". `spearman_method` is still `"not_estimable"` for every
  untested pair.
* **`ptmqtl_trait_effects()` reports its `scale`** (mzLib #1464), `"natural_log_odds"` at mzLib's
  default, the unit of every effect, standard error and variance; a seventh caveat names it. The
  numbers are unchanged unless you pass `scale = "log2_odds"` (above).
* **`isobaric_kits()` gains an `sdrf_label` column** (mzLib #1467): each channel's SDRF
  `comment[label]` by mzLib's `SdrfIsobaricLabels`, such as `TMT127N`, `ITRAQ114`, `DiLeu115a`, and
  `TMT131` for TMT10's 131N. A sixth caveat says what the label is.
* **ProForma and Unimod output for a terminal modification that shares its name** (mzLib #1442).
  Entries such as `Methyl on X` and `Less Common:Methylation on X` exist once for a side chain and
  once for a terminus. The terminal one now resolves, so `peptidoform_convert(..., target =
  "ProForma")` and `target = "Unimod"`, and `readers_read_records()`'s `pro_forma` column, write its
  `UNIMOD:` accession where they wrote its typed name or left it unresolved.
* **`peptidoform_fragments()` gives a peptide with a permanently charged modification a mass one
  proton lower per charge.** The bridge now takes UniProt's modifications from mzLib's own
  corrected list (`Mods.UniprotModifications`), as MetaMorpheus does. Trimethyllysine counts as
  C3H6, 42.046951 Da (Unimod's Trimethyl), not UniProt's charged 43.054227; the other modifications
  PSI-MOD marks charged move the same way (trimethylalanine, dimethylproline,
  O-(2-cholinephosphoryl)serine). This changes `monoisotopic_mass`, the `modifications` table's
  `mass` and every fragment's `neutral_mass` that includes it. `peptide_mz()` is unchanged, because
  `fixed_charges` drops from 1 to 0 at the same time. N,N,N-trimethylglycine is the exception: it
  stays at UniProt's mass, one hydrogen heavy, until the bridge builds from an mzLib that includes
  mzLib #1415.
* **`fixed_charges` (in `peptides`) and `formal_charge` (in `modifications`) are deprecated** and
  always 0. No peptide carries a charge on top of its mass any more; both columns stay until the
  bridge's next protocol version removes them. `?peptide_mz` and `?peptidoform_fragments` say so,
  and a fragment's m/z at charge z is now simply `(neutral_mass + z * 1.00727646677) / z`.

## Statistics, designs and PTM-QTL (the pyMzLib 0.6.0 bridge)

Ported from pyMzLib 0.6.0, verb for verb; each *needs the pyMzLib 0.6.0 bridge*. Every number is
mzLib's, and every recording the help pages and articles replay is pyMzLib's, byte for byte. New
articles: `vignette("diff")` and `vignette("ptmqtl")`; the differential-abundance article gains a
section for each new statistic.

* **`stats_fit(contrasts = )`: any comparison, with a confidence interval** (mzLib #1452,
  `EmpiricalBayes.ModerateContrasts`). Give the contrasts file `diff_design(out = )` writes beside
  the design (`contrasts_<i>.tsv`), or your own: a contrast id column, then exactly the design's
  coefficient columns, one weight each. Each contrast is moderated with the same prior as the
  coefficients, has its own Benjamini-Hochberg family, and carries `confidence_low` and
  `confidence_high`, as limma's `contrasts.fit()`, `eBayes()` and `topTable(confint = )`. The answers
  are a second table, `fit$contrasts` (an `mzlibr_moderated_contrasts` with `confidence_level`,
  `tested_contrasts`, `weights` and `records`), so the coefficient table is unchanged;
  `coefficients` may now be empty and defaults to `character(0)`. `confidence` sets the level, as
  a fraction; left `NULL`, it is mzLib's default (0.95). `fit$contrasts_file` names the file.
* **`stats_fit(estimator = )`: mzLib's marginal-likelihood variance prior** (mzLib #1368).
  `"marginal_likelihood"` fits the prior with each feature's own residual df, so the unequal df that
  missing values cause are handled exactly; `"moments_legacy"` (limma's `eBayes(legacy = TRUE)`)
  stays mzLib's default, and the new `fit$prior$estimator` reports which was used. On the MALAT1
  dilution series it lowers the prior df from 19.6 to 11.0 and drops the false FLuc call at half
  load.
* **An older bridge is named, not trusted.** A bridge before pyMzLib 0.6.0 lists `stats fit` but
  ignores `--contrasts` and `--estimator`; `stats_fit()` then raises an `mzlib_usage_error` naming
  the 0.6.0 bridge instead of returning a fit that did not do what was asked.
* **`stats_fit()`'s caveats name the variance-prior estimator.** The first caveat now begins
  `moments_legacy:` (or `marginal_likelihood:`), the `residual_df_differ` caveat names
  `estimator = "marginal_likelihood"`, and the `bh_adjusted` caveat says "per coefficient and per
  contrast". Every number is unchanged; only code that matches caveat text needs a look. On fewer
  than 3 features `prior$scale` is `NA` and `prior$df` 0 (mzLib #1398); `?stats_fit` now says
  bridges before pyMzLib 0.6.0, not 0.5.1, reported the first feature's variance there.
* **`stats_combine_p()`**: one p-value per feature from several independent studies, by Fisher's
  method or Stouffer's, weighted or not (mzLib #1367, `PValueCombination`). `method` is required.
  Give each test's `effect` and every p-value is read as two-sided and converted to the one-sided
  p-value for effect > 0 first, so opposite effects cancel. The tests hold the recordings to R's
  `pchisq()` and `qnorm()`, recorded and computed here.
* **`stats_correlate()`**: Spearman's rho of every feature with a per-sample trait, optionally
  partial given covariates (`given`; mzLib #1367, #1446, `SpearmanCorrelation`). Each row carries
  rho, a two-sided p-value, the samples used and how the p-value was computed (`"exact"`,
  `"asymptotic"` or `"not_estimable"`, with rho and p `NA`). The tests hold every recorded row to
  `stats::cor.test()` run here.
* **`stats_fit_logistic()`**: a logistic regression per feature of a 0/1 table, as R's
  `glm(family = binomial)` (mzLib #1367): is a feature's detection related to the condition? A
  feature the design separates is `"separated"`, with `NA` statistics, instead of R's huge
  coefficient. `max_iterations` and `tolerance` are mzLib's when left `NULL`. The tests hold the
  recording to R 4.6.1's `glm()` output and to `glm.fit()` run here.
* **`stats_fit_mixed()`**: a linear mixed model per feature with a random intercept per group of
  samples, as `nlme::lme(random = ~ 1 | group)` (mzLib #1367), tested by Wald t on nlme's
  containment df. `estimator = "reml"` or `"ml"`; left `NULL`, mzLib's default (REML), which the
  result reports. The tests hold the recordings to nlme 3.1.169's output and to `nlme::lme()` run
  here.
* **`stats_entrapment_fdp()`**: the lower bound, combined and paired estimates of a search's false
  discovery proportion at each q-value threshold (mzLib #1385, Wen et al. 2025, as FDRBench).
  `ratio` is required; `paired` is reported only at a ratio of exactly 1; `thresholds` left `NULL`
  is mzLib's default sweep, 0.001 to 0.100.
* **`diff_design()`, with `factor_spec()` and `covariate_spec()`**: a sample sheet to mzLib's design
  matrix per stratum (mzLib #1450, `AnalysisDesign`), equal to R's `model.matrix()` (the tests run
  it here), with its rank, aliased terms, warnings and the contrasts it cannot run; one row per
  default contrast with **mzLib's verdict on whether these samples can estimate it at all**. Each
  stratum's `matrix` is an R matrix named by sample id, and `strata` is a list named by stratum key.
  With `out` it writes `design_<i>.tsv` and `contrasts_<i>.tsv`, which `stats_fit()` reads.
* **`diff_observations()`**: every peptide in every biological sample from a FlashLFQ or
  MetaMorpheus peptide table (mzLib #1455): log2 intensity, fractions summed and technical
  replicates averaged, or **`NA`, never 0, with the reason in `state`**. `runs` is a data.frame of
  `file` and `sample`, and optionally `fraction` and `technical_replicate`; `basis` has no default.
  mzLib's shared-peptide normalization is always measured and reported, and subtracted with
  `normalize = TRUE`. `responses_out` writes the table `stats_fit()` reads.
* **`ptmqtl_occupancy()`, `ptmqtl_pairs()` and `ptmqtl_trait_effects()`**: PTM site occupancy per
  run with each cell's state (a censored cell has no `fraction`, never a 0), physical and
  co-varying PTM pairs (`type` is required; `remove_run_level` takes out a batch shift), and a
  per-run trait's effect on logit occupancy with a random intercept per biological replicate
  (mzLib #1369, #1430, #1446). **Positions are 1-based, unlike `readers_read_occupancy()`**.
* The vendored specs are re-synced from bridge `60b1d259`. `scripts/sync-specs.R` now reads a bare
  `y` or `n` in a spec as the string the bridge's own YAML reader takes it for, where R's `yaml`
  package read it as a boolean (`stats correlate`'s and `ptmqtl pairs`' `n` column).

## Data: checksums, alignment, orthology, spectral libraries and SDRF scoping (the pyMzLib 0.6.0 bridge)

Ported from pyMzLib 0.6.0, verb for verb; each *needs the pyMzLib 0.6.0 bridge* and is refused,
naming that release, by an older one. Every option left out is not sent, so mzLib's default
applies and the result says what was used. Every recording the help pages and articles replay is
pyMzLib's, byte for byte. A new article, `vignette("library")`; the PRIDE, protein-database and
SDRF articles gain a section for each new verb.

* **`pride_file_checksums()`: the MD5 and exact size PRIDE recorded for each deposited file**
  (mzLib #1427, wire `pride checksums`). A data.frame of `file_name`, `md5`, `size_bytes` and
  `project_accession`, its row names the file names, so `checksums["README.txt", ]` looks a file
  up by the name it was downloaded under. `size_bytes` is the exact download length, which neither
  other listing gives (PXD000001's `.xml.gz`: 10,668,000 bytes, where `pride_list_files()` says
  10,677,205), and `md5` is not the manifest's `checksum` (a SHA-1 where present). It is a third
  list: PRIDE's `generated/` files have no row. Zero rows is an answer; an unknown accession, or a
  project not yet published, raises `mzlib_project_not_found`; a list mzLib cannot read stays an
  `mzlib_bridge_error` (`MzLibException`). No timeout of mzLibR's own. `pride_download()` and
  `pride_download_files()` still check no bytes, and now say so; the PRIDE article shows the check
  with `file.size()` and `tools::md5sum()`.
* **`proteins_align()`: map a residue onto an isoform, an ortholog or another release** (mzLib
  #1428, wire `proteins align`). A global affine alignment scored with BLOSUM62, for a list of
  `c(source, target)` pairs in one bridge call: `records` (one row per pair: score in half-bits,
  identical and substituted residues, residues facing a gap, the aligned rows) and `residue_map`
  (one row per residue, its counterpart's position and letter, `NA` facing a gap, never the
  nearest residue). `pair_index` is 1-based. `proteins_map_source_to_target()` and
  `proteins_map_target_to_source()` look positions up, and `proteins_alignments()` splits a result
  into one alignment per pair. `gap_open`, `gap_extend`, `end_gaps` and `max_cells` default to
  `NULL`, mzLib's (11, 1, `"free"`, 200,000,000); `residue_map = FALSE` sends `--no-residue-map`.
  One pair is `proteins_align(list(c(source, target)))`, as the spec spells it: there is no
  separate one-pair function, which is pyMzLib's `proteins.align()`.
* **`proteins_parse_accessions()`: entry, isoform or version, and applied variants** (mzLib #1382,
  wire `proteins parse-accession`). `P04637-2_S70N` is entry `P04637`, isoform 2, with `S70N`;
  `NP_000537.3_R72P` is RefSeq `NP_000537`, version 3. It parses and never repairs: a FASTA header,
  a `DECOY_` prefix, a `|`-joined group or mzLib's `P04637_2` collision counter is `unrecognized`
  and kept verbatim, and every accession is sent exactly as given, not trimmed.
* **`proteins_build_orthology()`: an Ensembl Compara orthology snapshot** (mzLib #1381, wire
  `genes orthology`). From each species' GTF (or gene table), every species' homology dump - mzLib
  refuses a call without one - and the gene-tree file, all local; `gtfs` and `gene_sets` are named
  character vectors, production name to path. It reports what went in (`gene_sets`, `gene_trees`,
  `dumps`), the relationships per species pair (`pairs`; `proteins_orthology_pair()` finds one in
  either order), a `snapshot_id` to record with any result, with `pair_status` every gene's status
  in one other species (`has_ortholog`, `no_edge_in_shared_tree`, `tree_lacks_target_species`,
  `not_in_any_tree`), and with `out` the Parquet snapshot it wrote (`written`). mzLibR reads no
  Parquet; the proteins article queries it with `arrow`. mzLib's three-species test world is
  copied into `tests/orthology/`, kept byte-exact by `.gitattributes` because the
  recordings' sha256s are over its bytes.
* **`library_read()` and `library_read_many()`: mzLib's `.msl` spectral libraries** (wire `library
  read`). One row per precursor, or per fragment ion with `fragments = TRUE`, in the library's
  stored order; `retention_time_unit` is per row (iRT or minutes, never converted), and mzLib's
  "not available" values (ion mobility 0.0, a NaN q-value, NCE 0, an empty protein) are `NA`.
  `precursor_index` is 1-based. `index_only_threshold` defaults to `NULL`, mzLib's 1 GiB;
  `limit` and `offset` count precursors. The many-libraries form returns one long table and a
  `files` data.frame, in one bridge call.
* **`sdrf_isobaric_design()`: a TMT, TMTpro, iTRAQ or DiLeu design from an SDRF** (wire `sdrf
  isobaric-design`). One row per file and annotated channel, or every reason there is none. It
  guesses neither the kit nor the plex: declare exactly one of `plex_column`,
  `plex_file_name_pattern` and `single_plex`. **Its numbers are 1-based, as the SDRF writes
  them**, where `sdrf_design()`'s are 0-based; the split is mzLib's, and mzLibR converts neither.
  `out` names a folder for MetaMorpheus's `TmtDesign.txt`.
* **`sdrf_restrict()`: an SDRF cut to the files one search read** (wire `sdrf restrict`). mzLib
  joins `X-calib.mzML` to the row naming `X.raw` by stem, adds `comment[searched data file]`, and
  reports what it could not place (`searched_without_row`, `dropped_data_files`, `ambiguous`,
  `shared_stem`). Name an element of `searched_files` by its acquired file to override the
  suffix rule. The result is an `mzlibr_sdrf`; `sdrf_searched_data_files()` reads the new column.
  Restrict first, then read a design: both design readers key on that column.

*Same file, different value.*

* **`sdrf_pool(out = )`'s `written$path` names the file mzLib wrote.** mzLib adds `.sdrf.tsv` to a
  name that does not already end in it, so `out = "x.tsv"` writes `x.tsv.sdrf.tsv`. With the
  pyMzLib 0.6.0 bridge `written$path` says `x.tsv.sdrf.tsv`; older bridges reported `x.tsv`, a
  file that does not exist. The file written is the same, and an `out` already ending in
  `.sdrf.tsv` is unchanged. `sdrf_restrict(out = )` reports its file the same way.

*For contributors.* The replay bridge (`R/replay.R`) gains pyMzLib 0.6.0's rules for these verbs:
a recording answers only a call with the same pairs, accessions, searched files or species on
stdin, the same `--sequences`, `--no-residue-map`, `--pair-status` and `--restrict-gene-trees`
choice, and the options it was recorded with; and, as in pyMzLib, a recording whose echo of an
option is null never answers a call that sends it. The test harness gains
`recorded_bridge_call()`, which runs a public function through the real transport and verb check
against a recording and returns what it sent.

## Defaults come from mzLib (the pyMzLib 0.5.0 bridge)

* **`peptidoform_fragments()` takes its defaults from mzLib, so a call that names no protease or
  dissociation type gives different peptides and ions.** Every digestion argument (`protease`,
  `dissociation`, `missed_cleavages`, `min_length`, `max_length`, `max_modifications`,
  `max_isoforms`, `terminus`) now defaults to `NULL`, mzLibR sends nothing for it, and mzLib's
  `DigestionParams` default applies. Two of the old defaults were mzLibR's own and change:
  - protease: mzLib's default is `"trypsin"`, which also cuts before proline. The old default was
    `"trypsin|P"`, which does not. Albumin's default digest goes from 303 to 310 peptidoforms, and
    from 195 to 202 distinct sequences.
  - dissociation: mzLib names no default, so the bridge uses HCD (b and y ions) until it does
    (smith-chem-wisc/mzLib#1488). The old default was ETD (c and z-dot ions).

  The other defaults (missed cleavages, the length window, the modification and isoform caps,
  terminus) have the same values as before, now set by mzLib. **To get the old answer, pass
  `protease = "trypsin|P", dissociation = "ETD"`.** A digest records `min_length`, `max_length`
  and `missed_cleavages` as `NA` when they were left to mzLib, since the bridge does not report
  those three back, and prints them as mzLib's default. A `max_length` below `min_length` is
  refused only when both are given. The help page's counts are re-measured at mzLib `d10017ce`.
* **`pride_list_files()` and `pride_search()` leave the page size to mzLib.** `page_size` now
  defaults to `NULL` and is sent only when given, so mzLib's own paging default applies (100
  today, as before). The results never depended on it. A given page size is checked as before.
* **`pride_list_files()`, `pride_list_ftp_files()` and `pride_search()` add no timeout of their
  own.** They defaulted to `timeout = 300`, which turned a PRIDE outage that hangs into an
  `mzlib_timeout`. mzLib already bounds each request (#1412: three retries after 5, 20 and 60 s,
  each attempt up to 100 s), so they now default to `NULL` and a hang is reported as
  `mzlib_service_unavailable`, after about 8 minutes per request. Pass `timeout = 300` for the old
  limit; a number still cuts the call off with `mzlib_timeout`. The PRIDE article's outage
  handling now waits minutes and tries once more, since mzLib has already retried, and
  `?mzlib_error` says so.
* `?stats_fit` no longer says the bridge's spline default is 4. With `trend = TRUE` and no
  `spline_basis`, mzLibR has always sent none; the 0.5.0 bridge now asks mzLib for none either,
  so mzLib sizes the trend as limma does and a trend fit on fewer than 30 features gets limma's
  numbers. `prior$spline_basis_count` says what was used.
* The help pages are re-rendered from the bridge's specs, re-checked against mzLib `d10017ce`.
  PXD000001's REST manifest listed 9 of the 14 files in its FTP tree on 2026-10-21, and
  `?pride_list_files` says so rather than 8 of 13.

## mzLib 1.0.595

Built from mzLib 1.0.595 (tag commit `d10017ce`), which carries mzLib 1.0.594 as well.

* **DIA-NN 2.x's `report.parquet` reads** (mzLib #1384, file type `DiaNnReportParquet`). DIA-NN
  2.x writes its main report only as parquet. mzLib reads it into the same record as the 1.8/1.9
  `report.tsv`, so it is a fifth quantifiable type: `readers_read_results()`,
  `readers_read_records()` and `flashlfq_quantify()` take it. **`readers_formats()` lists 39
  types, up from 38, and 5 offer the `quantifiable` view, up from 4**; 19 still offer none. Like
  the TSV report it is recognised by its columns (`Precursor.Id`, `Stripped.Sequence`, `Run`), not
  its name, and any other `.parquet` is an `mzlib_usage_error`. **Its `retention_time_unit` is
  `"minutes"`**, as for the TSV report. `absent_fields` is always empty for it, because there is no
  text header to judge from, so the columns 2.x stopped writing (`File.Name`, `PG.Quantity`,
  `MS2.Scan` and others) are `NA` in `readers_read_records()` without being named there, and
  `ms2_scan_number` reads 0.
* **`peptidoform_convert(target = "ProForma")` names UniProt modifications by UNIMOD accession**
  (mzLib 1.0.594, #1403). `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]-SEQK`, where it
  was written back by name; the caveat that said so is retired, so a result carries five caveats,
  not six. The `pro_forma` column `readers_read_records()` computes for a `.psmtsv` written
  without one gains the same accessions.

Every other result mzLib 1.0.594 and 1.0.595 changed, listed after the fact: pyMzLib measured
each with the bridge at `d10017ce`, and the conversions, the `stats_fit()` prior and the
`entrapment_fraction` column were re-checked from R against the pyMzLib 0.5.0 bridge. The help
pages carry them too, rendered from the bridge's specs.

*Calls that used to succeed are now refused.*

* **ProForma input that lost part of itself** (mzLib #1462).
  `peptidoform_convert(..., source = "ProForma")` converted `PEPTIDE-[UNIMOD:2][UNIMOD:35]` with
  one C-terminal modification dropped, and read `PEPTIDE-[UNIMOD:2]K` as `PEPTIDEK`. Stacked
  terminal modifications now fail with `failure_reason` `"IncompatibleModifications"`
  (`incompatible_items`: "Multiple C-terminal modifications (2 found)."), and a residue after a
  C-terminal modification, a leading `-` or a modification with no sequence fail with
  `"InvalidSequence"`.
* **SUMO remnants named `... on D`** (mzLib #1432). mzLib's human SUMO-1 and SUMO-2/3 remnants
  (`DVIEVYQEQTGG (SUMO-1 Site human)`, `DVFQQQTGG (SUMO-2/3 Site human)`) now target K, where
  they targeted D. A full sequence that names the `on D` form no longer converts; the `on K` form
  does.
* **TMTpro on S, T, Y or H** (mzLib #1419). mzLib used to turn these into TMT, 224 Da, so
  `PEPS[Common Fixed:TMTpro on S]K` converted at the wrong mass. They now fail to convert. TMTpro
  on K is unchanged (+304.2071 Da).

*Calls that used to fail now succeed.*

* **`flashlfq_quantify()` with gapped replicate numbers** (mzLib #1422). A design whose
  biological replicates skip a number, with fractions, raised an error; it now quantifies. So does
  the numbering the bridge uses when `spectra` gives no `biological_replicate`, which skips at
  every condition after the first. `bayesian_protein_quant = TRUE` on such a design now writes
  `BayesianProteinQuant.tsv` to `output_directory`.
* **A condition whose name ends in digits** (mzLib #1424). Condition `A1` biorep 0 beside
  condition `A` biorep 10 raised an `mzlib_bridge_error` in `flashlfq_quantify()` and
  `flashlfq_median_polish()`, because the two spell the same text. It now quantifies.

*Same file, different value.*

* **Entrapment and contaminant flags read every joined label** (mzLib #1371). In
  `readers_read_records()` and `readers_read_protein_groups()`, `is_decoy`, `is_contaminant` and
  `is_entrapment` are `TRUE` when any `|`-joined parent is D, C or E. A protein group labelled
  `T|C` is now a contaminant and `T|ET` entrapment; on a `.psmtsv` the "E" in MetaMorpheus's
  "Output too long for Excel" no longer reads as entrapment. Labels MetaMorpheus writes itself
  read as before.
* **`readers_read_records()` of a `.psmtsv` or `.osmtsv` gains an `entrapment_fraction` column**
  (mzLib #1371): the share of each PSM that is an entrapment discovery, from 0 to 1. A shared
  `T|ET` counts half when the sequences differ and nothing when both proteins carry the same
  peptide. Use it, not `is_entrapment`, to count entrapment discoveries. It is `NA` where mzLib
  cannot line the label up with the full sequence (mzLib's `NaN`).
* **Lactylation on K weighs one hydrogen more** (mzLib #1417). mzLib's formula is now C3H4O2,
  72.021129 Da, where it was C3H3O2, 71.013304 Da. `peptidoform_convert(..., target =
  "MassShift")` of `PEPK[Less Common:Lactylation on K]R` gives `PEPK[+72.021129]R`, and the
  `"ProForma"` and `"Unimod"` targets write `UNIMOD:2114`. A `.psmtsv`'s `pro_forma` column changes
  the same way.
* **`flashlfq_quantify(normalize = TRUE)`** (mzLib #1422). The reference is now the first
  condition's lowest biological replicate, and each fraction's lowest technical replicate, taken
  in numeric order. A design numbered from 1 is now normalized; before, biorep normalization did
  nothing for it. A gapless design numbered from 0 and listed in order is unchanged.
* **`stats_fit()` on fewer than 3 features** (mzLib 1.0.594, #1398). mzLib fits no prior, as
  limma does: `prior$df` is 0, each feature keeps its own variance, and t and p are ordinary
  per-feature t-tests. They used to borrow a prior fitted to the two features. The 0.5.0 bridge
  still reports the first feature's own variance as `prior$scale` there; bridges from pyMzLib
  0.5.1 report `NA`, and each feature's value is in the `prior_variance` column. `?stats_fit` says
  both. Tables of 3 or more features are unchanged.
* **`peptidoform_fragments()` with a `|P` protease** (mzLib #1320). `protease = "trypsin|P"` cut a
  protein that begins `RP` or `KP` after its first residue; it no longer does, and the other `|P`
  proteases change the same way. mzLib's default `"trypsin"`, and any protein that begins with M,
  are unaffected.
* **PRIDE: a failure that can pass is retried** (mzLib 1.0.594, #1412). Every `pride_` function
  retries such a failure three times, after 5, 20 and 60 seconds, so an outage takes longer to
  report. An FTP listing that comes back empty, which `pride_list_ftp_files()` used to return as
  no files (and so `mzlib_project_not_found`), is retried and then raised as an
  `mzlib_bridge_error` with `error_type` `"HttpRequestException"`. A failed download leaves
  `<name>.partial` and `<name>.partial.validator` behind, so the next call resumes.
* **PRIDE downloads are not checked yet.** mzLib can check a download against PRIDE's checksums
  and catch an EBI error page served at the right size (#1427, #1449), but not on the path
  `pride_download()` and `pride_download_files()` use, so a file can be reported as downloaded
  when its bytes are wrong. `?pride_download` and the PRIDE article say to check large files
  yourself.
* `?stats_fit` no longer says contrasts are not implemented: mzLib has them
  (`EmpiricalBayes.ModerateContrasts`, #1452), but `stats_fit()` does not expose them yet.
* The psmtsv caveat on `readers_read_results()` cites `SpectrumMatchFromTsv.cs:122` and `:197`;
  mzLib #1371 moved the lines it cited. The claim is unchanged.

## mzLib 1.0.593

* `readers_read_protein_groups()`, `readers_read_quantified_peptides()` and
  `readers_read_occupancy()` read an RNA search's `AllQuantifiedTranscriptGroups.tsv` and
  `AllQuantifiedOligos.tsv` (mzLib #1388). mzLib reads them with subclasses of the protein-group
  and peptide readers, so the columns keep their protein names: `protein_group_name` is the
  transcript group, `sequence` the oligonucleotide, and occupancy names RNA modifications.
* MetaMorpheus protein-group tables written without quantification now read (mzLib #1365):
  `AllProteinGroups.tsv` and each file's `<file>_ProteinGroups.tsv` failed with `Tsv file type not
  supported`, and now read as `MetaMorpheusQuantifiedProteinGroups`, with `intensity` in
  `absent_fields`.
* Four shipped modifications now write their Unimod accession in `pro_forma` (mzLib #1328):
  `GG (Ubiquitination Site)` as `[UNIMOD:121]`, both `Myristoylation` entries as `[UNIMOD:45]` and
  `EQIGG` as `[UNIMOD:846]`, instead of by name. Masses are unchanged.
* A read fault on an existing `.mzid` is an `mzlib_bridge_error` naming the file (mzLib #1362),
  where it was an `IOException`. A missing file is still an `mzlib_usage_error`.
* PRIDE download errors name the file and the host, never the URL (mzLib #1350), so a reviewer
  token in a query string cannot reach a log. They are still `mzlib_service_unavailable`.
* The help pages no longer state how many file types mzLib recognises: `readers_formats()` asks
  mzLib, and a number written down goes stale at the next release.

## Reading files

* `readers_formats()` and `readers_identify()` say what a file is and which cross-format views
  it offers.
* `readers_read_spectra()` reads scans from mzML, Thermo `.raw`, Bruker `.d`, timsTOF `.d`, MGF
  and msalign: headers always, peaks on request.
* `readers_read_results()`, `readers_read_features()` and `readers_read_matches()` read the three
  cross-format views of search results; `readers_read_records()` reads any recognised file into
  that format's own fields.
* `readers_retention_time_in_minutes()` converts using each file's declared unit, and refuses
  when there is none - for a `_many` batch, each file's rows by that file's unit.
* Many files in one call: `readers_identify_many()` and a `_many` form of every reader hand the
  whole list to one bridge process with `threads` and `on_error` stated, and return one long
  `records` table plus a `files` data.frame of per-file facts, one row per input
  (*needs the pyMzLib 0.2.0 bridge*).
* Every read reports the per-file block: `reader`, `rows_not_read`, `retention_time_unit`,
  `caveats`, and the three ways a column can be empty kept apart - `absent_fields` (the format has
  no source for it, so it is `NA` in every row), `failed_fields` and `excluded_fields`.
* `readers_read_spectra()` reports the run's `source`: instrument model and its PSI-MS accession,
  serial number, and acquisition start time with whether it is UTC (mzLib #1349).
* `readers_read_matches()` reads mzIdentML (`.mzid`, `.mzid.gz`) with each match's `q_value`,
  `rank` and `pass_threshold`, and with `scores = TRUE` one row per engine score; mzIdentML items
  mzLib does not represent are listed in `skipped` with the reason.
* New `readers_read_protein_groups()`, `readers_read_quantified_peptides()` and
  `readers_read_occupancy()` read MetaMorpheus and FlashLFQ quantification tables long, one row per
  record per sample, and per-site PTM occupancy (mzLib #1347; *needs the pyMzLib 0.2.0 bridge*).
* `mzlibr_bridge_version()` reports `verbs`, every command the bridge dispatches.
* mzLib now recognises Pytheas, mzIdentML and its gzipped form, and the MetaMorpheus and
  FlashLFQ quantification tables.
* `limit = 0` is accepted by every reader, as the bridge allows: the envelope alone.

## SDRF-Proteomics

* `sdrf_read()` reads one experimental-design file and `sdrf_pool()` pools several into one
  analysis table with provenance; `sdrf_value()`, `sdrf_all()`, `sdrf_records()` and the other
  `sdrf_` accessors read them.
* `sdrf_validate()`, `sdrf_lint()`, `sdrf_assess()`, `sdrf_samples()` and `sdrf_parse_ages()` ask
  mzLib whether an SDRF document is well-formed, whether several agree on their vocabulary,
  whether it describes its samples, and what its ages are in years; `sdrf_validate_many()`,
  `sdrf_assess_many()` and `sdrf_samples_many()` read a corpus in one bridge call
  (mzLib #1207, #1325, #1326, #1333, #1335).
* New `sdrf_design()` reads the label-free experimental design MetaMorpheus and FlashLFQ take out
  of an SDRF, with mzLib's `SdrfLabelFreeDesign`, or every reason it refuses to - all at once,
  and as a result rather than an error. `sdrf_design_spectra()` and `sdrf_design_run_design()`
  hand it to `flashlfq_quantify()` and `flashlfq_median_polish()`; `out` writes MetaMorpheus's
  `ExperimentalDesign.tsv`. Its replicate and fraction coordinates stay 0-based, as the quant
  functions take them (mzLib #1363; *needs the pyMzLib 0.3.0 bridge*).

## Isobaric kits

* New `isobaric_kits()` lists mzLib's TMT, TMTpro, iTRAQ and DiLeu kits with every channel's
  label, theoretical reporter-ion m/z and matching window, or one kit by MetaMorpheus's name for
  it; `channel_index` is 1-based (mzLib #1375; *needs the pyMzLib 0.3.0 bridge*).

## PRIDE, peptidoforms and quantification

* `pride_list_files()`, `pride_list_ftp_files()`, `pride_download()` and
  `pride_download_files()` list and fetch PRIDE Archive projects.
* `pride_search()` finds PRIDE projects by keyword, one row per project, with dates as `Date`.
* `peptidoform_fragments()` digests a UniProt entry and fragments its peptidoforms.
* New `peptidoform_convert()` rewrites full sequences in another notation with mzLib's
  `SequenceConversionService`. The main use is MetaMorpheus full sequences to Unimod accessions:
  `[UniProt:N-acetylserine on S]SEQK` becomes `[UNIMOD:1]SEQK`, and dimethyllysine becomes
  `UNIMOD:36`. One row per input, in order, in `records`, with mzLib's own verdict (`converted`,
  `converted_with_warnings` or `failed`), the modifications it could not write, and its warnings.
  `source` and `target` take any notation mzLib has registered; every result lists them. `mode` is
  mzLib's `SequenceConversionHandlingMode`. The ProForma target names UniProt modifications by
  UNIMOD accession too, since mzLib #1403 (see mzLib 1.0.595 above; at mzLib 1.0.593 it wrote them
  back by name, mzLib#1401). (pyMzLib #75; *needs the pyMzLib 0.4.0 bridge*.)
* `flashlfq_quantify()` runs FlashLFQ label-free quantification with match-between-runs.
* `flashlfq_median_polish()` re-runs FlashLFQ's protein median polish from a
  `QuantifiedPeptides.tsv` under a new experimental design, without re-reading spectra.
* `flashlfq_quantify()` no longer warns about mzLib#1111 when `max_threads` is not 1: mzLib#1155
  fixed it inside the pinned bridge. The default stays 1 here (pyMzLib and mzLibRust use -1).
* `flashlfq_median_polish()` with no design no longer hangs in terminal R. A call that sends the
  bridge nothing on stdin now gives it an empty stdin rather than R's own, which a terminal never
  closes. Every verb goes through the same transport, so none can wait on it (pyMzLib #73).

## Differential abundance

* New `stats_fit()` fits one linear model per feature of a feature-by-sample table and tests each
  named coefficient with limma's empirical-Bayes moderated t, Benjamini-Hochberg adjusted:
  `lmFit()` then `eBayes(legacy = TRUE)`, computed by mzLib's `LinearModel` and `EmpiricalBayes`.
  A feature that cannot be fitted is reported with its reason, never dropped, and
  `residual_df_differ` says when default limma would use a different prior estimator.
* New `stats_adjust()` Benjamini-Hochberg adjusts p-values from anywhere, keeping every position;
  `NA` is untested and not counted in m.
* New `stats_meta()` pools one estimate per study into a DerSimonian-Laird random-effects estimate
  per feature, with direction agreement and leave-one-out sensitivity.
* The tests hold all three to limma, metafor and `stats::p.adjust()` to 1e-8 relative - against
  their recorded output always, and against limma and metafor run in the test when they are
  installed (both are now in `Suggests`). The reference tables ship in `extdata/stats/` so the
  examples and the new article read real files (mzLib #1341, #1357; *needs the pyMzLib 0.3.0
  bridge*).

## Protein databases

* New `proteins_read()` reads UniProt XML or FASTA databases into one row per protein - organism,
  NCBI taxon, gene names, length, mass - with GO terms and Ensembl gene links on request;
  `proteins_resolve_genes()` resolves proteins to stable Ensembl gene ids against a gene set you
  pin; `proteins_classify_peptides()` calls each peptide `Unique`, `SharedWithinGene`,
  `SharedAcrossGenes` or `NotInDatabase`, with I and L one residue. One database or many, in one
  bridge call (mzLib #1336, #1338, #1348; *needs the pyMzLib 0.2.0 bridge*).
* New `proteins_annotate_go()` annotates a stored MetaMorpheus protein-group table with Gene
  Ontology terms: one row per (group, term) that any member holds, directly or through an
  ancestor, naming the members that carry it, so consensus and direct-only views are filters on
  the rows. Every non-decoy group gets a row, a term-less one saying why. It reads the go.obo you
  name and never downloads one; `proteins_update_go()` is the one function that fetches the
  current release, keeping the previous file beside it (mzLib #1353, #1366; *needs the pyMzLib
  0.3.0 bridge*).

## Documentation

* Every help page for a wire verb now carries that verb's facts, generated from the bridge's
  per-verb spec rather than retyped: what it wraps in mzLib, every parameter with its unit, range
  and default, every returned field and column with its unit and what `NA` means, the errors it
  raises as condition classes, caveats, the same verb in pyMzLib and mzLibRust, references and
  `since`. `scripts/name-parity.R` checks the hand-written `@param` and `@return` text against the
  same specs, and CI fails on any difference.
* Every exported function has an example, and `R CMD check` runs them all: a stand-in bridge
  answers each call from a recording of the real one, so the examples need no bridge, .NET or
  network.
* New `?mzlib_error` documents the condition classes and how to handle each.
* `mzlibr_install_bridge()`, its help page and the README give the bridge's real size: about
  65 MB to download and about 110 MB unpacked, measured on pyMzLib 0.7.0's assets. They said
  about 140 MB, a figure from before the bridge stopped carrying libtorch.
* A pkgdown site, with an article for each module.
* The articles are now vignettes, so `R CMD check` builds and runs every one of them, on every
  platform, against the replay bridge. Each opens with a question -> function -> mzLib table,
  ends with what to cite (rendered from the specs' DOIs, each checked to resolve), and states no
  count of formats or verbs and no mzLib version in its prose; `scripts/docs-lint.R` holds all of
  that in CI.
* The articles teach on real data recorded from the real bridge: a FlashLFQ run with
  match-between-runs on mzLib's K562 pair (and median polish reproducing its protein intensities
  for every group), the whole albumin digest, PXD000001's live FTP listing, an RNA search's
  transcript groups, and the MALAT1 dilution series for `stats_fit()`.
* Corrected from those runs: the K562 peptide roll-up shows 21 of the 140 MBR transfers, not 52;
  `use_pep_q_value` filters nothing, it changes the q-value FlashLFQ carries; PXD000001's FTP
  tree now holds 14 files.
* `knitr` and `rmarkdown` join `Suggests`, for the vignettes. The package still imports nothing
  but base R.
