# Find PRIDE Archive projects by keyword

Find PRIDE Archive projects by keyword, with every page fetched and no
accession repeated: the discovery step that produces the accessions
every other pride verb takes.

**The discovery entry point.** Every other `pride_` function takes an
accession you already have; this is the one that produces them, so you
can go from a subject to a dataset without leaving R. Paging is handled
for you: however many pages the result set spans, you get one
data.frame, with no accession repeated.

## Usage

``` r
pride_search(keyword, page_size = 100, timeout = 300)
```

## Arguments

- keyword:

  What to search for, e.g. `"phosphoproteome"`, as one string of 1 to
  1000 characters. PRIDE matches it across titles, descriptions,
  keywords, organisms, references and more; the `matched_fields` column
  says which fields actually matched each hit. A keyword beginning with
  `-` is refused, since the bridge would read it as an option.

- page_size:

  How many projects to request per underlying API call. Changes how many
  requests the fetch takes, never what comes back.

- timeout:

  Seconds to allow for the whole fetch, or `NULL` to wait indefinitely.

## Value

A data.frame with one row per matching project, in PRIDE's ranking
order. **Zero rows is a real answer** - PRIDE reports no hits as an
empty result, and unlike
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
this does not raise `mzlib_project_not_found`, because no accession here
could have been a typo.

Text columns: `accession`, `title`, `project_description`,
`sample_processing_protocol`, `data_processing_protocol`, `doi`,
`submission_type` and `sdrf`. Dates: `submission_date`,
`publication_date` and `updated_date`, as `Date` (bare calendar dates,
never timestamps), `NA` when PRIDE reported none. Popularity:
`download_count` in downloads, `avg_downloads_per_file` in downloads per
file, `percentile` a download-popularity percentile within PRIDE, and
`bot_count`, `hub_count` and `organic_count` in downloads split by
traffic kind. List columns, one character vector per hit:
`project_tags`, `keywords`, `submitters`, `lab_pis`, `affiliations`,
`instruments`, `softwares`, `quantification_methods`,
`sample_attributes`, `organisms`, `organism_parts`, `diseases`,
`references`, `experiment_types`, `project_file_names` and
`other_omics_links`; `highlights`, a named list per hit (PRIDE field to
matched snippets, the terms wrapped in `<em>`); `matched_fields`, its
sorted names; and `yearly_downloads`, a data.frame per hit of `year` and
`count`.

## A hit is not a project's metadata

PRIDE serves search from a separate index in which every
controlled-vocabulary field is **flattened to a display string**:
instruments arrive as `"Q Exactive"` with no accession, contacts as
names, publications as pre-formatted citation strings, and `sdrf` as one
space-joined bag of term values - not a file, name or URL. Follow the
`accession` for the vocabulary. `project_file_names` is not the manifest
either: use
[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md)
or
[`pride_list_ftp_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_ftp_files.md)
to act on files.

## Zero and empty mean not reported

PRIDE sends no nulls, so an absent count arrives as **0** and an absent
list as empty. A `download_count` of 0 means "not reported", never
"nobody downloaded it", and several fields are genuinely sparse: sampled
over 1,600 hits, `project_tags` was populated on 2.6%, `sdrf` on 2.4%,
and the bot/hub/organic counts on under half. `keywords` may hold empty
and whitespace-only strings, about 9% of hits; they are passed through,
so filter before joining.

## A live index

PRIDE pages a live index with no stable cursor. A project published
mid-fetch is deduplicated, but one removed mid-fetch can fall between
two pages and be missed. A result set that fits on one page cannot be
affected.

## Wraps

Wire verb `pride search`. Generated from the bridge's verb spec
`pride.search.yaml` (bridge commit `5db922d4cfe1`) by
`scripts/build-man.R`; the spec owns these facts, and all three bindings
render the same ones.

- [`PrideArchiveClient.SearchProjectsAsync`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideArchiveClient.MaxKeywordLength`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideArchiveClient.cs` at mzLib
  `23c2490e`

- [`PrideProjectSearchResult`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideProjectSearchResult.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideProjectSearchResult.cs` at
  mzLib `23c2490e`

- [`PrideYearlyDownloadCount`](https://github.com/smith-chem-wisc/mzLib/blob/23c2490e/mzLib/UsefulProteomicsDatabases/PrideYearlyDownloadCount.cs)
  in `mzLib/UsefulProteomicsDatabases/PrideYearlyDownloadCount.cs` at
  mzLib `23c2490e`

## Parameters: units, ranges and defaults

- `keyword`:

  string; required; range `1 to 1000 characters`. What to search for,
  matched by PRIDE across titles, descriptions, keywords, organisms,
  references and more. Longer than 1000 characters (mzLib
  MaxKeywordLength) is a usage error, because PRIDE answers one with
  HTTP 500, which cannot be told apart from an outage.

- `page_size` (wire `--page-size`):

  int; in **projects**; default `100`; range `>= 1`. Hits requested per
  underlying API call. Changes how many requests the fetch takes, never
  what comes back.

## Returned fields

Each field with its type, its unit, and what `NA` means when it is `NA`.

- `keyword`:

  string; never `NA`. The keyword, echoed.

- `result_count`:

  int; in **projects**; never `NA`. Entries in results. 0 is a real
  answer (no hits), not an error.

- `results`:

  object\[\]; never `NA`. One entry per matching project, deduplicated
  by accession, in PRIDE's ranking order; fields under result.columns.

## Columns

- `accession`:

  string; never `NA`. The project accession: usually PXD, but also
  legacy PRD and affinity PAD. Treat it as an opaque key.

- `title`:

  string; never `NA`. Project title.

- `project_description`:

  string; never `NA`. Submitter's free-text description.

- `sample_processing_protocol`:

  string; never `NA`. Free text.

- `data_processing_protocol`:

  string; never `NA`. Free text.

- `doi`:

  string; never `NA`. Dataset DOI, or "" when PRIDE has minted none.

- `submission_type`:

  string; never `NA`. COMPLETE, PARTIAL, AFFINITY or PRIDE (legacy); not
  a closed set.

- `sdrf`:

  string; never `NA`. The project's SDRF term VALUES flattened by the
  search index into one space-joined string. Not a file, name or URL; ""
  on most hits.

- `submission_date`:

  string; `NA` when PRIDE reported no date (mzLib's default DateTime).
  Bare calendar date yyyy-MM-dd: no time, no offset, because none was
  sent.

- `publication_date`:

  string; `NA` when PRIDE reported no date. As submission_date.

- `updated_date`:

  string; `NA` when PRIDE reported no date. As submission_date.

- `project_tags`:

  string\[\]; never `NA`. PRIDE's coarse classification tags. Sparse
  (2.6% of 1,600 sampled hits).

- `keywords`:

  string\[\]; never `NA`. Submitter keywords. May contain empty or
  whitespace-only strings (about 9% of hits); passed through, not
  filtered.

- `submitters`:

  string\[\]; never `NA`. Display names, flattened from contact objects.

- `lab_pis`:

  string\[\]; never `NA`. Display names.

- `affiliations`:

  string\[\]; never `NA`. Display strings.

- `instruments`:

  string\[\]; never `NA`. Display names only; the CV accessions are not
  sent by this endpoint.

- `softwares`:

  string\[\]; never `NA`. Display names.

- `quantification_methods`:

  string\[\]; never `NA`. Display names, e.g. TMT.

- `sample_attributes`:

  string\[\]; never `NA`. Sample characteristic VALUES; which
  characteristic each describes is not recoverable from a hit.

- `organisms`:

  string\[\]; never `NA`. Display names, e.g. 'Homo sapiens (human)'.

- `organism_parts`:

  string\[\]; never `NA`. Display names.

- `diseases`:

  string\[\]; never `NA`. Display names.

- `references`:

  string\[\]; never `NA`. Each a pre-formatted citation string; no
  separate PubMed id or DOI.

- `experiment_types`:

  string\[\]; never `NA`. e.g. Data-independent acquisition.

- `project_file_names`:

  string\[\]; never `NA`. File names only: not the manifest (no sizes,
  categories or locations). Use pride files or pride ftp-files to act on
  files.

- `other_omics_links`:

  string\[\]; never `NA`. Links to related datasets in other omics
  repositories.

- `highlights`:

  map\<string,string\[\]\>; never `NA`. Why the project matched: PRIDE
  field name -\> matched snippets with the terms wrapped in \<em\>. Keys
  vary per hit and per query and cross unchanged (not snake_cased).

- `yearly_downloads`:

  object\[\]; never `NA`. {year, count}: year a string (e.g. "2025"),
  count downloads in that year. Empty when not reported.

- `download_count`:

  int; in **downloads**; never `NA`. Total downloads. 0 means not
  reported, never a measured zero.

- `avg_downloads_per_file`:

  float; in **downloads per file**; never `NA`. Mean downloads per file.
  0 means not reported.

- `percentile`:

  int; in **percentile**; never `NA`. Download-popularity percentile
  within PRIDE. 0 means not reported.

- `bot_count`:

  int; in **downloads**; never `NA`. Downloads attributed to crawlers. 0
  means not reported (populated on under half of hits).

- `hub_count`:

  int; in **downloads**; never `NA`. Downloads attributed to
  institutional or aggregating hubs. 0 means not reported.

- `organic_count`:

  int; in **downloads**; never `NA`. Downloads attributed to ordinary
  human traffic. 0 means not reported.

## Errors

Each is an R condition carrying the class shown and `mzlib_error`; see
[`mzlib_error`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlib_error.md).

- `mzlib_usage_error` (usage):

  keyword missing or blank; keyword longer than 1000 characters;
  page-size not an integer or \<= 0 (both checked by the bridge before
  mzLib)

- `mzlib_service_unavailable` (service_unavailable):

  the bridge's own classification (Program.ClassifyError; mzLib#1350 is
  not merged): a timeout or cancellation, a socket failure, a request
  that never got a response (refused connection, DNS, TLS), an HTTP 408,
  429 or 5xx, or a body cut off in transit

- `mzlib_bridge_error` (correctness):

  any other HTTP status; a paging contract violation (MzLibException: an
  identical page while total_records says more remain); the page limit
  (mzLib MaxPages) exceeded

## Caveats

- A hit is NOT a project's metadata: PRIDE serves search from an
  Elasticsearch projection in which every controlled-vocabulary field is
  flattened to a display string, contacts to names and references to
  citation strings. Follow the accession for the vocabulary.

- A zero or an empty list means 'not reported', never a measured zero.
  Several fields are genuinely sparse (sampled over 1,600 hits:
  project_tags 2.6%, sdrf 2.4%, other_omics_links 18%, the
  bot/hub/organic trio under half).

- Dates are bare calendar dates, unlike pride files' timestamps with
  offsets; a binding must not widen them to a midnight timestamp.

- PRIDE pages a live index with no stable cursor: a project published
  mid-fetch is deduplicated, but one removed mid-fetch can fall between
  two pages and be missed. A result set that fits on one page cannot be
  affected.

## Performance

Every call starts one bridge process, which costs a .NET start-up before
any work.

## Same verb in other bindings

- Python (pyMzLib): `pymzlib.pride.search`

- Rust (mzLibRust): `mzlib::pride::search_with` with `SearchOptions`

- R (mzLibR): `pride_search`

## Since

Wire protocol 1; pyMzLib 0.1.0; mzLibRust not yet shipped; mzLibR not
yet shipped.

## Not yet verified

The spec records these as open. They are listed rather than hidden:

- Not yet in mzLibRust or mzLibR. The names are intended, mirroring each
  binding's pride conventions: Rust mzlib::pride::search (defaults) and
  search_with + SearchOptions {page_size, timeout} like
  list_files_with + ListOptions; R pride_search(keyword, page_size =
  100, timeout = 300) like pride_list_files. The checker cannot verify
  them.

- The bindings also refuse a keyword beginning with '-' (it would be
  parsed as an option); the bridge cannot, since a value starting with
  '–' reaches it as a flag and the keyword then reads as missing.

- The ranking order of results is PRIDE's and is not promised stable
  across calls.

## References

- [doi:10.1093/nar/gkae1011](https://doi.org/10.1093/nar/gkae1011): The
  PRIDE Archive and its search service (Perez-Riverol et al., The PRIDE
  database at 20 years: 2025 update, NAR)

## See also

[`pride_list_files`](https://smith-chem-wisc.github.io/mzLibR/reference/pride_list_files.md),
which is usually what you want next.

## Examples

``` r

hits <- pride_search("plasmodium falciparum schizont")
hits[, c("accession", "submission_type", "publication_date")]
#>   accession submission_type publication_date
#> 1 PXD070842        COMPLETE       2026-06-01
#> 2 PXD020210        COMPLETE       2021-01-25
#> 3 PXD020189        COMPLETE       2021-01-25
#> 4 PXD008250         PARTIAL       2018-02-13
#> 5 PXD001684         PARTIAL       2015-04-23
#> 6 PXD000070        COMPLETE       2014-04-24
hits$organisms[[1]]
#> [1] "Homo sapiens (human)"                "Plasmodium falciparum (isolate 3d7)"
hits$matched_fields[[1]]
#> [1] "references" "title"     
```
