# Download a bridge executable into a local cache

mzLibR needs a bridge executable and does not ship one: it is about 140
MB, which no CRAN package may carry. This fetches the one built for your
platform, verifies it against a recorded SHA-256, and unpacks it where
[`mzlibr_bridge_path`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_path.md)
will find it.

## Usage

``` r
mzlibr_install_bridge(version = MZLIB_BRIDGE_VERSION, destination = NULL, consent = NA,
  overwrite = FALSE, url = NULL, sha256 = NULL, quiet = FALSE, timeout = 1800)
```

## Arguments

- version:

  The pyMzLib release to take the payload from. Every pyMzLib release
  publishes the raw bridge for each platform as
  `mzlib-bridge-<rid>.tar.gz`, together with a `SHA256SUMS` manifest;
  this fetches the one matching your platform.

- destination:

  Directory to unpack into. Defaults to R's per-user cache directory for
  this package, in the same location
  [`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) would
  choose.

- consent:

  Set `TRUE` to confirm the download in a non-interactive session. Left
  as `NA` an interactive session asks, and a non-interactive one
  refuses.

- overwrite:

  Whether to replace a bridge already present at the destination.

- url:

  An explicit URL, overriding `version`. When you pass this, pass
  `sha256` too.

- sha256:

  The expected SHA-256 of the download, lowercase hex. Required with
  `url`.

- quiet:

  Suppress progress output.

- timeout:

  Seconds to allow for the download. The default is generous because the
  payload is large and R's own default of 60 seconds is not enough for
  it on most connections.

## Details

**It is never called for you.** CRAN policy forbids a package
downloading anything at install time or writing outside the session's
temporary directory without consent, and both are right — so this asks
in an interactive session and requires `consent = TRUE` otherwise.

If you already have a bridge, you do not need this at all: set
`MZLIB_BRIDGE` or `options(mzlibr.bridge=)` to point at it. That is also
how you would relink a modified mzLib, which LGPL section 4 requires
this package to permit.

## Value

The path of the installed bridge executable, invisibly.

## See also

[`mzlibr_bridge_path`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_path.md),
[`mzlibr_bridge_version`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_bridge_version.md)

## Examples

``` r

# \donttest{
# Downloads the bridge (about 60 MB) after asking; needs the network.
if (interactive()) mzlibr_install_bridge()
# }
```
