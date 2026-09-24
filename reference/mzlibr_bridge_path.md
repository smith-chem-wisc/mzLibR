# Path of the bridge executable mzLibR will use

Resolution order: the `mzlibr.bridge` option, then the `MZLIB_BRIDGE`
environment variable, then a bridge downloaded by
[`mzlibr_install_bridge`](https://smith-chem-wisc.github.io/mzLibR/reference/mzlibr_install_bridge.md).

## Usage

``` r
mzlibr_bridge_path()
```

## Value

A single file path.

## Examples

``` r

# Where mzLibR will look; an error naming three remedies when there is no bridge yet.
try(mzlibr_bridge_path())
#> Error : No mzLib bridge executable found.
#> 
#> Three ways to fix it, cheapest first:
#>   1. mzlibr_install_bridge() - downloads one into /home/runner/.cache/R/mzLibR/linux-x64/mzlib-bridge
#>   2. Sys.setenv(MZLIB_BRIDGE = "/path/to/mzlib-bridge") - for a bridge you already have,
#>      for example the one pyMzLib stages under pkg/python/src/pymzlib/_dotnet/<rid>/.
#>   3. options(mzlibr.bridge = "/path/to/mzlib-bridge") - same thing, scoped to this session.
#> 
#> Overriding the bridge is also how you relink a modified mzLib without rebuilding this package, which LGPL section 4 requires mzLibR to allow.
```
