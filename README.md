# TargetSearch

[![platforms](https://bioconductor.org/shields/availability/devel/TargetSearch.svg)](https://bioconductor.org/packages/devel/bioc/html/TargetSearch.html#archives)
[![rank](https://bioconductor.org/shields/downloads/devel/TargetSearch.svg)](http://bioconductor.org/packages/stats/bioc/TargetSearch/)
[![in Bioc](https://bioconductor.org/shields/years-in-bioc/TargetSearch.svg)](https://bioconductor.org/packages/devel/bioc/html/TargetSearch.html#since)
[![build](https://bioconductor.org/shields/build/devel/bioc/TargetSearch.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TargetSearch/)
[![updated](https://bioconductor.org/shields/lastcommit/devel/bioc/TargetSearch.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TargetSearch/)
[![dependencies](https://bioconductor.org/shields/dependencies/devel/TargetSearch.svg)](https://bioconductor.org/packages/devel/bioc/html/TargetSearch.html#since)

This packages provides a flexible, fast and accurate method for targeted
pre-processing of GC-MS data. The user provides a (often very large) set of GC
chromatograms and a metabolite library of targets. The package will
automatically search those targets in the chromatograms resulting in a data
matrix that can be used for further data analysis.

This software is part of the [Bioconductor](https://bioconductor.org) project.

https://bioconductor.org/packages/release/bioc/html/TargetSearch.html

This git repository mirrors the `TargetSearch` git repository on `Bioconductor`,
which is available here: https://code.bioconductor.org/browse/TargetSearch/

## Installation

Please follow the instructions on `TargetSearch`'s page on `Bioconductor` site
linked above. In short:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("TargetSearch")
```

## Citation

Cuadros-Inostroza A, Caldana C, Redestig H, Lisec J, Pena-Cortes H, Willmitzer
L, Hannah MA (2009). "TargetSearch - a Bioconductor package for the efficient
pre-processing of GC-MS metabolite profiling data." *BMC Bioinformatics*,
**10**, 428. doi:[10.1186/1471-2105-10-428](https://doi.org/10.1186/1471-2105-10-428).

## Bugs

Please report bugs at https://github.com/acinostroza/TargetSearch/issues

## License

The code is licensed under the GNU Public License version 2.0 or greater
(SPDX License ID: `GPL-2.0-or-later`).

To obtain a copy of the license, see: https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
