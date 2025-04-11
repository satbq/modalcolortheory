# modalcolortheory

## Overview

This repository hosts R code relevant to the article "Modal Color Theory" in the *Journal of Music Theory* 69/1 (2025): 1-49.
The code is frozen as of April 2025 and should correspond to the comments made in the print version of that article.

However, I'd strongly encourage you to use the package [musicMCT](https://github.com/satbq/musicMCT) instead. That package 
includes documentation, explanatory vignettes, and ongoing development of new features.

The main value of this repository (`modalcolortheory`) aside from compatibility with the *JMT* article's text is that it
hosts the files representative_scales.rds, representative_signvectors.rds, and color_adjacencies.rds. Because these are large, they need
to remain separate from the package `musicMCT` and so they're staying in this repository.

Otherwise, `musicMCT` is better to use in pretty much every way. And do check out its [intro vignette](https://satbq.github.io/musicMCT/articles/musicMCT.html),
too!


## Changes between `modalcolortheory` and `musicMCT` that affect code in my *JMT* article

* Rather than following the instructions of footnote 7, follow the installation instructions in `musicMCT`'s readme file.
* Footnote 43 refers to `ineqmats.rds` as a separate file, but the inequality matrices are contained inside the package `musicMCT` and don't need to be downloaded separately.
* Footnote 49 uses intervals such as `just_wt` and `just_p4`. Rather than defining these as global variables, `musicMCT` uses a function `j()` to access them. Instead of `just_wt`, use `j(t)` or `j(2)`; instead of `just_p4`, use `j(4)`. The entire 5-limit just diatonic scale is available as `j(dia)`.
* The function `colornum()` returns `NULL` by default unless you've downloaded and loaded `representative_signvectors`. In `modalcolortheory` it simply returns an error if it doesn't have access to the right file.
* Footnote 74 points you to a Google Drive link for `color_adjacencies.rds`. That file is now available here on GitHub (in this repository, `modalcolortheory`, not in `musicMCT`).