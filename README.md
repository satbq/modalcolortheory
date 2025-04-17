# modalcolortheory

## Overview

This repository hosts R code relevant to the article "Modal Color Theory" in the *Journal of Music Theory* 69/1 (2025): 1-49.
The code is frozen as of April 2025 and should correspond to the comments made in the print version of that article.

However, I'd strongly encourage you to use the package [musicMCT](https://satbq.github.io/musicMCT/) instead. That package 
includes documentation, explanatory vignettes, and ongoing development of new features.

The main value of this repository (`modalcolortheory`) aside from compatibility with the *JMT* article's text is that it
hosts the files representative_scales.rds, representative_signvectors.rds, and color_adjacencies.rds. Because these are large, they need
to remain separate from the package `musicMCT` and so they're staying in this repository.

Otherwise, `musicMCT` is better to use in pretty much every way. And do check out its [intro vignette](https://satbq.github.io/musicMCT/articles/musicMCT.html),
too!


## Changes between `modalcolortheory` and `musicMCT` that affect code in my *JMT* article

* Follow the [installation instructions in `musicMCT`'s readme file](https://github.com/satbq/musicMCT#readme) rather than those in footnote 7.
* Footnote 43 refers to `ineqmats.rds` as a separate file, but the inequality matrices are contained inside the package `musicMCT` and don't need to be downloaded separately.
* Footnote 49 uses intervals such as `just_wt` and `just_p4`. Rather than defining these as global variables, `musicMCT` uses a function `j()` to access them. Instead of `just_wt`, use `j(t)` or `j(2)`; instead of `just_p4`, use `j(4)`. The entire 5-limit just diatonic scale is available as `j(dia)`.
* The function `musicMCT::colornum()` returns `NULL` by default unless you've downloaded and loaded `representative_signvectors`. In `modalcolortheory` it simply returns an error if it doesn't have access to the right file.
* Footnote 74 points you to a Google Drive link for `color_adjacencies.rds`. That file is now available here on GitHub (in this repository, `modalcolortheory`, not in `musicMCT`).

## Format of the datasets

### color_adjacencies.rds

This file is a [list](https://cran.r-project.org/doc/manuals/r-release/R-intro.html#index-list) of [adjacency lists](https://en.wikipedia.org/wiki/Adjacency_list). That is, each entry `color_adjacencies[[i]]` is a list of the adjacencies between colors of scales with i notes, ordered according to color number. Thus `color_adjacencies[[i]][[j]]` tells you which colors are adjacent to the i-note scale with color number j. Adjacency to color number 0, the perfectly even scale, is not listed because it's true of all scales.

One- and two-note scales are too simple to have much structure, so `color_adjacencies[[1]]` and `color_adjacencies[[2]]` are both `NULL`. For monads, there is only a single scale structure. For dyads, there are only three colors, arranged in the following graph:

    1 <-> 0 <-> 2

Thus the only adjacencies involve the perfectly even scale (0, 6), which is excluded from the list.

The smallest adjacency list that's of interest is the one for trichords, which is small enough to display in full:

```
[[1]]
[1] 2 6

[[2]]
[1] 1 3

[[3]]
[1] 2 4

[[4]]
[1] 5 3

[[5]]
[1] 7 4

[[6]]
[1] 1 8

[[7]]
[1]  5 12

[[8]]
[1] 9 6

[[9]]
[1] 10  8

[[10]]
[1]  9 11

[[11]]
[1] 10 12

[[12]]
[1]  7 11
```

This tells us that trichordal color 1, which includes the scale (0, 3, 7), is adjacent to colors 2 and 6. Color 2 includes the scale (0, 2, 7) and color 6 includes the scale (0, 3.5, 7). These adjacencies make sense if we consider the scales in terms of their step sizes (**S** for small, **M** for medium, and **L** for large): color 1 has the step pattern SML, color 2 has the step pattern SLL, and color 6 has the step pattern SSL. Thus the adjacency between 1 and 2 involves enlarging **M** so that it equals **L**. Similarly, the adjacency between 1 and 6 involves shrinking **M** until it equals **S**. 

If we look at entry [[6]], we see that color 6 is adjacent to 1 (as we just learned) and also to 8. Color 8 is represented by the scale (0, 4, 7), which has the step word MSL. Thus moving from color 6 to color 8 involves distinguishing color 6's two **S** steps such that the first is larger than the second.

In musical terms, color 1 includes the minor triad, color 8 includes the major triad, and color 6 includes the neutral triad (C, E half-flat, G). The major and minor triads differ as scale structures because they subdivide their perfect fifth differently (i.e. with their larger their above or below, respectively). The neutral triad mediates between the two of them by dividing the perfect fifth exactly in half; thus the neutral triad is adjacent to both major and minor.

The file color_adjacencies so far includes adjacency lists for scales up through heptachords. Future work will hopefully create lists for scales of higher cardinality, but for now color_adjacencies[[8]] and above don't exist.

Since the adjacency list for heptachords is ~1.8 million entries long, you don't want to look at it directly. Instead, you'll want to explore it using tools designed for the study of network structures. I recommend the following:

```
library(igraph)
hepta_network <- igraph::graph_from_adj_list(color_adjacencies[[7]]) # May take a minute or two to construct
```

You can then use the tools of the igraph package, like `igraph::bfs()` to explore the structure of the space.


### representative_signvectors.rds

Like `color_adjacencies`, `representative_signvectors` is structured as a list whose nth entry corresponds to scales with n notes. Thus, for instances, representative_signvectors[[3]] contains the sign vectors for all 12 non-white trichordal scales. Each sign vector is stored as a comma-separated character string. For instance, trichord color 1 has the sign vector (-1, -1, -1) so if we access it in R we should see:

```
representative_signvectors[[3]][1]
#> "-1, -1, -1"
```

For each cardinality, the sign vectors are sorted lexicographically, i.e. those with -1 in the first position come before all those with 0 or 1 in the first position; then those with -1 in the second position come before 0 or 1 in the second position; and so on. Therefore the first sign vector in every cardinality is a list of -1s and the last sign vector is always a list of 1s.

For now, `representative_signvectors` only includes scales of 7 notes or smaller. Functions that depend on `representative_signvectors` (in particular `colornum()`) won't work for larger scales.

### representative_scales.rds

Like the previous two files, `representative_scales` is structured as a list whose nth entry corresponds to scales with n notes. This means that representative_scales[[n]] is an n-by-k matrix whose columns represent the k distinct colors in the space. (The perfectly even scale isn't included--it can be created by calling `edoo(n)`.) The first column represents color number 1, and so on. Thus the order of scales in `representative_scales` matches the order of sign vectors in `representative_signvectors`. It would be possible (but potentially time consuming!) to rederive `representative_signvectors` from `representative_scales` with the following:

```
# Don't try this at home for large n!
representative_signvectors[[n]] <- apply(apply(representative_scales[[n]], 2, signvector), 2, toString)
```

Additionally, all the representative scales have been normalized so that they're all equally even (in the sense that `evenness()` returns 1 when applied to every representative scale).

As with `color_adjacencies` and `representative_signvectors`, `representative_scales` only covers scales with 7 notes or fewer.