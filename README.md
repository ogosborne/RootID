# RootID

This package is part of a pipeline that uses GBS data to match unknown root samples to individual plants. It processes the output of the STACKS pipeline, identifies markers and haplotypes which are unique to each species and individual within a dataset and then matches these to root samples of unknown origin. The package also contains functions to validate and plot these results as a 3d 'root map'. 

To install, use:

```R
devtools::install_github("ogosborne/RootID")
```

See https://github.com/ogosborne/Caatinga_RootID for a full example of the pipeline, and see the package man pages for detailed usage information.

