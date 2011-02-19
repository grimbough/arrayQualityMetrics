## maxNumberOfArraysForDrawingDendrogram: draw dendrogram alongside the heatmap
## matrix only if the number of arrays is less than or equal this number.
## When the number of arrays is large, dendrograms become unreadable; also,
## there seems to be a bug in \code{\link[latticeExtra:dendrogramGrob]{dendrogramGrob}} which
## apparently can be circumvented by not drawing dendrograms that are too large. 

arrayQualityMetricsGlobalParameters = list(
  dpi = 72,
  maxNumberOfArraysForShowingArrayMetadataByDefault = 20,
  maxNumberOfArraysForDrawingDendrogram = 20
)
