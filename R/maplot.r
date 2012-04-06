aqm.maplot = function(x, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4, ...)
{
  stopifnot(length(maxNumArrays)==1, maxNumArrays>=1,
            length(Dthresh)==1, Dthresh>0,
            length(subsample)==1, subsample>0)

  maxNumArrays = nrColumns*ceiling(maxNumArrays/nrColumns)  ## make a multiple of 4

  if(x$nchannels==1)
    {
      stopifnot(identical(x$M, x$A))
      medArray = rowMedians(x$M, na.rm=TRUE)
      M =  x$M - medArray
      A = (x$M + medArray)/2
    } else {
      M = x$M
      A = x$A
      stopifnot(identical(dim(M), dim(A)))
    }

  ## For each array j, compute the D statistic from Hoeffding's test for independence
  ## and sort / detect outlier arrays by the value of the test statistic. We do subsampling
  ## to save time
  if(nrow(M)>subsample)
    {
      sel = sample(nrow(M), subsample)
      sM = M[sel, ]
      sA = A[sel, ]
    } else {
      sM = M
      sA = A
    }
  stat = sapply(seq_len(x$numArrays), function(j)
    {
      hoeffd(sA[,j], sM[,j])$D[1,2]
    })
  out = new("outlierDetection",
    statistic = stat,
    threshold = Dthresh,
    which = which(stat > Dthresh),
    description = c("Hoeffding's statistic <i>D<sub>a</sub></i>", "fixed"))

  ## Plot maximally maxNumArrays scatterplots
  if(x$numArrays <= maxNumArrays)
    {
      whj = seq_len(x$numArrays)
      lay = c(nrColumns, ceiling(x$numArrays/nrColumns))  # number of columns, number of rows
      legOrder = ""
    } else {
      first = seq_len(maxNumArrays/2)
      last  = x$numArrays + 1 - rev(first)
      whj = order(stat, decreasing=TRUE)[c(first, last)]
      lay = c(nrColumns, ceiling(maxNumArrays/nrColumns))
      legOrder = sprintf("Shown are first% d arrays with the highest values of <i>D<sub>a</sub></i>, then the %d arrays with the lowest values. ", length(first), length(last))
    }

  xlim = quantile(A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  ylim = quantile(M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  panelNames = sprintf("array %d (D=%4.2f)", whj, stat[whj])

  i = seq(along=whj)
  df = data.frame(
    i = factor(i, levels = i),
    px = i,
    py = i)

  ma = xyplot(py ~ px | i,
    data = df,
    xlim = xlim,
    ylim = ylim,
    xlab = "A",
    ylab = "M",
    panel = function(x, y, ...) panel.smoothScatter(x=A[, whj[x]], y=M[, whj[y]], raster=TRUE, nbin=250, ...),
    as.table = TRUE,
    layout = lay,
    asp = "iso",
    strip = function(..., bg, factor.levels) strip.default(..., bg ="#cce6ff", factor.levels = panelNames))

  vv = if(length(out@which)==1)
    c("One array", "was", "") else
    c(paste(length(out@which), "arrays"), "were", "s")
  outliertext = sprintf("%s had <i>D<sub>a</sub></i>&gt;%g and %s marked as outlier%s. ",
                         vv[1],                    Dthresh, vv[2],             vv[3])

  legend = paste("The figure <!-- FIG --> shows MA plots. M and A are defined as:<br>",
    "M = log<sub>2</sub>(I<sub>1</sub>) - log<sub>2</sub>(I<sub>2</sub>)<br>",
    "A = 1/2 (log<sub>2</sub>(I<sub>1</sub>)+log<sub>2</sub>(I<sub>2</sub>)),<br>",
    if(x$nchannels == 1)
      paste("where I<sub>1</sub> is the intensity of the array studied,",
            "and I<sub>2</sub> is the intensity of a \"pseudo\"-array that consists of the median across arrays.") else
      "where I<sub>1</sub> and I<sub>2</sub> are the intensities of the two channels.",
    " Typically, we expect the mass of the distribution in an MA plot to be concentrated along the M = 0 axis, ",
    "and there should be no trend in M as a function of A. If there is a trend in the lower range of A, this often ",
    "indicates that the arrays have different background intensities; this may be addressed by background correction. ",
    "A trend in the upper range of A can indicate saturation of the measurements; in mild cases, this may be addressed ",
    "by non-linear normalisation (e.g. quantile normalisation).<br>",
    "Outlier detection was performed by computing Hoeffding's statistic <i>D<sub>a</sub></i> on the joint distribution ",
    "of A and M for each array. ",
    legOrder,
    "The value of <i>D<sub>a</sub></i> is shown in the panel headings. ",
    outliertext,
    "For more information on Hoeffing's <i>D</i>-statistic, please see the manual page of the function ",
    "<tt>hoeffd</tt> in the <tt>Hmisc</tt> package.",
    sep="")

  new("aqmReportModule",
      plot = ma,
      section = "Individual array quality",
      title = "MA plots",
      id = "ma",
      legend = legend,
      outliers = out,
      size = c(w=10, h=6),
      colors  = x$arrayColors)
}

