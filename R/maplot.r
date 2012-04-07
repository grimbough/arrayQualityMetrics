## This function processes the parameters 'maxNumArrays' and
## 'nrColumns', and it is used by aqm.maplot as well as aqm.spatial

processMaxNumArrays = function(numArrays, maxNumArrays, nrColumns, stat)
{
  stopifnot(length(maxNumArrays)==1, maxNumArrays>=1,
            length(nrColumns)==1, nrColumns>=1)

  ## make a multiple of nrColumns
  maxNumArrays = nrColumns*ceiling(maxNumArrays/nrColumns)

  if(numArrays <= maxNumArrays)
  {
      j = seq_len(numArrays)
      lay = c(nrColumns, ceiling(numArrays/nrColumns))  # number of columns, number of rows
      legOrder = ""
  } else {
      first = seq_len(maxNumArrays/2)
      last  = numArrays + 1 - rev(first)
      j = order(stat, decreasing=TRUE)[c(first, last)]
      lay = c(nrColumns, ceiling(maxNumArrays/nrColumns))
      legOrder = sprintf("Shown are first the %d arrays with the highest values of %s, then the %d arrays with the lowest values. ",
          length(first), attr(stat, "name"), length(last))
    }

  list(j=j, lay=lay, legOrder=legOrder)
}

aqm.maplot = function(x, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4, ...)
{
  stopifnot(length(Dthresh)==1, Dthresh>0,
            length(subsample)==1, subsample>0)

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
  attr(stat, "name") = "<i>D<sub>a</sub></i>"

  out = new("outlierDetection",
    statistic = stat,
    threshold = Dthresh,
    which = which(stat > Dthresh),
    description = c(attr(stat, "name"), "fixed"))

  selected = processMaxNumArrays(x$numArrays, maxNumArrays, nrColumns, stat)

  xlim = quantile(A, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  ylim = quantile(M, probs=1e-4*c(1,-1)+c(0,1), na.rm=TRUE)
  panelNames = sprintf("array %d (D=%4.2f)", selected$j, stat[selected$j])

  i = seq(along=selected$j)
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
    panel = function(x, y, ...) panel.smoothScatter(x=A[, selected$j[x]], y=M[, selected$j[y]], raster=TRUE, nbin=250, ...),
    as.table = TRUE,
    layout = selected$lay,
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
    "Outlier detection was performed by computing Hoeffding's statistic <i>D<sub>a</sub></i> on the joint distribution of A and M for each array. ",
    selected$legOrder,
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
      size = with(selected, c(w=3*lay[1], h=2.5*lay[2])),
      colors  = x$arrayColors)
}

