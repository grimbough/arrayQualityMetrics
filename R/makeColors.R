##----------------------------------------
## 'intgroupColors' returns a list with:
##   arrayColors:  color code for each array
##   key: a key explaining the mapping of factor values to colors
##----------------------------------------
intgroupColors = function(x)
{
  
  if (length(x$intgroup)>0) {
    colors = brewer.pal(9, "Set1")
    fac  = as.factor(x$pData[[x$intgroup[1]]])
    fac  = maximumLevels(fac, n = length(colors)) ## make sure that factor has at most n levels
    colors = colors[seq_len(nlevels(fac))]
    ac = colors[as.integer(fac)]

    key = list(
      rect = list(col = colors),
      text = list(levels(fac)),
      rep = FALSE)
    
  } else {
    key = NULL
    ac = rep("#1F78B4", x$numArrays)
  }

  list(arrayColors = ac, key = key)
}

maximumLevels = function(f, n)
  {
    if(nlevels(f) >= n) {
      warning(sprintf("A factor was provided with %d levels, but only the first %d were used for coloring.",
                       nlevels(f), n))
      wipe = (as.integer(f) > n)
      f[wipe] = levels(f)[n] = "other"
      f = factor(f) ## remove unneeded factor levels
    }
    return(f)
  }
