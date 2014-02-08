##----------------------------------------
## 'intgroupColors' returns a list with:
##   arrayColors:  color code for each array
##   key: a key explaining the mapping of factor values to colors
##----------------------------------------
intgroupColors = function(x)
{
  
  if (length(x$intgroup)>0) {
    colors = c(brewer.pal(9, "Set1"), brewer.pal(8, "Dark2"))
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
    if(nlevels(f) > n) {
      warning(sprintf("A factor was provided with %d levels, but the colour map used here has only %d colours. For the purpose of colouring, levels %d ('%s') to %d ('%s') are being collapsed. Please consider grouping together some of the levels of your factor of interest to reduce the number of levels, this might improve the legibility of the plots.\n",
                       nlevels(f), n, n, levels(f)[n], nlevels(f), levels(f)[nlevels(f)]))
      wipe = (as.integer(f) > n)
      f[wipe] = levels(f)[n] = "other"
      f = factor(f) ## remove extra factor levels
    }
    return(f)
  }
