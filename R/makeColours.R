##----------------------------------------
## 'intgroupColours' returns a list with:
##   arrayColors:  color code for each array
##   key: a key explaining the mapping of factor values to colours
##----------------------------------------
intgroupColours = function(x)
{

  n = x$numArrays
  cols = rep("#1F78B4", n)
  key = NULL
  
  if (length(x$intgroup)>0) {
    groups  = as.factor(x$pData[[x$intgroup[1]]])
    igroups = as.integer(groups)
    colours = brewer.pal(9, "Set1")
    if(nlevels(groups) > length(colours)) {
      warning(sprintf("'intgroup[1]' has %d levels, but only the first 9 are used for colouring.",
                      nlevels(groups)))
      igroups[igroups > length(colours)] = length(colours)+1
      colours = c(colours, "#101010")
    } else {
      colours = colours[seq_len(nlevels(groups))]
    }
    cols = colours[igroups]
    key = list(
      rect = list(col = colours),
      text = list(levels(groups)),
      rep = FALSE)
  } # if

  list(arrayColours = cols, key = key)
}

