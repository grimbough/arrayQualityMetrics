## ---------------------------------------------------------------
## Functions for defining colours for boxplots, heatmap colour
## side bar, density plots etc.
## -----------------------------------------------------------------


##----------------------------------------
## 'intgroupColours' returns a list with:
##   arrayColors:  color code for each array
##   key: a key explaining the mapping of factor values to colours
##----------------------------------------
intgroupColours = function(x, withOpacity = FALSE)
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

  if(withOpacity)
    cols = addOpacity(cols, n)

  list(arrayColours = cols, key = key)
}

addOpacity = function(cols, n){
  stopifnot(all(nchar(cols)==7))  ## expecting something like #80d020
  opacity = if(n > 50) 0.125 else (0.8-n*0.0135)  ## the more arrays, the more transparency
  opacity = as.hexmode(as.integer(255 * opacity))
  paste(cols, opacity, sep="")
}

##----------------------------------------
## 'outlierColours' returns a list with:
##   arrayColors:  color code for each array
##----------------------------------------
outlierColours = function(outliers, n, withOpacity = FALSE)
{

  colours = c("outliers"="#ff0000", "non-outliers"="#000000")
  if(withOpacity) colours[2] = addOpacity(colours[2], n) 
  cols = ifelse(seq_len(n) %in% outliers, colours[1], colours[2])
    
  key = list(
    rect = list(col = colours),
    text = list(names(colours)),
    rep = FALSE)
  
  list(arrayColours = cols, key = key)
}
