setClassUnion("aqmInputObj", c("ExpressionSet", "AffyBatch", "NChannelSet", "BeadLevelList"))

## set the method for arrayQualityMetrics
setGeneric("arrayQualityMetrics",
           function(expressionset,
                    outdir = getwd(),
                    force = FALSE,
                    do.logtransform = FALSE,
                    split.plots = FALSE,
                    intgroup = "Covariate",
                    grouprep = FALSE)
           standardGeneric("arrayQualityMetrics"))

## aqmInputObj
setMethod("arrayQualityMetrics",signature(expressionset = "aqmInputObj"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            olddir = getwd()
            on.exit(setwd(olddir))
            dircreation(outdir, force)
                
            arg = as.list(match.call(expand.dots = TRUE))
            name = arg$expressionset
            f = 1
            obj = list()
    
            datap = aqm.prepdata(expressionset, do.logtransform)

            obj$maplot = try(aqm.maplot(dataprep = datap))
            if(inherits(obj$maplot,"try-error"))
              warning("Cannot draw MA plots \n") else {
                obj$maplot$legend = gsub("<!-- FIG -->",f,obj$maplot$legend)
                f = f+1
              }

            if(inherits(expressionset, 'BeadLevelList') || ("X" %in% rownames(featureData(expressionset)@varMetadata) && "Y" %in% rownames(featureData(expressionset)@varMetadata)))
              {            
                obj$spatial =  try(aqm.spatial(expressionset = expressionset, dataprep = datap))
                if(inherits(obj$spatial,"try-error"))
                  warning("Cannot draw spatial distribution of intensities \n") else {
                    obj$spatial$legend = gsub("<!-- FIG -->",f,obj$spatial$legend)
                    f = f+1
                  }
                
                if((inherits(expressionset,'NChannelSet') && ("Rb" %in% colnames(dims(expressionset)) && "Gb" %in% colnames(dims(expressionset)))) || inherits(expressionset,'BeadLevelList'))
                  {
                    obj$spatialbg =  try(aqm.spatialbg(expressionset = expressionset, dataprep = datap))
                    if(inherits(obj$spatialbg,"try-error"))
                      warning("Cannot draw spatial distribution of background intensities \n") else {
                        obj$spatialbg$legend = gsub("<!-- FIG -->",f,obj$spatialbg$legend)
                        f = f+1
                      }
                  }                
              }
                        
            obj$boxplot = try(aqm.boxplot(dataprep = datap))
            if(inherits(obj$boxplot,"try-error"))
              warning("Cannot draw boxplots \n") else {
                obj$boxplot$legend = gsub("<!-- FIG -->",f,obj$boxplot$legend)
                f = f+1
              }

            obj$density = try(aqm.density(dataprep = datap))
            if(inherits(obj$density,"try-error"))
              warning("Cannot draw density plots \n") else {
                obj$density$legend = gsub("<!-- FIG -->",f,obj$density$legend)
                f = f+1
              }

            obj$heatmap = try(aqm.heatmap(expressionset = expressionset, dataprep = datap, intgroup))
            if(inherits(obj$heatmap,"try-error"))
              warning("Cannot draw heatmap \n") else {
                obj$heatmap$legend = gsub("<!-- FIG -->",f,obj$heatmap$legend)
                f = f+1
              }

            obj$meansd = try(aqm.meansd(dataprep = datap))
            if(inherits(obj$meansd,"try-error"))
              warning("Cannot draw Mean vs Standard Deviation \n") else {
                obj$meansd$legend = gsub("<!-- FIG -->",f,obj$meansd$legend)
                f = f+1
              }

            if(inherits(expressionset,'BeadLevelList'))
              warning("Cannot plot the probes mapping densities on a BeadLevelList object.")
            if(!inherits(expressionset,'BeadLevelList'))
              {
                if("hasTarget" %in% rownames(featureData(expressionset)@varMetadata))
                {
                  obj$probesmap = try(aqm.probesmap(expressionset = expressionset, dataprep = datap))
                  if(inherits(obj$probesmap,"try-error"))
                    warning("Cannot draw probes mapping plot \n") else {
                      obj$probesmap$legend = gsub("<!-- FIG -->",f,obj$probesmap$legend)
                      f = f+1
                    }
                }
              }
            

            if(inherits(expressionset, "AffyBatch"))
              {
                obj$rnadeg = try(aqm.rnadeg(expressionset))
                if(inherits(obj$rnadeg,"try-error"))
                  warning("Cannot draw the RNA degradation plot \n") else {
                    obj$rnadeg$legend = gsub("<!-- FIG -->",f,obj$rnadeg$legend)
                    f = f+1
                  }
        
                affyproc = aqm.prepaffy(expressionset, datap$sN)
                obj$rle = try(aqm.rle(affyproc))
                if(inherits(obj$rle,"try-error"))
                  warning("Cannot draw the RLE plot \n") else {
                    obj$rle$legend = gsub("<!-- FIG -->",f,obj$rle$legend)
                    f = f+1
                  }
                
                obj$nuse = try(aqm.nuse(affyproc))
                if(inherits(obj$nuse,"try-error"))
                  warning("Cannot draw the NUSE plot \n") else {
                    obj$nuse$legend = gsub("<!-- FIG -->",f,obj$nuse$legend)
                    f = f+1
                  }
        
                obj$qcstats =  try(aqm.qcstats(expressionset))
                if(inherits(obj$qcstats,"try-error"))
                  warning("Cannot draw the QCStats plot \n") else {
                    obj$qcstats$legend = gsub("<!-- FIG -->",f,obj$qcstats$legend)
                    f = f+1
                  }
                
                obj$pmmm = try(aqm.pmmm(expressionset))
                if(inherits(obj$pmmm,"try-error"))
                  warning("Cannot draw the Perfect Match versus MisMatch plot \n") else {
                    obj$pmmm$legend = gsub("<!-- FIG -->",f,obj$pmmm$legend)
                    f = f+1
                  }
              }

            aqm.writereport(name, expressionset, obj)
            return(invisible(obj))
          })



## RGList
setMethod("arrayQualityMetrics",signature(expressionset = "RGList"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a RGList and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method RGList

## marrayRaw
setMethod("arrayQualityMetrics",signature(expressionset = "marrayRaw"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "NChannelSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayRaw and cannot be converted automatically in a NChannelSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method marrayRaw

## MAlist
setMethod("arrayQualityMetrics",signature(expressionset = "MAList"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "ExpressionSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a MAList and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method MAList

## marrayNorm
setMethod("arrayQualityMetrics",signature(expressionset = "marrayNorm"), function(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          {
            expressionset = try(as(expressionset, "ExpressionSet"))
            if(inherits(expressionset,'try-error'))
              stop("The expressionset is a marrayNorm and cannot be converted automatically in a ExpressionSet. Try to convert it manually.\n")
           arrayQualityMetrics(expressionset, outdir, force, do.logtransform, split.plots, intgroup, grouprep)
          })####end set method marrayNorm

