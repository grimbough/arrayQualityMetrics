rg2ncs = function(RG, type = "complete")
  {
    if(type == "complete")
       ad = with(RG, assayDataNew(R=R, G=G, Rb=Rb, Gb=Gb))
    if(type == "nobkg")
      ad = with(RG, assayDataNew(R=R, G=G))
    
    ncs = new("NChannelSet", assayData = ad)
    
    if(length(RG$genes) != 0)
      {
        fd = new("AnnotatedDataFrame", data = RG$genes)
        featureData(ncs) = try(fd)

      }
    if(length(RG$genes) != 0)
      {
        pd = new("AnnotatedDataFrame", data = RG$targets)
        phenoData(ncs) = try(pd)
      }

  }

