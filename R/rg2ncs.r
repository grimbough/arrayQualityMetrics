rg2ncs = function(RG, type = "complete")
  {
    if(type == "complete")
       ad = with(RG, assayDataNew(R=R, G=G, Rb=Rb, Gb=Gb))
    if(type == "nobkg")
      ad = with(RG, assayDataNew(R=R, G=G))

    ncs = new("NChannelSet", assayData = ad)
  }

