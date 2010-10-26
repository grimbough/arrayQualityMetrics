setClass("aqmReportModule",
   representation(
     plot     = "ANY",
     section  = "character",
     title    = "character",
     legend   = "character",
     shape    = "list",
     outliers = "integer",
     svg      = "list"),
    prototype(
     plot      = new("namedList"),
     section   = NA_character_,
     title     = NA_character_,
     legend    = NA_character_,
     shape     = list(),
     outliers  = integer(0),       
     svg       = list()))

