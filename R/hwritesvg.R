## I am not yet sure what is the best way to include SVG
## into HTML pages. One problem is lack of support from IE <= 8
##  (but it is announced for IE 9). Also, there seem to be at
## least four different tags that can be used
## <embed>, <object>, <iframe>, <img> (the latter is currently 
## reported to be only supported by Opera, Safari and Google Chrome).
##
## Here, based on http://www.carto.net/papers/svg/samples/svg_html.shtml
## I choose <object>
##
## The function tries to be clever and guesses the image format from the extension.

aqm.hwriteImage = function (image.url, page = NULL, ..., image.border = 0, width = NULL, 
    height = NULL, capture = FALSE) 
{
    if (capture) {
        if (is.null(width)) 
            width = 400
        if (is.null(height)) 
            height = 400
        dev.print(png, width = width, height = height, image.url)
    }

    getExtension = function(x) { s = strsplit(x, split=".", fixed=TRUE)[[1]] ; s[length(s)] }
    alt = paste(image.url, "appears to be missing or not renderable by this browser.")
    
    str = switch(getExtension(image.url),
      "svg" = {
        alt = paste(alt, "To render SVG, consider using Firefox, Safari, Opera or IE >=9.")
        ## This is ugly: in order to avoid confusion with hmakeTag's own 'data' argument
        ##   I use DATA, using the fact that R is case-sensitive while HTML is not.
        ## I also put the 'alt' text both in the alt attribute and as text between <object..> and </object>,
        ##   as IE does not seem to honour the alt attribute.
        hwriter::hmakeTag("object", type="image/svg+xml", DATA = image.url, border = image.border, 
                          alt = alt, data = alt, width = width, height = height)
      }, {
      ## default:
        hwriter::hmakeTag("img", src = image.url, border = image.border, 
                          alt = alt, width = width, height = height)
      })
    hwriter::hwrite(str, page, ...)
}

