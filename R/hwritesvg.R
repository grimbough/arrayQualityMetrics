# from hwriter

aqm.hwriteImage = function (image.url, page = NULL, ..., image.border = 0, width = NULL, 
    height = NULL, capture = FALSE, tagname = "img") 
{
    if (capture) {
        if (is.null(width)) 
            width = 400
        if (is.null(height)) 
            height = 400
        dev.print(png, width = width, height = height, image.url)
    }
    str = hwriter::hmakeTag(tagname, border = image.border, src = image.url, 
        alt = image.url, width = width, height = height)
    hwriter::hwrite(str, page, ...)
}

