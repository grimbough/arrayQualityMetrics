

// (C) WH, 16.8.2010, from an example by DTL

// Global variable: content of a message window for displaying object names.
var messageText;

function init(evt) {
  // Get reference to child (content) of the text-Element
  messageText = document.getElementById("annotationtext").firstChild;
}

// This function is called upon 'onmouseover' (on=TRUE) and 
// 'onmouseout' (on=FALSE) events. 'which' is a vector of IDs 
// of the plot elements to be toggled, 'title' the text to be 
// displayed in the message text window in the case of 'onmouseover'.
function toggleSeries(which, title, on)
{

    var el;
    var oldwidth;
    var newwidth;
    var factor = 5;

    if (on) {
	messageText.nodeValue = title;
    } else {
	messageText.nodeValue = " ";
    }

    for( var i = 0; i < which.length; i++ ) {

	el = document.getElementById(which[i]);
	if(!el) { 
	    throw new Error("Did not find 'which[i]'.");
	}

	oldwidth = el.getAttribute('stroke-width');
	if(on) {
	    newwidth = oldwidth * factor; 
            if(!el.getAttribute('original-stroke-opacity'))
              el.setAttribute('original-stroke-opacity', el.getAttribute('stroke-opacity'));
	    el.setAttribute('stroke-opacity', 1);
	} else {
	    newwidth = oldwidth / factor;
	    el.setAttribute('stroke-opacity', el.getAttribute('original-stroke-opacity'));
	}
	el.setAttribute('stroke-width', newwidth);


    }
   return(true);
}


