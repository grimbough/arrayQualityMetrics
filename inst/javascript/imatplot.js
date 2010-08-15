// global variables
var messageText;

function init(evt) {
  // get reference to child (content) of the text-Element
  messageText = document.getElementById("annotationtext").firstChild;
}


function toggleSeries(which, title, on)
{

    var el;
    var cur;
    var val;
    var factor = 4;

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

	cur = el.getAttribute('stroke-width');
	
	if(on) {
	    val = cur * factor; 
	} else {
	    val = cur / factor;
	}
	el.setAttribute('stroke-width', val);
    }
   return(true);
}
