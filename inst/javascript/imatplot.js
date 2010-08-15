// global variables
var messageText;

function init(evt) {
  // get reference to child (content) of the text-Element
  messageText = document.getElementById("annotationtext").firstChild;
}


function toggleSeries(which, on)
{
   var el = document.getElementById(which);
   if(!el) {
      return(false);
   }

   var cur = el.getAttribute('stroke-width');
   var val;
   var factor = 4;

   if(on) {
      val = cur * factor; 
      messageText.nodeValue = which;

   } else {
      val = cur / factor;
      messageText.nodeValue = " ";
   }

   el.setAttribute('stroke-width', val);
   
   return(true);
}
