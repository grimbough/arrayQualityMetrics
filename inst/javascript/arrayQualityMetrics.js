
// (C) Wolfgang Huber 1.11.2010

// script parameters - these are set by 'writeReport' when copying the 
//   template from arrayQualityMetrics/inst/scripts into the report.
var svgObjectIds = [ @SVGOBJECTIDS@ ];
var highlightArraysInitial = [ @HIGHLIGHTARRAYSINITIAL@ ];
var strokeOpacity = [ @STROKEOPACITY@ ];
var strokeWidth   = [ @STROKEWIDTH@ ];

// var svgObjectIds = ["svg1", "svg2"];
// var highlightArraysInitial = [ false, true ];
// var strokeOpacity = [ [0.4, 1], [0.4, 1] ];
// var strokeWidth   = [ [1, 3], [1, 6] ];
// var strokeColor   = ["rgb(0%,0%,0%)", "rgb(100%,0%,0%)"];


var latestArray;        // info about the most recently selected array
var svgObjects;         // array of all the SVG objects on the page
var checkboxes;         // location of the checkboxes
var tipObject;          // tooltips


function reportinit() {
 
    var a, i;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ArrayCheckBoxes");
    if(checkboxes.length != highlightArraysInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightArraysInitial.length="+ highlightArraysInitial.length);
    
    /*--------find SVG objects and cache their locations------*/
    svgObjects = new Array(svgObjectIds.length);
    for(i=0; i<svgObjectIds.length; i++) 
    {
        svgObjects[i] = document.getElementById(svgObjectIds[i]);
        if(svgObjects[i]==null)
            throw new Error("Id "+ svgObjectIds[i] + " not found.");
    }

    // checkboxes[a] is (expected to be) of class HTMLInputElement
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightArraysInitial[a];
        setAllPlots(a, checkboxes[a].checked);
    }

    tipObject = document.getElementById("arraytooltip");
    if(tipObject==null)
        throw new Error("Id 'arraytooltip' not found.");
 
}

// array - numeric (integer) index of the array to be updated
function checkboxEvent(array)
{
    var status;
    status = checkboxes[array].checked;
    setAllPlots(array, status);
}

// array - numeric (integer) index of the array to be updated
function setAllPlots(array, status)
{
    var i, idx_status;
    var el;
    var id;

    idx_status = (0+status); // convert from logical (FALSE, TRUE) to integer (0, 1)

    for(i=0; i<svgObjectIds.length; i++) 
    {
	id = "aqm_" + (array+1);
	el = svgObjects[i].contentDocument.getElementById(id);
	if(!el) 
	{ 
            throw new Error("Did not find Id '" + id + "'");
	}
	// el.setAttribute('stroke',         strokeColor[idx_status]); 
	el.setAttribute('stroke-opacity', strokeOpacity[i][idx_status]); 
	el.setAttribute('stroke-width',   strokeWidth[i][idx_status]); 

    }
}

function clickPlotElement(array)
{
    var status;
    status = !checkboxes[array].checked;
    checkboxes[array].checked = status;
    setAllPlots(array, status);
}



function showTip(array) 
{
    var curX =100;
    var curY = 100;

    var offsetxpoint = -60; // Customize x offset of tooltip
    var offsetypoint =  20;  // Customize y offset of tooltip
 
    thetext = "blabla";
    tipObject.innerHTML = thetext;
    tipObject.style.left = curX + offsetxpoint + "px";
    tipObject.style.top = curY + offsetypoint+"px"
    tipObject.style.visibility = "visible"

    return false
}

function hideTip()
{
    tipObject.style.visibility = "hidden";
}





