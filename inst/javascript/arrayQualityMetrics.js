// (C) Wolfgang Huber 13.11.2010

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.
var highlightInitial = [ @HIGHLIGHTINITIAL@ ];
var arrayMetadata    = [ @ARRAYMETADATA@ ];
var svgObjectNames   = [ @SVGOBJECTNAMES@ ];
var idFuns           = [ @IDFUNS@ ];
var strokeOpacity    = [ @STROKEOPACITY@ ];
var strokeWidth      = [ @STROKEWIDTH@ ];

// Global variables - these are set up below by 'reportinit'

var svgObjects;         // array of all the SVG objects on the page
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes

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
        svgObjects[i] = safeGetElementById("Fig:"+svgObjectNames[i]);
    }

    /*--------find associated tables and cache their locations------*/
    tables = new Array(tableIds.length);
    for(i=0; i<tableIds.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    // checkboxes[a] is (expected to be) of class HTMLInputElement
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightArraysInitial[a];
        setAllPlots(a, checkboxes[a].checked);
    }
 
}


safeGetElementById = function(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

// Callback for 'onchange' events of the checkboxes.
// Arguments:
// array:  numeric (integer) index of the array to be updated
function checkboxEvent(array)
{
    var status;
    status = checkboxes[array].checked;
    setAllPlots(array, status);
}

// Callback for 'onclick' events of the plot elements
// Arguments:
// array:  numeric (integer) index of the array to be updated
function clickPlotElement(array)
{
    var status;
    status = !checkboxes[array].checked;
    checkboxes[array].checked = status;
    setAllPlots(array, status);
}

// This function iterates over all SVG objects and updates the plot element
// with id 'aqm'+array. It also attempts to update the plots elements
// with ids 'aqm'+(array+n), 'aqm'+(array+2*n), ..., this is used e.g.
// in the density plots for two-colour arrays, where we have three lines
// per array. 

// Arguments:
// array: numeric (integer) index of the array to be updated
// status: logical, indicating whether to highlight
function setAllPlots(array, status)
{
    var i, idx_status;
    var el;
    var id;

    idx_status = (0+status); // convert from logical (FALSE, TRUE) to integer (0, 1)

    for(i=0; i<svgObjects.length; i++) 
    {
	id = "aqm" + array;
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


// From table ID (e.g. 'Tab:density'), find its numeric index in the 'tables' array.
function geti(tabID) 
{
    var i;
    for(i=0; i<tableIds.length; i++)
    {
        if(tableIds[i] == tabID)
	{
	    return(i);
	}
    }
    throw new Error("Did not find '" + tabID + "'.");
}

// From table ID (e.g. 'Tab:density'), find the table rows
function getrows(tabID)
{
    var rows = tables[geti(tabID)].rows;
    return(rows);
}

function showTip(tabID, array) {   

    rows = getrows(tabID);
    if(rows.length != arrayMetadata[array].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[array].length);

    for(i=0; i<rows.length; i++) 
    {
	rows[i].cells[1].innerHTML = arrayMetadata[array][i];
    }
}

function hideTip(tabID)
{
    rows = getrows(tabID); 
    for(i=0; i<rows.length; i++) 
    {
	rows[i].cells[1].innerHTML = "";
    }
}

