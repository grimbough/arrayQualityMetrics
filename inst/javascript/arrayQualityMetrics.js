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
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find SVG objects and cache their locations------*/
    svgObjects = new Array(svgObjectNames.length);
    for(i=0; i<svgObjects.length; i++) 
    {
        svgObjects[i] = safeGetElementById("Fig:"+svgObjectNames[i]);
    }

    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    // checkboxes[a] is (expected to be) of class HTMLInputElement
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        setAllPlots("r:"+(a+1), checkboxes[a].checked);
    }
 
}


safeGetElementById = function(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
 This function iterates over all SVG objects and updates the 
 plot elements corresponding to report object 'ro'. 
 Arguments:
 reportObjId: character with the id of the report object to be updated
 status: logical, indicating whether to highlight
 ---------------------------------------------------------------*/
function setAllPlots(reportObjId, status)
{
    var i, idx_status, el, id;

    idx_status = (0+status); // convert from logical (FALSE, TRUE) to integer (0, 1)

    for(i=0; i<svgObjects.length; i++) 
    {
	id = idFuns[i][1](reportObjId);
        for(j=0; j<id.length; j++)
	{
	    el = svgObjects[i].contentDocument.getElementById(id[j]);
	    if(!el) 
	    { 
		throw new Error("Did not find Id '" + id[j] + "'");
	    }
	    el.setAttribute('stroke-opacity', strokeOpacity[i][idx_status]); 
	    el.setAttribute('stroke-width',   strokeWidth[i][idx_status]); 
	} // for j

        showTipTable(i, reportObjId);
    } // for i

}


/*------------------------------------------------------------
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    if(rows.length != arrayMetadata[array].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[array].length);

    var a = parseInt(reportObjId.replace('^r:', ''));
    if(isNan(a)) throw new Error('Invalid report object id ' + x);

    for(i=0; i<rows.length; i++) 
    {
	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
    }
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
    {
	rows[i].cells[1].innerHTML = "";
    }
}


/*------------------------------------------------------------
  From 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function geti(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
    {
        if(svgObjectNames[i] == name)
	{
	    return i;
	}
    }
    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  ------------------------------------------------------------*/
function showTip(plotObjId, name) 
{   
    var i = geti(name);
    var reportObjId = idFuns[i][2](plotObjId);
    showTipTable(i, reportObjId);
}

function hideTip(plotObjId, name)
{
    hideTipTable(geti(name));
}

