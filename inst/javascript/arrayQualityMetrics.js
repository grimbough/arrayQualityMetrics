// (C) Wolfgang Huber 14.11.2010

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

function reportinit() 
{
 
    var a, i, status;

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

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setAllPlotObjsInAllPlots("r:"+(a+1), status);
    }
}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Plot Objects 
 ---------------------------------------------------------------*/
function setAllPlotObjsInAllPlots(reportObjId, status)
{
    var i, j, plotObjIds;

    for(i=0; i<svgObjects.length; i++) 
    {
	plotObjIds = idFuns[i](reportObjId);
        for(j=0; j<plotObjIds.length; j++) 
	    setOnePlotObjInOnePlot(i, plotObjIds[j], status)

        showTipTable(i, reportObjId);
    } 
}

function setOnePlotObjInOnePlot(i, plotObjId, status)
{
    var el = svgObjects[i].contentDocument.getElementById(plotObjId);
    if(!el) 
	throw new Error("Did not find Id '" + plotObjId + "'");
  
    // '0+' converts integer to boolean
    el.setAttribute('stroke-opacity', strokeOpacity[i][0+status]); 
    el.setAttribute('stroke-width',   strokeWidth[i][0+status]); 
}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/

function getIndexFromReportObjId(reportObjId)
{
   var a = parseInt(reportObjId.replace(/^r:/, ''));
   if(isNaN(a)) 
       throw new Error('Invalid report object id ' + reportObjId);
   return (a-1);
}

function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = getIndexFromReportObjId(reportObjId);

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, plotObjId, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = getIndexFromReportObjId(reportObjId)
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setAllPlotObjsInAllPlots(reportObjId, status);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = getIndexFromReportObjId(reportObjId);
    var status = checkboxes[a].checked;
    setAllPlotObjsInAllPlots(reportObjId, status);
}

