# ClassifyDamageClaims
Identify flood damage claims caused by surface water flooding or fluvial flooding with an ArcGIS Python Toolbox 

## Background
Surface water flooding is a considerable threat for our society. However, the process is difficult to study because very often suitable data is lacking. An interesting source of information can be flood claims of insurance companies. For instance, flood claims to buildings can point to the dates and locations where flooding must have occured. However, as insurance companies rarely differentiate between claims caused be fluvial flooding or surface water flooding (SWF), the claims need to be classified first. For that purpose, [Bernet et al. (2017)](https://doi.org/10.5194/nhess-17-1659-2017) developped a classification scheme. The scheme has then been implemented in an ArcGIS Python Toolbox so that it can easily be used by others. A detailed description of the method can be found in [Bernet et al. (2017)](https://doi.org/10.5194/nhess-17-1659-2017).

Note that the method has been developed for claims in Switzerland. When using it elsewhere, the method might need to be adapted accordingly.

## Prerequisites
### Software
This toolbox was developed and tested with *ArcGIS 10.6*, which is shipped with *Python 2.7*. In addition, the toolbox requires the *Spatial Analyst* extention. The toolbox may also work with other versions, however, this has not been tested. Refer to the following links for detailed information regarding these software packages, their dependencies and how to install them:

* [Arcpy](http://desktop.arcgis.com/en/arcmap/latest/analyze/arcpy/what-is-arcpy-.htm)
* [ArcGIS Spatial Analyst extension](http://desktop.arcgis.com/en/arcmap/latest/extensions/spatial-analyst/what-is-the-spatial-analyst-extension.htm)
* [Python](https://www.python.org/)

### Data
* Damage claim locations (point or polygon feature class)
* Flood hazard map and/or additional flood maps (polygon feature class)
* Digital elevation model (raster)
* River network (line feature class)

### Getting started
* Download the repository to a specific, known folder on your local computer. 
* Start **ArcMap** or **ArcCatalog**
* Using the **Catalog Tree** navigate to the folder, where the repository is stored. If the folder is not visible, a folder connection might first need to be added.
* Unfold the toolbox **DamageClassification.pyt** and start the tool's dialog window by double clicking on **Classify damage claims**
* Follow the instructions and start the tool be clicking on **OK**

## References
Bernet, D. B., Prasuhn, V., and Weingartner, R.: Surface water floods in Switzerland: what insurance claim records tell us about the damage in space and time, Nat. Hazards Earth Syst. Sci., 17, 1659-1682, doi: [10.5194/nhess-17-1659-2017](https://doi.org/10.5194/nhess-17-1659-2017), 2017.
