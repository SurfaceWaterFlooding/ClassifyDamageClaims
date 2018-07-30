# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Claim classification toolbox
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Classify flood damage claims caused by surface water flooding or fluvial
# flooding

# Author:   Daniel Bernet
# Date:     30.07.2018
# License:  GNU General Public License v3.0


import os
import shutil
import time
import arcpy

class Toolbox(object):
    def __init__(self):
        self.label = "Claim classification toolbox"
        self.alias = "classification"

        # List of tool classes associated with this toolbox
        self.tools = [ClassifyClaims]


class ClassifyClaims(object):
    
    # First level functions
    def __init__(self):
        self.label = "ClassifyDamageClaims"
        self.description = "Classify damage claims according to the \
                            triggering processes of surface water flooding \
                            or fluvial flooding."
        self.canRunInBackground = False

    def getParameterInfo(self):
        # Parameter definitions
        
        in_claims = arcpy.Parameter(
            displayName = "Damage claim locations",
            name = "in_claims",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        in_claims.filter.list = ["Point", "Polygon"]

        in_hmap_peri = arcpy.Parameter(
            displayName = "Hazard map perimeters",
            name = "in_hmap_peri",
            datatype = "GPFeatureLayer",
            parameterType = "Optional",
            direction = "Input")
        in_hmap_peri.filter.list = ["Polygon"]

        in_hmap_hzone = arcpy.Parameter(
            displayName = "Hazard map flood zones",
            name = "in_hmap_hzone",
            datatype = "GPFeatureLayer",
            parameterType = "Optional",
            direction = "Input")
        in_hmap_hzone.filter.list = ["Polygon"]

        in_fmap_fzone = arcpy.Parameter(
            displayName="Flood map flood zones",
            name="in_fmap_fzone",
            datatype="GPFeatureLayer",
            parameterType="Optional",
            direction="Input")
        in_fmap_fzone.filter.list = ["Polygon"]

        in_d25 = arcpy.Parameter(
            displayName="Reference distance: d25 (25th percentile) in m",
            name="in_d25",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        in_d25.value = 60.0

        in_d50 = arcpy.Parameter(
            displayName="Reference distance: d50 (50th percentile) in m",
            name="in_d50",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        in_d50.value = 138.0

        in_d75 = arcpy.Parameter(
            displayName="Reference distance: d75 (75th percentile) in m",
            name="in_d75",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        in_d75.value = 305.0

        in_d99 = arcpy.Parameter(
            displayName="Reference distance: d99 (99th percentile) in m",
            name="in_d99",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        in_d99.value = 1148.0

        in_searchdist = arcpy.Parameter(
            displayName="Search distance for closest watercourse",
            name="in_searchdist",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input")
        in_searchdist.value = 5000

        in_dem = arcpy.Parameter(
            displayName = "Digital elevation model",
            name = "in_dem",
            datatype = "GPRasterLayer",
            parameterType = "Required",
            direction = "Input")

        in_rivers = arcpy.Parameter(
            displayName = "River network",
            name = "in_rivers",
            datatype = "GPFeatureLayer",
            parameterType = "Required",
            direction = "Input")
        in_rivers.filter.list = ["Polyline"]

        in_outdir = arcpy.Parameter(
            displayName="Output file directory",
            name="in_outdir",
            datatype="DEFolder",
            parameterType="Required",
            direction="Input")

        in_outshp = arcpy.Parameter(
            displayName = "Name of output feature class",
            name = "in_outshp",
            datatype = "GPString",
            parameterType = "Required",
            direction = "Input")
        in_outshp.value = "ClassifiedClaims"

        in_detailed = arcpy.Parameter(
            displayName = "Detailed (check) or simple (uncheck) output format",
            name = "in_detailed",
            datatype = "GPBoolean",
            parameterType = "Required",
            direction = "Input")
        in_detailed.value = True

        in_cleanup = arcpy.Parameter(
            displayName = "Remove (check) temporary folder and its " + \
                          "content afterwards?",
            name = "in_cleanup",
            datatype = "GPBoolean",
            parameterType = "Required",
            direction = "Input")
        in_cleanup.value = False

        parameters = [in_claims, 
                      in_hmap_peri, in_hmap_hzone, in_fmap_fzone,
                      in_d25, in_d50, in_d75, in_d99,
                      in_searchdist, in_dem, in_rivers,
                      in_outdir, in_outshp, in_detailed, in_cleanup]

        return parameters

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        try:
            if arcpy.CheckExtension("Spatial") != "Available":
                raise Exception
        except Exception:
                return False  # tool cannot be executed

        arcpy.CheckOutExtension("Spatial")
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        if parameters[7].altered or parameters[8].altered:
            if parameters[8].value < 1.5 * parameters[7].value:
                parameters[8].value = 1.5 * parameters[7].value

        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        if not (parameters[4].value < parameters[5].value):
            parameters[4].setErrorMessage("d25 must be smaller than d50")

        if not (parameters[5].value > parameters[4].value):
            parameters[5].setErrorMessage("d50 must larger than d25")

        if not (parameters[5].value < parameters[6].value):
            parameters[5].setErrorMessage("d50 must be smaller than d75")
        
        if not (parameters[6].value > parameters[5].value):
            parameters[6].setErrorMessage("d75 must larger than d50")

        if not(parameters[6].value < parameters[7].value):
            parameters[6].setErrorMessage("d75 must be smaller than d99")

        if not(parameters[7].value > parameters[6].value):
            parameters[7].setErrorMessage("d99 must be larger than d75")

        return

    def initializeParameters(self, parameters, exe_type):
        if exe_type == "Debug":
            # Hardcoded input for debugging and testing
            in_root = os.path.join("D:\\05_Analyses", "Classification", 
                                   "CL_CEucl_1.1_Toolbox")
            in_gis = os.path.join("D:\\08_Data", "GIS")
            in_claims = os.path.join(in_root, "01_Input",
                        "20160911_DDat_PICB_geo_just_maincols_GR_1sample.shp")
            in_hmap_peri = os.path.join(in_gis, "GFK", "v201711", "Derivate", 
                                        "CH_SGK_W_Peri.shp")
            in_hmap_hzone = os.path.join(in_gis, "GFK", "v201711", "Derivate",
                                         "CH_SGK_W_gbr.shp")
            in_fmap_fzone = os.path.join(in_gis, "Aquaprotect", "Derivate", 
                                         "v201711", "A250_kompGFK.shp")
            in_d25 = 60
            in_d50 = 138
            in_d75 = 305
            in_d99 = 1148
            in_searchdist = 2000
            in_outdir = os.path.join(in_root, "02_Toolbox")
            in_dem = os.path.join(in_gis, "swissALTI3D", "swissALTI3D_tif", 
                                  "10m", "swissALTI3D_v2013_10m.tif")
            in_rivers = os.path.join(in_gis, "swissTLM3D", "Derivate",
                                     "1_x_2016", "TLM3D_FGew_n_o.shp")

        elif exe_type == "Dialog":
            # Input from dialog box 
            in_claims = parameters[0].valueAsText
            in_hmap_peri = parameters[1].valueAsText
            in_hmap_hzone = parameters[2].valueAsText
            in_fmap_fzone = parameters[3].valueAsText
            in_d25 = parameters[4].value
            in_d50 = parameters[5].value
            in_d75 = parameters[6].value
            in_d99 = parameters[7].value
            in_searchdist = parameters[8].value
            in_dem = parameters[9].valueAsText
            in_rivers = parameters[10].valueAsText
            in_outdir = parameters[11].valueAsText
            in_outshp = parameters[12].valueAsText
            in_detailed = parameters[13].value
            in_cleanup = parameters[14].value

        parameters = [in_claims, 
                      in_hmap_peri, in_hmap_hzone, in_fmap_fzone,
                      in_d25, in_d50, in_d75, in_d99,
                      in_searchdist, in_dem, in_rivers,
                      in_outdir, in_outshp, in_detailed, in_cleanup]
        arcpy.AddMessage(parameters)
        return parameters

    def initializeWorkspace(self, in_outdir):
        if os.path.exists(in_outdir):
            paths = os.listdir(in_outdir)
            paths = [os.path.join(in_outdir, s) for s in paths]
            for path in paths:
                if os.path.isfile(path):
                    os.remove(path)
                elif os.path.isdir(path):
                    try:
                        shutil.rmtree(path)
                    except WindowsError:
                        time.sleep(5)
                        shutil.rmtree(path)
        else:
            os.mkdir(in_outdir)

        subfolders = ["Temp"]
        subfolders_dirs = [os.path.join(in_outdir, s) for s in subfolders]
        for subfolders_dir in subfolders_dirs:
            if not os.path.exists(subfolders_dir):
                os.mkdir(subfolders_dir)

        arcpy.Delete_management("in_memory")
        arcpy.env.overwriteOutput = True
        # arcpy.env.workspace = "in_memory"
        
    def createOutputDataset(self, in_feature, out_feature):
        arcpy.CopyFeatures_management(in_features = in_feature,
                                      out_feature_class = out_feature)

        fields_drop = []
        fields_orig = arcpy.ListFields(dataset = out_feature)
        for field_orig in fields_orig:
            if field_orig.name in ["FID", "Shape"]:
                next
            else:
                fields_drop.append(field_orig.name)

        if fields_drop != []:
            arcpy.DeleteField_management(in_table = out_feature,
                                         drop_field = fields_drop)

        fields_add = ["fid_orig",   # FID of input claim dataset
                      "z_fromdem",  # elevation of each claim's location
                      "ed_hmap_p",  # Eucl. distance to hazard map perimeter
                      "lg_hmap_p",  # is claim within within hazard map peri.?
                      "ed_hmap",    # Eucl. dist. to hazard map's hazard zone
                      "lg_hmap",    # is claim within hazard map's hazard zone?
                      "ed_fmap",    # Eucl. dist. to flood map's flood zone
                      "lg_fmap",    # is claim within flood map's flood zone?
                      "x_wc_near",  # x-coordinate of closest wc section (AC)
                      "y_wc_near",  # y-coordinate of closest wc section (AC)
                      "z_wc_near",  # z-coordinate of closest wc section (AC)
                      "ed_wc",      # Eucl. dist to closest wc section
                      "aced_wc",    # altitude constr. Eucl. dist (ACED) to wc
                      "dh_wc",      # elevation diff. between claim and wc
                      "aced_perc",  # ACED's percentile relative to ref dist.
                      "class"]      # Claim's class
        field_types = ["LONG",      # fid_orig
                       "FLOAT",     # z_fromdem
                       "FLOAT",     # ed_hmap_p
                       "SHORT",     # lg_hmap_p
                       "FLOAT",     # ed_hmap
                       "SHORT",     # lg_hmap
                       "FLOAT",     # ed_fmap
                       "SHORT",     # lg_fmap
                       "FLOAT",     # x_wc_near
                       "FLOAT",     # y_wc_near
                       "FLOAT",     # z_wc_near
                       "FLOAT",     # ed_wc
                       "FLOAT",     # aced_wc
                       "FLOAT",     # dh_wc
                       "SHORT",     # aced_perc
                       "TEXT"]      # class
        for i in range(0, len(fields_add)):
            # Note: field_precision and field_scale are set, but ignored due to 
            # an unknown reason, thus, they are set to the default values
            arcpy.AddField_management(in_table = out_feature,
                                      field_name = fields_add[i], 
                                      field_type = field_types[i],
                                      field_precision = 3,
                                      field_scale = 15,
                                      field_length = 1)
            if field_types[i] != "TEXT":
                arcpy.CalculateField_management(in_table = out_feature,
                                                field = fields_add[i],
                                                expression = -9999,
                                                expression_type = "PYTHON_9.3")
        
        fid_orig = []
        with arcpy.da.SearchCursor(in_table = in_feature,
                                   field_names = ["FID"]) as c_search_row:
            for search_row in c_search_row:
                fid_orig.append(search_row[0])

        count_i = 0
        with arcpy.da.UpdateCursor(in_table = out_feature,
                                   field_names = ["fid_orig"]) as c_update_row:
            for update_row in c_update_row:
                update_row[0] = fid_orig[count_i]
                c_update_row.updateRow(update_row)
                count_i += 1

        return fid_orig

    def inferZfromPointsAndDEM(self, in_points, in_dem, in_field_name):
        out_points_path = os.path.join("in_memory", "ExtractedPointZ")
        arcpy.sa.ExtractValuesToPoints(in_point_features = in_points, 
                                       in_raster = in_dem,
                                       out_point_features = out_points_path)
        fid_extracted = []
        inferred_z = []
        search_fields = ["fid_orig", "RASTERVALU"]
        with arcpy.da.SearchCursor(in_table = out_points_path,
                                   field_names = search_fields) as c_search_row:
            for search_row in c_search_row:
                fid_extracted.append(search_row[0])
                inferred_z.append(search_row[1])

        update_fields = ["fid_orig", in_field_name]
        with arcpy.da.UpdateCursor(in_table = in_points,
                                   field_names = update_fields) as c_update_row:
            for update_row in c_update_row:
                fid_orig_i = update_row[0]
                if fid_orig_i not in fid_extracted:
                    update_row[1] = -9999.0
                else:
                    ind = fid_extracted.index(fid_orig_i)
                    update_row[1] = inferred_z[ind]
                c_update_row.updateRow(update_row)

    def inferZfromPolygonsAndDEM(self, in_polygons, in_dem, in_field_name):
        out_table_path = os.path.join("in_memory", "ExtractedMeanZ")
        arcpy.sa.ZonalStatisticsAsTable (in_zone_data = in_polygons,
                                         zone_field = "fid_orig",
                                         in_value_raster = in_dem,
                                         out_table = out_table_path,
                                         statistics_type= "MEAN")
        fid_extracted = []
        inferred_z = []
        search_fields = ["fid_orig", "MEAN"]
        with arcpy.da.SearchCursor(in_table = out_table_path,
                                   field_names = search_fields) as c_search_row:
            for search_row in c_search_row:
                fid_extracted.append(search_row[0])
                inferred_z.append(search_row[1])

        update_fields = ["fid_orig", in_field_name]
        with arcpy.da.UpdateCursor(in_table = in_polygons,
                                   field_names = update_fields) as c_update_row:
            for update_row in c_update_row:
                fid_orig_i = update_row[0]
                if fid_orig_i not in fid_extracted:
                    update_row[1] = -9999
                else: 
                    ind = fid_extracted.index(fid_orig_i)
                    update_row[1] = inferred_z[ind]
                c_update_row.updateRow(update_row)  

    def computeED(self, in_feature_from, in_feature_to, in_field_name):
        fields_orig = arcpy.ListFields(in_feature_from)
        field_names_orig = [field.name for field in fields_orig]
        if in_field_name not in field_names_orig:
            arcpy.AddField_management(in_table = in_feature_from,
                                      field_name = in_field_name,
                                      field_type = "DOUBLE")

        if in_feature_to is None: #parameter has not been set
            msg_tmp = "     Assessment of " + in_field_name + \
                      " is skipped and set to a value of 999999" + \
                      " because input feature class was not specified"
            arcpy.AddMessage(msg_tmp)
            print msg_tmp
            arcpy.CalculateField_management(in_table = in_feature_from,
                                            field = in_field_name,
                                            expression = 999999.0,
                                            expression_type = "PYTHON_9.3")
            next

        else:
            arcpy.Near_analysis(in_features = in_feature_from, 
                                near_features = in_feature_to,
                                method = "PLANAR")
            arcpy.CalculateField_management(in_table = in_feature_from,
                                            field = in_field_name,
                                            expression = "!NEAR_DIST!",
                                            expression_type = "PYTHON_9.3")
            arcpy.DeleteField_management(in_table = in_feature_from, 
                                         drop_field = ("NEAR_FID", "NEAR_DIST"))

    def copyDamageLocation(self, in_claims, out_claim, fid, fid_field_name):
        clause = '"%s" = %s' % (fid_field_name, str(fid))
        arcpy.Select_analysis(in_features = in_claims,
                              out_feature_class = out_claim,
                              where_clause = clause)

    def clipRasterDem(self, in_dem, out_dem, xcoord_center, ycoord_center,
                      search_dist):
        xmin = round(xcoord_center - search_dist, 0)
        ymin = round(ycoord_center - search_dist, 0)
        xmax = round(xcoord_center + search_dist, 0)
        ymax = round(ycoord_center + search_dist, 0)
        coord_extent = str(xmin) + " " + \
                       str(ymin) + " " + \
                       str(xmax) + " " + \
                       str(ymax)

        arcpy.Clip_management(in_raster = in_dem,
                              rectangle = coord_extent, 
                              out_raster = out_dem)

    def createFeatureClipMask(self, in_dem, out_mask, zcoord_claim):
        dem_greater_z = arcpy.sa.GreaterThanEqual(in_dem, zcoord_claim)
        dem_mask = arcpy.sa.SetNull(in_conditional_raster = dem_greater_z, 
                                    in_false_raster_or_constant = dem_greater_z,
                                    where_clause = "VALUE = 0")
        arcpy.RasterToPolygon_conversion(in_raster = dem_mask,
        								 out_polygon_features = out_mask,
        								 simplify = "NO_SIMPLIFY")

    def clipFeature(self, in_feature, in_mask, out_feature):
        arcpy.Clip_analysis(in_features = in_feature,
                            clip_features = in_mask,
                            out_feature_class = out_feature)

    def computeACED(self, in_claim, in_constr_feat, search_dist):
        out_near_table = os.path.join("in_memory", "NearTable")
        near_table = arcpy.GenerateNearTable_analysis(
            in_features = in_claim,
            near_features = in_constr_feat,
            out_table = out_near_table,
            location = "LOCATION")
        num_rows = arcpy.GetCount_management(in_rows = near_table)
        num_rows_int = int(num_rows.getOutput(0))

        if num_rows_int == 0: #no features within crop window
            aced = 999999.0
            near_x = -9999.0
            near_y = -9999.0
        else:
            near_field_names = ["NEAR_DIST", "NEAR_X", "NEAR_Y"]
            with arcpy.da.SearchCursor(
                in_table = near_table, 
                field_names = near_field_names) as c_search_row:
                search_row = c_search_row.next()
                aced = search_row[0]
                if aced > search_dist: #features outside of search distance
                    aced = 999999.0
                    near_x = -9999.0
                    near_y = -9999.0
                else:
                    near_x = search_row[1]
                    near_y = search_row[2]

        if aced == 999999.0:
            msg_tmp = "     ATTENTION: No watercourse found within search " +\
                      " radius -> consider increasing search distance!"
            arcpy.AddMessage(msg_tmp)
            print msg_tmp
        
        aced_vals = [aced, near_x, near_y]
        return aced_vals

    def assessPercentileValues(self, in_feature, aced_reference):
        with arcpy.da.UpdateCursor(
            in_table = in_feature,
            field_names = ["aced_wc", "aced_perc"]) as c_update_row:
            for update_row in c_update_row:
                update_row[1] = sum([update_row[0] > x for x in 
                                       aced_reference])
                c_update_row.updateRow(update_row)


    # Second level functions
    def assessEuclideanDists(self, in_feature, in_hmap_peri, in_hmap_hzone,
                             in_fmap_fzone, in_wc, in_dem, d_tol):

        # Compute Euclidean distances to flood maps and watercourses
        in_fmap = [in_hmap_peri, in_hmap_hzone, in_fmap_fzone, in_wc]
        field_names = ["ed_hmap_p", "ed_hmap", "ed_fmap", "ed_wc"]

        for i in range(0, len(in_fmap)):
            msg_tmp = "  -> Processing near calculation: " + field_names[i]
            arcpy.AddMessage(msg_tmp)
            print msg_tmp

            self.computeED(in_feature_from = in_feature,
                           in_feature_to = in_fmap[i],
                           in_field_name = field_names[i])

        # Postprocess Euclidean distance assessment
        fields = ["ed_hmap_p", "lg_hmap_p", "ed_hmap", "lg_hmap", "ed_fmap", 
                  "lg_fmap", "ed_wc"]
        with arcpy.da.UpdateCursor(in_table = in_feature,
                                   field_names = fields) as c_update_row:
            for update_row in c_update_row:
                if update_row[0] != -9999.0:
                    update_row[1] = int(update_row[0] <= d_tol)
                if update_row[2] != -9999.0:
                    update_row[3] = int(update_row[2] <= d_tol)
                if update_row[4] != -9999.0:
                    update_row[5] = int(update_row[4] <= d_tol)
                c_update_row.updateRow(update_row)
        
        # Extract altitude of claim locations based on indicated DEM
        feature_type = arcpy.Describe(in_feature).shapeType
        if feature_type == "Point":
            self.inferZfromPointsAndDEM(in_points = in_feature,
                                        in_dem = in_dem,
                                        in_field_name = "z_fromdem")
        elif feature_type == "Polygon":
            self.inferZfromPolygonsAndDEM(in_polygons = in_feature,
                                          in_dem = in_dem,
                                          in_field_name = "z_fromdem") 

    def assessACEDist2WatCourses(self, in_feature, in_dem, in_rivers, 
                                 claim_fid, search_dist):

        acedist = [-9999.0] * len(claim_fid)
        near_x = [-9999.0] * len(claim_fid)
        near_y = [-9999.0] * len(claim_fid)

        # Compute altidude constrained Euclidean distances to next water course
        for i in range(0, len(claim_fid)):
            msg = "  -> Processing claim i = " + str(i)
            arcpy.AddMessage(msg)
            print msg
            claim_copy = os.path.join("in_memory", "Claim_copy")
            self.copyDamageLocation(in_claims = in_feature,
                                    out_claim = claim_copy, 
                                    fid = claim_fid[i],
                                    fid_field_name = "fid_orig")

            search_fields = ["SHAPE@X", "SHAPE@Y", "z_fromdem"]
            with arcpy.da.SearchCursor(
                    in_table = claim_copy,
                    field_names = search_fields) as c_search_row:
                search_row = c_search_row.next()
                xcoord_i = search_row[0]
                ycoord_i = search_row[1]
                zcoord_i = search_row[2]

            if zcoord_i != -9999.0:
                dem_clip = os.path.join("in_memory", "Dem_clip")
                self.clipRasterDem(in_dem = in_dem, 
                                   out_dem = dem_clip, 
                                   xcoord_center = xcoord_i, 
                                   ycoord_center = ycoord_i, 
                                   search_dist = search_dist)

                dem_mask = os.path.join("in_memory", "Dem_mask")
                self.createFeatureClipMask(in_dem = dem_clip,
                                           out_mask = dem_mask,
                                           zcoord_claim = zcoord_i)

                rivers_clip = os.path.join("in_memory", "Feature_clip")
                self.clipFeature(in_feature = in_rivers,
                                 in_mask = dem_mask,
                                 out_feature = rivers_clip)

                aced_vals_i = self.computeACED(in_claim = claim_copy,
                                               in_constr_feat = rivers_clip,
                                               search_dist = search_dist)
                acedist[i] = aced_vals_i[0]
                near_x[i] = aced_vals_i[1]
                near_y[i] = aced_vals_i[2]
            else:
                next

        acedist_return = [acedist, near_x, near_y]

        return acedist_return

    def postprocessACEDists(self, in_feature, in_dem, in_wc, claim_fid, 
                            aced_results):
        acedist_inferred = aced_results[0]
        x_wc_near_inferred = aced_results[1]
        y_wc_near_inferred = aced_results[2]

        # Assess altitude of closest watercourse section (aced)
        # - When xy location lies on raster boundary, it is not clear which
        #   cell's value is assessed. Thus, assess four locations close to the 
        #   xy coordinates (i.e., 100th of cell resolution NE, SE, SW, NW of xy)  
        #   and take maximum value
        z_wc_near_inferred = [-9999.0] * len(x_wc_near_inferred)
        dem_res_x = arcpy.GetRasterProperties_management(
                                                in_raster = in_dem,
                                                property_type = "CELLSIZEX")
        dem_res_y = arcpy.GetRasterProperties_management(
                                                in_raster = in_dem,
                                                property_type = "CELLSIZEY")
        dem_shift_x = float(dem_res_x.getOutput(0)) / 100
        dem_shift_y = float(dem_res_y.getOutput(0)) / 100

        for i in range(0, len(x_wc_near_inferred)):
            xy_NE = str(x_wc_near_inferred[i] + dem_shift_x) + " " + \
                    str(y_wc_near_inferred[i] + dem_shift_y)
            z_NE = arcpy.GetCellValue_management(in_raster = in_dem, 
                                                 location_point = xy_NE)
            xy_SE = str(x_wc_near_inferred[i] + dem_shift_x) + " " + \
                    str(y_wc_near_inferred[i] - dem_shift_y)
            z_SE = arcpy.GetCellValue_management(in_raster = in_dem, 
                                                 location_point = xy_SE)
            xy_SW = str(x_wc_near_inferred[i] - dem_shift_x) + " " + \
                    str(y_wc_near_inferred[i] - dem_shift_y)
            z_SW = arcpy.GetCellValue_management(in_raster = in_dem, 
                                                 location_point = xy_SW)
            xy_NW = str(x_wc_near_inferred[i] - dem_shift_x) + " " + \
                    str(y_wc_near_inferred[i] + dem_shift_y)
            z_NW = arcpy.GetCellValue_management(in_raster = in_dem, 
                                                 location_point = xy_NW)
            z_cellvalues = [z_NE.getOutput(0), z_SE.getOutput(0), 
                            z_SW.getOutput(0), z_NW.getOutput(0)]
            z_cellnum = []
            for z_cellvalue in z_cellvalues:
                try:
                    z_num = float(z_cellvalue)
                    z_cellnum.append(z_num)
                except:
                    next
            if z_cellnum != []:
                z_wc_near_inferred[i] = max(z_cellnum)
            else: 
                next

        # Transfer values to in_feature
        fields = ["fid_orig", "z_fromdem", "x_wc_near", "y_wc_near", 
                  "z_wc_near", "aced_wc", "dh_wc"]
        with arcpy.da.UpdateCursor(in_table = in_feature,
                                   field_names = fields) as c_update_row:
            for update_row in c_update_row:
                fid_orig_i = update_row[0]
                ind_i = claim_fid.index(fid_orig_i)
                update_row[2] = x_wc_near_inferred[ind_i]
                update_row[3] = y_wc_near_inferred[ind_i]
                update_row[4] = z_wc_near_inferred[ind_i]
                update_row[5] = acedist_inferred[ind_i]
                if update_row[1] != -9999.0 and update_row[4] != -9999.0:
                    update_row[6] = update_row[4] - update_row[1]
                c_update_row.updateRow(update_row)

    def classifyClaims(self, in_feature, aced_reference):

        self.assessPercentileValues(in_feature = in_feature, 
                                    aced_reference = aced_reference)
        
        update_fields = ["lg_hmap_p", "lg_hmap", "lg_fmap", "aced_perc",
                         "class"]
        with arcpy.da.UpdateCursor(
            in_table = in_feature,
            field_names = update_fields) as c_update_row:
            for update_row in c_update_row:
                lg_hmap_p = update_row[0]
                lg_hmap = update_row[1]
                lg_fmap = update_row[2]
                aced_perc = update_row[3]

                # Class A: Most likely surface water flood
                lg_A_1a = lg_hmap_p & (not lg_hmap) & (aced_perc in [3, 4])
                lg_A_1b = (not lg_hmap_p) & (not lg_fmap) & (aced_perc == 4)

                # Class B: Likely surface water flood
                lg_B_2a = lg_hmap_p & (not lg_hmap) & (aced_perc == 2)
                lg_B_2b = (not lg_hmap_p) & (not lg_fmap) & (aced_perc == 3)

                # Class C: Surface water flood or fluvial flood
                lg_C_3a = lg_hmap_p & (not lg_hmap) & (aced_perc == 1)
                lg_C_3b = (not lg_hmap_p) & (not lg_fmap) & (aced_perc == 2)

                # Class D: Likely fluvial flood
                lg_D_4a = lg_hmap_p & (not lg_hmap) & (aced_perc == 0)
                lg_D_4b = (not lg_hmap_p) & (not lg_fmap) & (aced_perc == 1)
                lg_D_4c = (not lg_hmap_p) & lg_fmap

                # Class E: Most likely fluvial flood
                lg_E_5a = lg_hmap_p & lg_hmap
                lg_E_5b = (not lg_hmap_p) & (not lg_fmap) & (aced_perc == 0)

                claim_class = "NULL"
                if lg_A_1a or lg_A_1b:
                    claim_class = "A"
                elif lg_B_2a or lg_B_2b:
                    claim_class = "B"
                elif lg_C_3a or lg_C_3b:
                    claim_class = "C"
                elif lg_D_4a or lg_D_4b or lg_D_4c:
                    claim_class = "D"
                elif lg_E_5a or lg_E_5b:
                    claim_class = "E"
                else:
                    claim_class = "NA"

                update_row[4] = claim_class
                c_update_row.updateRow(update_row)

    # Top level function
    def execute(self, parameters, messages):
        is_licensed = self.isLicensed()
        # print(is_licensed)    

        msg_tmp = "Initializing parameters and preparing output"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp
        # "Debug" (run in VSC) or "Dialog" (run via ArcGIS Toolbox)
        parameters = self.initializeParameters(parameters, "Dialog")
        claims_input = parameters[0]
        hmap_peri = parameters[1]
        hmap_hzone = parameters[2]
        fmap_fzone = parameters[3]
        d25 = float(parameters[4])
        d50 = float(parameters[5])
        d75 = float(parameters[6])
        d99 = float(parameters[7])
        searchdist = float(parameters[8])
        dem = parameters[9]
        rivers = parameters[10]
        outdir = parameters[11]
        outshp = parameters[12]
        detailed = parameters[13]
        cleanup = parameters[14]

        self.initializeWorkspace(in_outdir = outdir)

        claims_temp = os.path.join("in_memory", "Claims_temp")
        claim_fid_orig = self.createOutputDataset(in_feature = claims_input,
                                                  out_feature = claims_temp)
        result_backup_00 = os.path.join(outdir, "Temp", "00_Output_Init.shp")
        arcpy.CopyFeatures_management(in_features = claims_temp,
                                      out_feature_class = result_backup_00)

        msg_tmp = "Assessing Euclidean distances to flood maps"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp
        self.assessEuclideanDists(in_feature = claims_temp, 
                                  in_hmap_peri = hmap_peri, 
                                  in_hmap_hzone = hmap_hzone,
                                  in_fmap_fzone = fmap_fzone,
                                  in_wc = rivers,
                                  in_dem = dem, 
                                  d_tol = 25)
        result_backup_01 = os.path.join(outdir, "Temp", "01_Output_EDist.shp")
        arcpy.CopyFeatures_management(in_features = claims_temp,
                                      out_feature_class = result_backup_01)

        msg_tmp = "Assessing alt. constr. eucl. dist. to watercourses"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp
        acedist_vals = self.assessACEDist2WatCourses(in_feature = claims_temp,
                                                     in_dem = dem,
                                                     in_rivers = rivers,
                                                     claim_fid = claim_fid_orig,
                                                     search_dist = searchdist)

        self.postprocessACEDists(in_feature = claims_temp,
                                 in_dem = dem,
                                 in_wc = rivers,
                                 claim_fid = claim_fid_orig,
                                 aced_results = acedist_vals)
        result_backup_02 = os.path.join(outdir, "Temp", "02_Output_ACEDist.shp")
        arcpy.CopyFeatures_management(in_features = claims_temp,
                                      out_feature_class = result_backup_02)

        msg_tmp = "Classify damage claims"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp

        aced_ref = [d25, d50, d75, d99]
        self.classifyClaims(in_feature = claims_temp, 
                            aced_reference = aced_ref)
        result_backup_03 = os.path.join(outdir, "Temp", "03_Output_Cl.shp")
        arcpy.CopyFeatures_management(in_features = claims_temp,
                                      out_feature_class = result_backup_03)

        msg_tmp = "Export classified damage claims"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp

        export_claims = os.path.join("in_memory", "Export_claims")
        arcpy.CopyFeatures_management(in_features = claims_input,
                                      out_feature_class = export_claims)
        arcpy.AddField_management(in_table = export_claims,
                                  field_name = "fid_orig", 
                                      field_type = "LONG")
        count_i = 0
        with arcpy.da.UpdateCursor(in_table = export_claims,
                                   field_names = ["fid_orig"]) as c_update_row:
            for update_row in c_update_row:
                update_row[0] = claim_fid_orig[count_i]
                c_update_row.updateRow(update_row)
                count_i += 1

        fields_data = arcpy.ListFields(dataset = claims_temp)
        fields_add = []
        if detailed:
            for field_data in fields_data:
                if field_data.name in ["FID", "Shape", "fid_orig"]:
                    next
                else:
                    fields_add.append(field_data.name)
        else:
            fields_add.append(["class"])

        arcpy.JoinField_management(in_data = export_claims,
                                   in_field = "fid_orig", 
                                   join_table = result_backup_03,
                                   join_field = "fid_orig",
                                   fields = fields_add)

        if outshp[-4] != ".shp":
            outshp = outshp + ".shp"
        final_output = os.path.join(outdir, outshp)
        arcpy.CopyFeatures_management(in_features = export_claims,
                                      out_feature_class = final_output)

        msg_tmp = "Clean-up"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp

        if cleanup:
            folder_temp = os.path.join(outdir, "Temp")
            if os.path.exists(folder_temp):
                paths = os.listdir(folder_temp)
                paths = [os.path.join(folder_temp, s) for s in paths]
                for path in paths:
                    if os.path.isfile(path):
                        os.remove(path)
                    elif os.path.isdir(path):
                        try:
                            shutil.rmtree(path)
                        except WindowsError:
                            time.sleep(5)
                            shutil.rmtree(path)
                os.remove(folder_temp)

        msg_tmp = "Classification process completed"
        arcpy.AddMessage(msg_tmp)
        print msg_tmp


# def main():
#     tbx = Toolbox()
#     tool = ClassifyClaims()
#     parameters = tool.getParameterInfo()
#     tool.execute(parameters = parameters, messages = None)


# if __name__ == "__main__":
#     main()
