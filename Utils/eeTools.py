import shapefile
import numpy as np
import ee
import functools
import tagee
import math
    
# def feature2ee(file):
#     gdf = shapefile.Reader(file)    
#     geom_dict = shapefile.Reader(file).__geo_interface__

#     features=[]
#     g = [i for i in gdf.shapes()]
#     for i in range(len(g)):
#         x,y = g[i].points[0]
#         cords = np.dstack((x,y)).tolist()
#         double_list = functools.reduce(lambda x,y: x+y, cords)
#         single_list = functools.reduce(lambda x,y: x+y, double_list)                    
#         geom_=ee.Geometry.Point(single_list)
#         feature = ee.Feature(geom_)
#         features.append(feature)
#     return ee.FeatureCollection(features) 



def getCovariates(ftrs,start,end,resol):
    
    EVAPOTRANSPIRATION = ee.ImageCollection("MODIS/006/MOD16A2")
    VEGETATION = ee.ImageCollection("MODIS/006/MOD13A1")
    MODIS = ee.ImageCollection("MODIS/061/MOD11A1")
    DEM = ee.Image("USGS/SRTMGL1_003")
    DEMAttributes = tagee.TAGEE.terrainAnalysis(DEM, ftrs.geometry().bounds())   
    landSurfaceTemperature = MODIS.select('LST_Day_1km')
    
    
    Locations = ftrs.map(lambda ftr: ftr.set('system:time_start',ee.Date('2013-01-01').millis()))
    
    

    ##### Temperature Fluctutation ####


    #Add Land surface temperature dataset//

    LST = landSurfaceTemperature\
          .filterDate(start,end)\
          .map(lambda img: img.clip(Locations.geometry().bounds()))\
          .select('LST_Day_1km')

    #Calculate Coldest day
    def getMinT(image):
        
        value_lst = ee.Image(image)\
        .reduceRegion(reducer=ee.Reducer.mean(),geometry=ee.Image(image).geometry(),bestEffort=True)\
        .get('LST_Day_1km')
    
        date_ = ee.Image(image)\
                .get('system:time_start')

        return ee.Feature(ee.Image(image).geometry()).set('LST_Day_1km',value_lst).set('system:time_start',date_)



    LST_list = LST.map(getMinT)
    LST_list = LST_list.filter(ee.Filter.gt('LST_Day_1km',0)) #some days have null values  
    MinT=ee.Date(LST_list.sort('LST_Day_1km',True).first().get('system:time_start'))
    


     ##neighbourhood statistics

    VEGETATION_VARIABILITY = VEGETATION.map(lambda img: img.reduceNeighborhood\
                                            (reducer= ee.Reducer.stdDev(),\
                                              kernel= ee.Kernel.circle(50))\
                                           )
    #Harmonic model  for Land Surface Temperature

    independents = ee.List(['constant', 't'])
    dependent = ee.String('LST_Day_1km')
    
    # for a harmonic model, let's create a function that adds Temp, time, cos(time), 
    # sin(time) and constant

    def addInfo(img):
        Temp = ee.Image(img)
        date0 = MinT
        time = ee.Image(ee.Date(img.get('system:time_start')).difference(date0, 'year')).float()
        timeRadians = time.multiply(2 * math.pi)
        cosT = timeRadians.cos()
        sinT = timeRadians.sin()
        return img.addBands(time.rename('t')).addBands(ee.Image(1).rename('constant'))\
                  .addBands(Temp.rename('Temp')).addBands(cosT.rename('cos')).addBands(sinT.rename('sin'))\
                  .copyProperties(img, ['system:time_start'])

    landSurfaceTemperature_Harmonic_analysis = LST.map(addInfo)

    #now our independents are 4 terms
    harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin'])
    #The output of the regression reduction is a 4x1 array image.
    LST_harmonicTrend = landSurfaceTemperature_Harmonic_analysis\
                        .select(harmonicIndependents.add(dependent))\
                        .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1))

     #Turn the array image into a multi-band image of coefficients.
    LST_harmonicTrendCoefficients = LST_harmonicTrend.select('coefficients')\
                                                     .arrayProject([0])\
                                                     .arrayFlatten([harmonicIndependents])

    #Compute fitted values.
    LST_fittedHarmonic = landSurfaceTemperature.map(lambda img: img.addBands(\
                                                                img.select(harmonicIndependents)\
                                                                   .multiply(LST_harmonicTrendCoefficients)\
                                                                   .reduce('sum')\
                                                                   .rename('fitted')))



    #From these, we can estimate the phase and amplitude of the collection
    LST_phase_Temp = LST_harmonicTrendCoefficients.select('cos').atan2(LST_harmonicTrendCoefficients.select('sin'))
            
    LST_amplitude_Temp = LST_harmonicTrendCoefficients.select('cos').hypot(LST_harmonicTrendCoefficients.select('sin'))




    
    ##Create and print the chart for a selected area.
    ##add Terrain dataset 

    aspect = ee.Terrain.aspect(DEM)
    hillshade = ee.Terrain.hillshade(DEM)
    slope = ee.Terrain.slope(DEM)

    #add secondary DEM attributes

    Hor_curv = DEMAttributes.select('HorizontalCurvature')
    Ver_curv = DEMAttributes.select('VerticalCurvature')


    ###add Vegetation dataset 


    ET = EVAPOTRANSPIRATION\
         .filterDate(start,end)\
         .map(lambda img: img.clip(ftrs.geometry().bounds()))\
         .select('ET')

    VEGETATION = VEGETATION\
         .filterDate(start,end)\
         .filterBounds(ftrs)\
         .map(lambda img: img.clip(ftrs.geometry().bounds()))\
          .select('EVI')


    # // // Generate multiband image from all features

    def getLayersbyFeature(ftr):
        newf = ee.Feature(ftr)
                                                    
        
        #local vegetation indices
        et = ET.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)
        evi = VEGETATION.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)
                                        
                                                    

        #   // local vegetation indices
        et = ET.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)
        evi = VEGETATION.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)
        
        #Temperature fluctuation
                                                    
        LST_phase = LST_phase_Temp.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol).rename(['cos'],['pha'])
        LST_amp = LST_amplitude_Temp.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol).rename(['cos'],['amp'])
  

        #  derivative vegetation index
        evi_sd = VEGETATION_VARIABILITY.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)
        # evi_sd = VEGETATION_VARIABILITY.median().reduceRegion(ee.Reducer.median(), newf.geometry(), resol)

        #   // local terrain attributes 
        hil = hillshade.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol)
        slp = slope.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol)

        #   // secondary terrain attributes 
        hcrv = Hor_curv.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol)
        vcrv = Ver_curv.reduceRegion(ee.Reducer.mean(), newf.geometry(), resol)

        newf = newf\
        .set(LST_phase)\
        .set(LST_amp)\
        .set(et)\
        .set(evi)\
        .set(evi_sd)\
        .set(hil)\
        .set(slp)\
        .set(hcrv)\
        .set(vcrv)

        # if the value is not null, set the values as a property of the feature. The name of the property will be the date
        return ee.Feature(newf)
   
    
    Extract = ftrs.map(lambda ftr: getLayersbyFeature(ftr))

    
    return Extract