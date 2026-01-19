# Direct Recorded Script
# PowerVIZ 6-2022-R4 ( 6.2.3 )
# Date: Thu Oct 24 10:49:17 2024

#import numpy as np


lower_r = 0 # Lower limit to get the images from
upper_r = 180 # Upper limit to get the images
spacing_r = 10 # Must be the same as in the coordinates creation script
variable_n = 'Vorticity Magnitude' # Name of Variable file as first argument after code name
variable = 'omega' # Name of variable to save file as first argument after code name
img_number = 0
timestep_i = 1855664
timestep_2i = 1856688
timestep_n = 2058416
angular_n = 180 
desire_rev = 1
timestep_rec = timestep_2i - timestep_i # To check in PowerViz
print('JR')

#radius_all = np.linspace(lower_r,upper_r+1,spacing_r)
#radius = np.array([50, 55, 60, 65, 85, 90, 95, 100])
radius_all = [  0,  10,  20,  30,  40,  50,  60,  70,  80,  90, 100, 110, 120,
       130, 140, 150, 160, 170, 180]
timesteps = [1855664, 1860784, 1865904, 1871024, 1876144, 1881264, 1886384,
       1891504, 1896624, 1901744, 1906864, 1911984, 1917104, 1922224,
       1927344, 1932464, 1937584, 1942704, 1947824, 1952944, 1958064,
       1963184, 1968304, 1973424, 1978544, 1983664, 1988784, 1993904,
       1999024, 2004144, 2009264, 2014384, 2019504, 2024624, 2029744,
       2034864, 2039984, 2045104]
#timesteps = np.arange(timestep_i,timestep_n+1,timestep_rec) # Till the end of the PowerViz file
#timesteps = np.arange(timestep_i,timestep_i+(timestep_rec*angular_n*desire_rev),timestep_rec) # Defined number of revolutions
#timesteps = [1638592]
cdi_path = 'Baseline_7e-5_nowalls_1e-4.cdi' # Name of CDI file as first argument after code name
fnc_path = '1rev-SMR-VR8.fnc' # Name of CDI file as second argument after code name
snc_path = ""
variables = ["Velocity Magnitude","X-Velocity","Y-Velocity","Z-Velocity","Vorticity Magnitude","Static Pressure","Density","Temperature","Total Pressure","Turb Kinetic Energy","Lambda2"]
image_path = ""

#Setting files
project1=app.newProject(CDIFilename=cdi_path, fluidFilename=fnc_path, name="Project1", variablesToLoad=variables)
project1.setUpdateEnabledAllViewers(0)

# Hidding walls and inlet geometries
baseAssembly1=project1.baseAssembly
modelView5=baseAssembly1.defaultModelView
part6=project1.getEntityByPath("/inlet1")
modelViewObject777=modelView5.getModelViewObject(part6)
modelViewObject777.displayMode="Hidden"
part7=project1.getEntityByPath("/inlet2")
modelViewObject778=modelView5.getModelViewObject(part7)
modelViewObject778.displayMode="Hidden"
part8=project1.getEntityByPath("/inlet3")
modelViewObject779=modelView5.getModelViewObject(part8)
modelViewObject779.displayMode="Hidden"
project1.setUpdateEnabledAllViewers(1)
project1.setUpdateEnabledAllViewers(0)
part12=project1.getEntityByPath("/outlet1")
modelViewObject783=modelView5.getModelViewObject(part12)
modelViewObject783.displayMode="Hidden" 
part13=project1.getEntityByPath("/outlet2")
modelViewObject784=modelView5.getModelViewObject(part13)
modelViewObject784.displayMode="Hidden"
part14=project1.getEntityByPath("/outlet3")
modelViewObject785=modelView5.getModelViewObject(part14)
modelViewObject785.displayMode="Hidden" 
project1.setUpdateEnabledAllViewers(1) 
project1.setUpdateEnabledAllViewers(0)

# Setting background color
project1.setUpdateEnabledAllViewers(1)
viewerBackground3=project1.get(name="ViewerBackground3", type="ViewerBackground")
viewerBackground3.color="#ffffff"

slice1=project1.get(name="Slice1", type="Slice")
slice1.setSelected()
scalarPropertySet9=project1.get(name=variable_n, type="ScalarPropertySet")
slice1.imageSPS=scalarPropertySet9
slice1.fitSPSRange(referenceSPS=scalarPropertySet9)

dt = 1
for t in timesteps:

       project1=app.currentProject
       project1.timeStep=t
       timeAnimation1=project1.timeAnimation
       timeAnimation1.timestep=t

       for radius in radius_all:

       # Setting variable to display
              referenceGeometry4=project1.newReferenceGeometry(name="ReferenceGeometry1")
              part56=project1.getEntityByPath("/plane_fnc_radial_{}deg".format(radius))
              referenceGeometry4.createGeometry(args={"entities":[part56]}, coordinateSystem='default_csys', definedVia="Entities", invertNormals=False, type="Copy of CDI Geometry")
              referenceGeometry4.visibility="Hide Reference Geometry"
              slice1.referenceGeometry=referenceGeometry4
              slice1.fitSPSRange(referenceSPS=scalarPropertySet9)
              slice1.referenceFrame=None
              # scalarPropertySet9.range=( (0.2374841, 35), "m/sec")
              # slice4.contourLinesRange=( (0.2374841, 35), "m/sec")
              slice1.saveAsciiData(filename="avg-radial-cut-{}deg-{}.txt".format(radius,variable), dataSamplesPerCellWidth=0.25, significantDigits=7)
       dt += 1