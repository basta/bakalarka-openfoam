# state file generated using paraview version 5.13.3
import paraview
paraview.compatibility.major = 5
paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1020, 663]
renderView1.AxesGrid = 'Grid Axes 3D Actor'
renderView1.CenterOfRotation = [0.0, 0.0, 0.004999999888241291]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [0.027529799707233366, -0.08892604597823897, 0.3202058853476469]
renderView1.CameraFocalPoint = [-1.4555109329408168e-18, 2.6332918645431033e-18, 0.0049999998882412945]
renderView1.CameraViewUp = [-0.3095052880658176, 0.907768034814641, 0.28313189793461996]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 0.08506468252434378
renderView1.LegendGrid = 'Legend Grid Actor'
renderView1.PolarGrid = 'Polar Grid Actor'
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'Grid Axes 3D Actor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1020, 663)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Open FOAM Reader'
expfoam = OpenFOAMReader(registrationName='exp.foam', FileName='/home/basta/Projects/bakalarka-openfoam/bakalarka-openfoam/cyl/exp.foam')
expfoam.MeshRegions = ['internalMesh', 'patch/coldBlockInterface', 'patch/defaultFaces', 'patch/hotBlockInterface', 'patch/wallBottomFluid', 'patch/wallTopFluid']
expfoam.CellArrays = ['F', 'T', 'U', 'acceleration', 'gradT', 'p', 'p_rgh', 'wallHeatFlux']

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=expfoam)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]
slice1.PointMergeMethod = 'Uniform Binning'

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, 0.0, 0.004999999888241291]
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.004999999888241291]

# create a new 'Cell Centers'
cellCenters1 = CellCenters(registrationName='CellCenters1', Input=slice1)

# create a new 'Stream Tracer'
streamTracer1 = StreamTracer(registrationName='StreamTracer1', Input=expfoam,
    SeedType='Line')
streamTracer1.Vectors = ['POINTS', 'U']
streamTracer1.MaximumStreamlineLength = 0.10000000149011612

# init the 'Line' selected for 'SeedType'
streamTracer1.SeedType.Point1 = [-0.05000000074505806, -0.05000000074505806, 0.0]
streamTracer1.SeedType.Point2 = [0.05000000074505806, 0.05000000074505806, 0.009999999776482582]

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=cellCenters1,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'U']
glyph1.ScaleArray = ['POINTS', 'U']
glyph1.ScaleFactor = 52.41520509848052
glyph1.GlyphTransform = 'Transform2'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'U'
uTF2D = GetTransferFunction2D('U')

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.TransferFunction2D = uTF2D
uLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 0.000127325765838697, 0.865003, 0.865003, 0.865003, 0.000254651531677394, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = ['POINTS', 'U']
slice1Display.LookupTable = uLUT
slice1Display.SelectNormalArray = 'None'
slice1Display.SelectTangentArray = 'None'
slice1Display.SelectTCoordArray = 'None'
slice1Display.TextureTransform = 'Transform2'
slice1Display.OSPRayScaleArray = 'p'
slice1Display.OSPRayScaleFunction = 'Piecewise Function'
slice1Display.Assembly = 'Hierarchy'
slice1Display.SelectedBlockSelectors = ['']
slice1Display.SelectOrientationVectors = 'U'
slice1Display.ScaleFactor = 0.010000000149011612
slice1Display.SelectScaleArray = 'p'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'p'
slice1Display.GaussianRadius = 0.0005000000074505806
slice1Display.SetScaleArray = ['POINTS', 'p']
slice1Display.ScaleTransferFunction = 'Piecewise Function'
slice1Display.OpacityArray = ['POINTS', 'p']
slice1Display.OpacityTransferFunction = 'Piecewise Function'
slice1Display.DataAxesGrid = 'Grid Axes Representation'
slice1Display.PolarAxes = 'Polar Axes Representation'
slice1Display.SelectInputVectors = ['POINTS', 'U']
slice1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [-0.0003568099928088486, 0.0, 0.5, 0.0, 0.0002314912562724203, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [-0.0003568099928088486, 0.0, 0.5, 0.0, 0.0002314912562724203, 1.0, 0.5, 0.0]

# show data from streamTracer1
streamTracer1Display = Show(streamTracer1, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'p'
pTF2D = GetTransferFunction2D('p')

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')
pLUT.TransferFunction2D = pTF2D
pLUT.RGBPoints = [-0.0007310601067729294, 0.231373, 0.298039, 0.752941, 6.468029459938407e-05, 0.865003, 0.865003, 0.865003, 0.0008604206959716976, 0.705882, 0.0156863, 0.14902]
pLUT.ScalarRangeInitialized = 1.0

# trace defaults for the display properties.
streamTracer1Display.Representation = 'Surface'
streamTracer1Display.ColorArrayName = ['POINTS', 'p']
streamTracer1Display.LookupTable = pLUT
streamTracer1Display.SelectNormalArray = 'None'
streamTracer1Display.SelectTangentArray = 'None'
streamTracer1Display.SelectTCoordArray = 'None'
streamTracer1Display.TextureTransform = 'Transform2'
streamTracer1Display.OSPRayScaleArray = 'p'
streamTracer1Display.OSPRayScaleFunction = 'Piecewise Function'
streamTracer1Display.Assembly = ''
streamTracer1Display.SelectedBlockSelectors = ['']
streamTracer1Display.SelectOrientationVectors = 'Normals'
streamTracer1Display.ScaleFactor = 0.00998695157468319
streamTracer1Display.SelectScaleArray = 'p'
streamTracer1Display.GlyphType = 'Arrow'
streamTracer1Display.GlyphTableIndexArray = 'p'
streamTracer1Display.GaussianRadius = 0.0004993475787341595
streamTracer1Display.SetScaleArray = ['POINTS', 'p']
streamTracer1Display.ScaleTransferFunction = 'Piecewise Function'
streamTracer1Display.OpacityArray = ['POINTS', 'p']
streamTracer1Display.OpacityTransferFunction = 'Piecewise Function'
streamTracer1Display.DataAxesGrid = 'Grid Axes Representation'
streamTracer1Display.PolarAxes = 'Polar Axes Representation'
streamTracer1Display.SelectInputVectors = ['POINTS', 'Normals']
streamTracer1Display.WriteLog = ''

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
streamTracer1Display.ScaleTransferFunction.Points = [-0.0007310601067729294, 0.0, 0.5, 0.0, 0.0008604206959716976, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
streamTracer1Display.OpacityTransferFunction.Points = [-0.0007310601067729294, 0.0, 0.5, 0.0, 0.0008604206959716976, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Position = [0.8862745098039215, 0.01809954751131222]
uLUTColorBar.Title = 'U'
uLUTColorBar.ComponentTitle = 'Magnitude'

# set color bar visibility
uLUTColorBar.Visibility = 1

# get color legend/bar for pLUT in view renderView1
pLUTColorBar = GetScalarBar(pLUT, renderView1)
pLUTColorBar.WindowLocation = 'Upper Right Corner'
pLUTColorBar.Title = 'p'
pLUTColorBar.ComponentTitle = ''

# set color bar visibility
pLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
streamTracer1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')
pPWF.Points = [-0.0007310601067729294, 0.0, 0.5, 0.0, 0.0008604206959716976, 1.0, 0.5, 0.0]
pPWF.ScalarRangeInitialized = 1

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Points = [0.0, 0.0, 0.5, 0.0, 0.000254651531677394, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.ViewModules = renderView1
animationScene1.Cues = timeAnimationCue1
animationScene1.AnimationTime = 1109.8365157655971
animationScene1.StartTime = 9.836516
animationScene1.EndTime = 1109.8365157655971
animationScene1.PlayMode = 'Snap To TimeSteps'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------