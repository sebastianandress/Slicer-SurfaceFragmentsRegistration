import os
import unittest
import logging
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import numpy as np

#
# SurfaceFragmentsRegistration
#

class SurfaceFragmentsRegistration(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Surface Fragments Registration"
    self.parent.categories = ["Registration"]
    self.parent.dependencies = []
    self.parent.contributors = ["Sebastian Andress (LMU Munich)"]
    self.parent.helpText = ""
    # TODO: replace with organization, grant and thanks
    self.parent.acknowledgementText = """
This file was originally developed by Sebastian Andress, LMU Munich.
"""

    # Additional initialization step after application startup is complete
    slicer.app.connect("startupCompleted()", registerSampleData)

#
# Register sample data sets in Sample Data module
#

def registerSampleData():
  """
  Add data sets to Sample Data module.
  """
  # It is always recommended to provide sample data for users to make it easy to try the module,
  # but if no sample data is available then this method (and associated startupCompeted signal connection) can be removed.

  import SampleData
  # To ensure that the source code repository remains small (can be downloaded and installed quickly)
  # it is recommended to store data sets that are larger than a few MB in a Github release.

  SampleData.SampleDataLogic.registerCustomSampleDataSource(
      sampleName='Surface Fragments Registration Models',
      category='SurfaceFragmentsRegistration',
      uris=['https://raw.githubusercontent.com/sebastianandress/Slicer-SurfaceFragmentsRegistration/master/SurfaceFragmentsRegistration/Testing/Assets/sourceModel.vtk', 'https://raw.githubusercontent.com/sebastianandress/Slicer-SurfaceFragmentsRegistration/master/SurfaceFragmentsRegistration/Testing/Assets/targetModel.vtk'],
      loadFiles=[True, True],
      fileNames=['sourceModel.vtk','targetModel.vtk'],
      nodeNames=['sourceModel','targetModel'],
      loadFileType=['ModelFile', 'ModelFile'],
)

#
# SurfaceFragmentsRegistrationWidget
#

class SurfaceFragmentsRegistrationWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    VTKObservationMixin.__init__(self)  # needed for parameter node observation
    self.logic = None
    self._parameterNode = None
    self._updatingGUIFromParameterNode = False

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer).
    # Additional widgets can be instantiated manually and added to self.layout.
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/SurfaceFragmentsRegistration.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    # Set scene in MRML widgets. Make sure that in Qt designer the top-level qMRMLWidget's
    # "mrmlSceneChanged(vtkMRMLScene*)" signal in is connected to each MRML widget's.
    # "setMRMLScene(vtkMRMLScene*)" slot.
    uiWidget.setMRMLScene(slicer.mrmlScene)

    # Create logic class. Logic implements all computations that should be possible to run
    # in batch mode, without a graphical user interface.
    self.logic = SurfaceFragmentsRegistrationLogic()
    self.logic.logCallback = self.addLog

    # Connections

    # These connections ensure that we update parameter node when scene is closed
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
    self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

    # These connections ensure that whenever user changes some settings on the GUI, that is saved in the MRML scene
    # (in the selected parameter node).
    self.ui.sourceModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.targetModelSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.initializationClusterRadiusSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.minimalClusterAreaSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.cutoffThresholdSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.preRegistrationCB.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.sourceLandmarkSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.targetLandmarkSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.updateParameterNodeFromGUI)
    self.ui.initializationIterationsSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.openingWidthSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.maximalIterationsSelector.connect("valueChanged(double)", self.updateParameterNodeFromGUI)
    self.ui.createTransformationsCB.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.markDeviationsCB.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.markFragmentsCB.connect("toggled(bool)", self.updateParameterNodeFromGUI)
    self.ui.fragmentSelectorSB.connect("valueChanged(int)", self.updateParameterNodeFromGUI)

    # Buttons
    self.ui.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.ui.defaultsButton.connect('clicked(bool)', self.onDefaultsButton)
    self.ui.fragmentSelectorSB.connect('valueChanged(int)', self.onFragmentSelector)

    # Make sure parameter node is initialized (needed for module reload)
    self.initializeParameterNode()

  def cleanup(self):
    """
    Called when the application closes and the module widget is destroyed.
    """
    self.removeObservers()

  def enter(self):
    """
    Called each time the user opens this module.
    """
    # Make sure parameter node exists and observed
    self.initializeParameterNode()

  def exit(self):
    """
    Called each time the user opens a different module.
    """
    # Do not react to parameter node changes (GUI wlil be updated when the user enters into the module)
    self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

  def onSceneStartClose(self, caller, event):
    """
    Called just before the scene is closed.
    """
    # Parameter node will be reset, do not use it anymore
    self.setParameterNode(None)

  def onSceneEndClose(self, caller, event):
    """
    Called just after the scene is closed.
    """
    # If this module is shown while the scene is closed then recreate a new parameter node immediately
    if self.parent.isEntered:
      self.initializeParameterNode()

  def initializeParameterNode(self):
    """
    Ensure parameter node exists and observed.
    """
    # Parameter node stores all user choices in parameter values, node selections, etc.
    # so that when the scene is saved and reloaded, these settings are restored.

    self.setParameterNode(self.logic.getParameterNode())

    # Select default input nodes if nothing is selected yet to save a few clicks for the user
    firstModelNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLModelNode")
    if firstModelNode:
      if not self._parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL):
        self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, firstModelNode.GetID())
      if not self._parameterNode.GetNodeReference(PARAMETER_TARGETMODEL):
        self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, firstModelNode.GetID())

    firstFiducialNode = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLMarkupsFiducialNode")
    if firstFiducialNode:
      if not self._parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL):
        self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, firstFiducialNode.GetID())
      if not self._parameterNode.GetNodeReference(PARAMETER_TARGETMODEL):
        self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, firstFiducialNode.GetID())

  def setParameterNode(self, inputParameterNode):
    """
    Set and observe parameter node.
    Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
    """

    if inputParameterNode:
      self.logic.setDefaultParameters(inputParameterNode)

    # Unobserve previously selected parameter node and add an observer to the newly selected.
    # Changes of parameter node are observed so that whenever parameters are changed by a script or any other module
    # those are reflected immediately in the GUI.
    if self._parameterNode is not None:
      self.removeObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)
    self._parameterNode = inputParameterNode
    if self._parameterNode is not None:
      self.addObserver(self._parameterNode, vtk.vtkCommand.ModifiedEvent, self.updateGUIFromParameterNode)

    # Initial GUI update
    self.updateGUIFromParameterNode()

  def updateGUIFromParameterNode(self, caller=None, event=None):
    """
    This method is called whenever parameter node is changed.
    The module GUI is updated to show the current state of the parameter node.
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    # Make sure GUI changes do not call updateParameterNodeFromGUI (it could cause infinite loop)
    self._updatingGUIFromParameterNode = True

    self.ui.sourceModelSelector.setCurrentNode(self._parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL))
    self.ui.targetModelSelector.setCurrentNode(self._parameterNode.GetNodeReference(PARAMETER_TARGETMODEL))
    self.ui.initializationClusterRadiusSelector.value = float(self._parameterNode.GetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS))
    self.ui.minimalClusterAreaSelector.value = float(self._parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA))
    self.ui.cutoffThresholdSelector.value = float(self._parameterNode.GetParameter(PARAMETER_CUTOFFTHRESHOLD))
    self.ui.preRegistrationCB.checked = (self._parameterNode.GetParameter(PARAMETER_PREREGISTRATION) == "true")
    self.ui.sourceLandmarkSelector.setCurrentNode(self._parameterNode.GetNodeReference(PARAMETER_SOURCELANDMARKS))
    self.ui.targetLandmarkSelector.setCurrentNode(self._parameterNode.GetNodeReference(PARAMETER_TARGETLANDMARKS))
    self.ui.initializationIterationsSelector.value = float(self._parameterNode.GetParameter(PARAMETER_INITIALIZATIONITERATIONS))
    self.ui.openingWidthSelector.value = float(self._parameterNode.GetParameter(PARAMETER_OPENINGWIDTH))
    self.ui.maximalIterationsSelector.value = float(self._parameterNode.GetParameter(PARAMETER_MAXIMALITERATIONS))
    self.ui.createTransformationsCB.checked = (self._parameterNode.GetParameter(PARAMETER_ADDTRANSFORMATION) == "true")
    self.ui.markDeviationsCB.checked = (self._parameterNode.GetParameter(PARAMETER_MARKDEVIATION) == "true")
    self.ui.markFragmentsCB.checked = (self._parameterNode.GetParameter(PARAMETER_MARKFRAGMENTS) == "true")

    # Update buttons states and tooltips
    if self._parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL) and self._parameterNode.GetNodeReference(PARAMETER_TARGETMODEL):
      self.ui.applyButton.toolTip = "Run the algorithm."
      self.ui.applyButton.enabled = True
    else:
      self.ui.applyButton.toolTip = "Select source and target models first."
      self.ui.applyButton.enabled = False
    if self._parameterNode.GetParameter(PARAMETER_PREREGISTRATION) == "true":
      self.ui.landmarksCollapibleBox.enabled = True
      self.ui.landmarksCollapibleBox.collapsed = False
      self.ui.sourceLandmarkSelector.enabled = True
      self.ui.targetLandmarkSelector.enabled = True
    else:
      self.ui.landmarksCollapibleBox.enabled = False
      self.ui.landmarksCollapibleBox.collapsed = True
      self.ui.sourceLandmarkSelector.enabled = False
      self.ui.targetLandmarkSelector.enabled = False

    # All the GUI updates are done
    self._updatingGUIFromParameterNode = False

  def updateParameterNodeFromGUI(self, caller=None, event=None):
    """
    This method is called when the user makes any change in the GUI.
    The changes are saved into the parameter node (so that they are restored when the scene is saved and loaded).
    """

    if self._parameterNode is None or self._updatingGUIFromParameterNode:
      return

    wasModified = self._parameterNode.StartModify()  # Modify all properties in a single batch

    self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, self.ui.sourceModelSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID(PARAMETER_TARGETMODEL, self.ui.targetModelSelector.currentNodeID)
    self._parameterNode.SetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS, str(self.ui.initializationClusterRadiusSelector.value))
    self._parameterNode.SetParameter(PARAMETER_MINIMALCLUSTERAREA, str(self.ui.minimalClusterAreaSelector.value))
    self._parameterNode.SetParameter(PARAMETER_CUTOFFTHRESHOLD, str(self.ui.cutoffThresholdSelector.value))
    self._parameterNode.SetParameter(PARAMETER_PREREGISTRATION, "true" if self.ui.preRegistrationCB.checked else "false")
    self._parameterNode.SetNodeReferenceID(PARAMETER_SOURCELANDMARKS, self.ui.sourceLandmarkSelector.currentNodeID)
    self._parameterNode.SetNodeReferenceID(PARAMETER_TARGETLANDMARKS, self.ui.targetLandmarkSelector.currentNodeID)
    self._parameterNode.SetParameter(PARAMETER_INITIALIZATIONITERATIONS, str(self.ui.initializationIterationsSelector.value))
    self._parameterNode.SetParameter(PARAMETER_OPENINGWIDTH, str(self.ui.openingWidthSelector.value))
    self._parameterNode.SetParameter(PARAMETER_MAXIMALITERATIONS, str(self.ui.maximalIterationsSelector.value))
    self._parameterNode.SetParameter(PARAMETER_ADDTRANSFORMATION, "true" if self.ui.createTransformationsCB.checked else "false")
    self._parameterNode.SetParameter(PARAMETER_MARKDEVIATION, "true" if self.ui.markDeviationsCB.checked else "false")
    self._parameterNode.SetParameter(PARAMETER_MARKFRAGMENTS, "true" if self.ui.markFragmentsCB.checked else "false")

    self._parameterNode.EndModify(wasModified)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    if self.ui.applyButton.text == 'Cancel':
      self.logic.requestCancel()
      self.addLog('Cancel requested...')
      return

    self.ui.applyButton.text = 'Cancel'
    errorMessage = None
    qt.QApplication.setOverrideCursor(qt.Qt.WaitCursor)
    
    try:
      self.logic.process(self._parameterNode)
    except Exception as e:
      import traceback
      traceback.print_exc()
      errorMessage = str(e)

    slicer.util.showStatusMessage("")
    self.ui.applyButton.text = 'Apply'
    qt.QApplication.restoreOverrideCursor()
    self.ui.fragmentSelectorSB.value = 1

    if errorMessage:
      slicer.util.errorDisplay("Registration failed: " + errorMessage)
    

  def onDefaultsButton(self):
    self.logic.setDefaultParameters(self._parameterNode, overwrite=True)

  def onFragmentSelector(self, value):
    source = self._parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL)
    if not source:
      return

    trfNodes = slicer.util.getNodesByClass("vtkMRMLLinearTransformNode")
    trfNrs = list()
    for trf in trfNodes:
      if trf.GetName().startswith(PREFIX_TRANSFORMATION):
        trfNrs.append(int(float(trf.GetName().split(PREFIX_TRANSFORMATION)[-1])))
    
    devNrs = list()
    for a in range(source.GetPolyData().GetPointData().GetNumberOfArrays()):
      if source.GetPolyData().GetPointData().GetArrayName(a).startswith(PREFIX_DEVIATION):
        devNrs.append(int(float(source.GetPolyData().GetPointData().GetArrayName(a).split(PREFIX_DEVIATION)[-1])))

    if len(devNrs) == 0:
      return

    newVal = min(devNrs, key=lambda x:abs(x-value))

    source.GetDisplayNode().SetScalarVisibility(True)
    source.GetDisplayNode().SetActiveScalarName(PREFIX_DEVIATION + str(newVal))

    if newVal in trfNrs:
      source.SetAndObserveTransformNodeID(slicer.util.getNode(PREFIX_TRANSFORMATION + str(newVal)).GetID())
    
    self.ui.fragmentSelectorSB.value = newVal


  def addLog(self, text):
    slicer.util.showStatusMessage(text)
    slicer.app.processEvents() # force update


#
# SurfaceFragmentsRegistrationLogic
#

class SurfaceFragmentsRegistrationLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self):
    """
    Called when the logic class is instantiated. Can be used for initializing member variables.
    """
    ScriptedLoadableModuleLogic.__init__(self)
    self.logCallback = None
    self.cancelRequested = False

  def setDefaultParameters(self, parameterNode, overwrite=False):
    """
    Initialize parameter node with default settings.
    """
    if overwrite or not parameterNode.GetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS):
      parameterNode.SetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS, str(DEFAULT_INITIALIZATIONCLUSTERRADIUS))
    if overwrite or not parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA):
      parameterNode.SetParameter(PARAMETER_MINIMALCLUSTERAREA, str(DEFAULT_MINIMALCLUSTERAREA))
    if overwrite or not parameterNode.GetParameter(PARAMETER_CUTOFFTHRESHOLD):
      parameterNode.SetParameter(PARAMETER_CUTOFFTHRESHOLD, str(DEFAULT_CUTOFFTHRESHOLD))
    if overwrite or not parameterNode.GetParameter(PARAMETER_PREREGISTRATION):
      parameterNode.SetParameter(PARAMETER_PREREGISTRATION, "true" if DEFAULT_PREREGISTRATION else "false")
    if overwrite or not parameterNode.GetParameter(PARAMETER_INITIALIZATIONITERATIONS):
      parameterNode.SetParameter(PARAMETER_INITIALIZATIONITERATIONS, str(DEFAULT_INITIALIZATIONITERATIONS))
    if overwrite or not parameterNode.GetParameter(PARAMETER_OPENINGWIDTH):
      parameterNode.SetParameter(PARAMETER_OPENINGWIDTH, str(DEFAULT_OPENINGWIDTH))
    if overwrite or not parameterNode.GetParameter(PARAMETER_MAXIMALITERATIONS):
      parameterNode.SetParameter(PARAMETER_MAXIMALITERATIONS, str(DEFAULT_MAXIMALITERATIONS))
    if overwrite or not parameterNode.GetParameter(PARAMETER_ADDTRANSFORMATION):
      parameterNode.SetParameter(PARAMETER_ADDTRANSFORMATION, "true" if DEFAULT_ADDTRANSFORMATION else "false")
    if overwrite or not parameterNode.GetParameter(PARAMETER_MARKDEVIATION):
      parameterNode.SetParameter(PARAMETER_MARKDEVIATION, "true" if DEFAULT_MARKDEVIATION else "false")
    if overwrite or not parameterNode.GetParameter(PARAMETER_MARKFRAGMENTS):
      parameterNode.SetParameter(PARAMETER_MARKFRAGMENTS, "true" if DEFAULT_MARKFRAGMENTS else "false")


  def process(self, parameterNode):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param parameterNode: node containing all parameters
    """
    
    # check if input nodes are valid
    if not parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL) or not parameterNode.GetNodeReference(PARAMETER_TARGETMODEL):
      raise ValueError("Source or target model is invalid")
    if parameterNode.GetParameter(PARAMETER_PREREGISTRATION) == "true" and (not parameterNode.GetNodeReference(PARAMETER_SOURCELANDMARKS) or not parameterNode.GetNodeReference(PARAMETER_TARGETLANDMARKS)):
      raise ValueError("Source or target landmarks is invalid")

    self.setDefaultParameters(parameterNode)

    self.cancelRequested = False

    # remove conflicting arrays by name
    ar_names = list()
    for a in range(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData().GetPointData().GetNumberOfArrays()):
      ar_names.append(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData().GetPointData().GetArrayName(a))
    for a in reversed(ar_names):
      if a == None: continue
      elif a.startswith(PREFIX_DEVIATION) or a.startswith(PREFIX_FRAGMENT) or a in [LET_INTERMEDIATEDISTANCES, LET_INTERMEDIATEIDS, LET_INTERMEDIATEREGISTERED]:
        parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData().GetPointData().RemoveArray(a)

    # remove conflicting transformations by name
    trfNodes = slicer.util.getNodesByClass("vtkMRMLLinearTransformNode")
    for node in reversed(trfNodes):
        if node.GetName().startswith(PREFIX_TRANSFORMATION):
          slicer.mrmlScene.RemoveNode(node)

    # initial landmark registration
    if parameterNode.GetParameter(PARAMETER_PREREGISTRATION) == "true":
      lmTrf = self._regLandmark(parameterNode.GetNodeReference(PARAMETER_SOURCELANDMARKS), parameterNode.GetNodeReference(PARAMETER_TARGETLANDMARKS))
    else:
      lmTrf = vtk.vtkTransform()

    sourcePD = vtk.vtkPolyData()
    sourcePD.DeepCopy(self._transformPD(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData(), lmTrf))
    # sourcePD.DeepCopy(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData())

    idArrayId = self._addArray(sourcePD, range(sourcePD.GetNumberOfPoints()), LET_INTERMEDIATEIDS)
    regArrayId = self._addArray(sourcePD, [0]*(sourcePD.GetNumberOfPoints()), LET_INTERMEDIATEREGISTERED)


    clusters = 1
    for it in range(int(float(parameterNode.GetParameter(PARAMETER_MAXIMALITERATIONS)))):
      self._log("Iteration %s/%s: %s clusters" % (it, int(float(parameterNode.GetParameter(PARAMETER_MAXIMALITERATIONS))), clusters-1))
      if self.cancelRequested:
        break
      smallestTrf = vtk.vtkTransform()
      smallestDev = -1
      smallestIds = list()

      # extract remain model
      remainModel = self._thresholdPD(sourcePD, LET_INTERMEDIATEREGISTERED, 0,0)

      # check for large enough clusters
      remPDList = self._connectivityPD(remainModel, largest=False)
      remPDSizes = list()
      for remPD in remPDList:
        remPDSizes.append(self._getSurfaceArea(remPD))
      
      for i, size in enumerate(remPDSizes):
        if size < float(parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA)):
          continue

        for _ in range(int(float(parameterNode.GetParameter(PARAMETER_INITIALIZATIONITERATIONS)))):
          # extract random piece
          vertexId = np.random.randint(0,remPDList[i].GetNumberOfPoints())
          extractPD = self._extractDistancePD(remPDList[i], vertexId, float(parameterNode.GetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS)))

          # register piece
          trf = self._regICP(extractPD, parameterNode.GetNodeReference(PARAMETER_TARGETMODEL).GetPolyData())
          extractTrfPD = self._transformPD(extractPD, trf)

          # get/save mean deviation/piece/trf
          dev = self._calculateDistances(extractTrfPD, parameterNode.GetNodeReference(PARAMETER_TARGETMODEL).GetPolyData())
          meanDev = np.mean(dev)

          if smallestDev == -1 or smallestDev > meanDev:
            smallestDev = meanDev
            smallestTrf = trf
            smallestIds = vtk.util.numpy_support.vtk_to_numpy(extractTrfPD.GetPointData().GetArray(self.getArrayNumber(extractPD, LET_INTERMEDIATEIDS)))

      # breaks loop if no piece remaining or all pieces smaller then float(parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA))
      print('smallest dev: ', smallestDev)
      if smallestDev == -1:
        break
      
      # register whole model, calculate deviations
      trfPD = self._transformPD(sourcePD, smallestTrf)
      dev = self._calculateDistances(trfPD, parameterNode.GetNodeReference(PARAMETER_TARGETMODEL).GetPolyData())
      self._addArray(trfPD, dev, LET_INTERMEDIATEDISTANCES)

      # extract all pieces with distance <x
      extractPD = self._thresholdPD(trfPD, LET_INTERMEDIATEDISTANCES, 0, float(parameterNode.GetParameter(PARAMETER_CUTOFFTHRESHOLD)))
      extractPD.GetPointData().RemoveArray(LET_INTERMEDIATEDISTANCES)

      for sub in range(int(float(parameterNode.GetParameter(PARAMETER_OPENINGWIDTH)))+1):
        sgPD = self._openingPD(extractPD, int(float(parameterNode.GetParameter(PARAMETER_OPENINGWIDTH)))-sub)
        extractPdList = self._connectivityPD(sgPD, largest=False)
        try:
          mcId = self._regionMostCommonIds(extractPdList, LET_INTERMEDIATEIDS, smallestIds, float(parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA)))
        except IndexError:
          continue

        sgPD = extractPdList[mcId]
        print('shrinkgrow: ', int(float(parameterNode.GetParameter(PARAMETER_OPENINGWIDTH)))-sub, self._getSurfaceArea(sgPD))
        
        if self._getSurfaceArea(sgPD) >= float(parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA)):
          break
          
      if self._getSurfaceArea(sgPD) < float(parameterNode.GetParameter(PARAMETER_MINIMALCLUSTERAREA)):
        continue
      

      trf = self._regICP(sgPD, parameterNode.GetNodeReference(PARAMETER_TARGETMODEL).GetPolyData())
      trfPD = self._transformPD(trfPD, trf)
      
      dev = self._calculateDistances(trfPD, parameterNode.GetNodeReference(PARAMETER_TARGETMODEL).GetPolyData())
      self._addArray(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData(), dev, PREFIX_DEVIATION + str(clusters))
      
      # mark island as registered
      eia = vtk.util.numpy_support.vtk_to_numpy(sgPD.GetPointData().GetArray(self.getArrayNumber(sgPD, LET_INTERMEDIATEIDS)))
      sra = vtk.util.numpy_support.vtk_to_numpy(sourcePD.GetPointData().GetArray(LET_INTERMEDIATEREGISTERED))
      
      np.put(sra, eia, 1)
      regArrayId = self._addArray(sourcePD, sra, LET_INTERMEDIATEREGISTERED)
      
      if parameterNode.GetParameter(PARAMETER_ADDTRANSFORMATION) == "true":
        t = self._multiplyTransforms([lmTrf, smallestTrf, trf])
        self._addTrf(t, name=PREFIX_TRANSFORMATION + str(clusters))

      if parameterNode.GetParameter(PARAMETER_MARKFRAGMENTS) == "true":
        ca = np.zeros(sourcePD.GetNumberOfPoints())
        np.put(ca, eia, 1)
        self._addArray(parameterNode.GetNodeReference(PARAMETER_SOURCEMODEL).GetPolyData(), ca, PREFIX_FRAGMENT + str(clusters))
      clusters += 1


  
  def requestCancel(self):
    logging.info("User requested cancelling.")
    self.cancelRequested = True

  def _log(self, message):
    if self.logCallback:
      self.logCallback(message)

  #
  # Method's functions
  #
  @staticmethod
  def _openingPD (polydata, width):
    baseArId = polydata.GetPointData().AddArray(vtk.util.numpy_support.numpy_to_vtk(range(polydata.GetNumberOfPoints())))
    basePd = vtk.vtkPolyData()
    basePd.DeepCopy(polydata)

    pd = vtk.vtkPolyData()
    pd.DeepCopy(polydata)

    # erode
    erodewidth = 0
    for w in range(width):

      arId = pd.GetPointData().AddArray(vtk.util.numpy_support.numpy_to_vtk(range(pd.GetNumberOfPoints())))

      featureEdges = vtk.vtkFeatureEdges()
      featureEdges.SetInputData(pd)
      featureEdges.BoundaryEdgesOn()
      featureEdges.FeatureEdgesOff()
      featureEdges.ManifoldEdgesOff()
      featureEdges.NonManifoldEdgesOff()
      featureEdges.Update()

      ar = vtk.util.numpy_support.vtk_to_numpy(featureEdges.GetOutput().GetPointData().GetArray(arId))
      pd.GetPointData().RemoveArray(arId)
      ids = vtk.vtkIdTypeArray()
      ids.SetNumberOfComponents(1)
      for i in ar:
        ids.InsertNextValue(i)


      selectionNode = vtk.vtkSelectionNode()
      selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
      selectionNode.SetContentType(4)
      selectionNode.SetSelectionList(ids)
      selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1)
      selectionNode.GetProperties().Set(vtk.vtkSelectionNode.INVERSE(),1)

      selection = vtk.vtkSelection()
      selection.AddNode(selectionNode)

      extractSelection = vtk.vtkExtractSelection()
      extractSelection.SetInputData(0, pd)
      extractSelection.SetInputData(1, selection)
      extractSelection.Update()

      geometryFilter = vtk.vtkGeometryFilter() 
      geometryFilter.SetInputData(extractSelection.GetOutput()) 
      geometryFilter.Update()

      cleanFilter = vtk.vtkCleanPolyData()
      cleanFilter.SetInputData(geometryFilter.GetOutput())
      cleanFilter.Update()

      if cleanFilter.GetOutput().GetNumberOfPoints() == 0:
        break
      
      erodewidth = w + 1
      pd.DeepCopy(cleanFilter.GetOutput())

    # dilate
    baseAr = set(vtk.util.numpy_support.vtk_to_numpy(pd.GetPointData().GetArray(baseArId)).tolist())
    for _ in range(erodewidth):

      featureEdges = vtk.vtkFeatureEdges()
      featureEdges.SetInputData(pd)
      featureEdges.BoundaryEdgesOn()
      featureEdges.FeatureEdgesOff()
      featureEdges.ManifoldEdgesOff()
      featureEdges.NonManifoldEdgesOff()
      featureEdges.Update()

      ar = vtk.util.numpy_support.vtk_to_numpy(featureEdges.GetOutput().GetPointData().GetArray(baseArId))
      for i in ar:
        conCellsIdList = vtk.vtkIdList()
        basePd.GetPointCells(i, conCellsIdList)

        for c in range(conCellsIdList.GetNumberOfIds()):
          cellsPointIds = basePd.GetCell(conCellsIdList.GetId(c)).GetPointIds()
        
          for p in range(cellsPointIds.GetNumberOfIds()):
            baseAr.add(cellsPointIds.GetId(p))
      
      ids = vtk.vtkIdTypeArray()
      ids.SetNumberOfComponents(1)
      for i in baseAr:
          ids.InsertNextValue(i)
      selectionNode = vtk.vtkSelectionNode()
      selectionNode.SetFieldType(vtk.vtkSelectionNode.POINT)
      selectionNode.SetContentType(4)
      selectionNode.SetSelectionList(ids)
      selectionNode.GetProperties().Set(vtk.vtkSelectionNode.CONTAINING_CELLS(), 1)

      selection = vtk.vtkSelection()
      selection.AddNode(selectionNode)

      extractSelection = vtk.vtkExtractSelection()
      extractSelection.SetInputData(0, basePd)
      extractSelection.SetInputData(1, selection)
      extractSelection.Update()

      geometryFilter = vtk.vtkGeometryFilter() 
      geometryFilter.SetInputData(extractSelection.GetOutput()) 
      geometryFilter.Update()

      cleanFilter = vtk.vtkCleanPolyData()
      cleanFilter.SetInputData(geometryFilter.GetOutput())
      cleanFilter.Update()

      pd.DeepCopy(cleanFilter.GetOutput())

    pd.GetPointData().RemoveArray(baseArId)
    return pd
  

  def _extractDistancePD(self, polydata, seedId, distance):
    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1)
    sphere.SetCenter(polydata.GetPoint(seedId))
    sphere.Update()
    
    if polydata.GetPointData().GetArray(LET_INTERMEDIATEIDS):
      seedId = polydata.GetPointData().GetArray(LET_INTERMEDIATEIDS).GetValue(seedId)
    else:
      self._addArray(polydata, range(polydata.GetNumberOfPoints()), LET_INTERMEDIATEIDS)

    dist = vtk.vtkDistancePolyDataFilter()
    dist.AddInputData(0, polydata)
    dist.AddInputData(1, sphere.GetOutput())
    dist.SignedDistanceOff()
    dist.NegateDistanceOff()
    dist.ComputeSecondDistanceOff()
    dist.Update()
    extractPd = self._thresholdPD(dist.GetOutput(), "Distance", upper=distance)
    
    extrSeedId = vtk.util.numpy_support.vtk_to_numpy(extractPd.GetPointData().GetArray(LET_INTERMEDIATEIDS)).tolist().index(seedId)
    
    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetInputData(extractPd)
    connectivityFilter.SetExtractionModeToPointSeededRegions()
    connectivityFilter.AddSeed(extrSeedId)
    connectivityFilter.Update()
    
    return connectivityFilter.GetOutput()

  @staticmethod
  def _thresholdPD(polydata, array, lower=0, upper=0):
    thresh = vtk.vtkThreshold()
    thresh.SetInputArrayToProcess(0, 0, vtk.vtkDataObject().FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes().SCALARS, array)
    thresh.SetInputData(polydata)
    thresh.ThresholdBetween(lower, upper)

    geom = vtk.vtkGeometryFilter()
    geom.SetInputConnection(thresh.GetOutputPort())

    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(geom.GetOutputPort())
    cleanFilter.Update()

    return cleanFilter.GetOutput()

  @staticmethod
  def _connectivityPD(polydata, largest=True):
    connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
    connectivityFilter.SetInputData(polydata)
    
    if largest:
      connectivityFilter.SetExtractionModeToLargestRegion()
      connectivityFilter.Update()

      cleanFilter = vtk.vtkCleanPolyData()
      cleanFilter.SetInputConnection(connectivityFilter.GetOutputPort())
      cleanFilter.Update()

      return cleanFilter.GetOutput()
    
    else:
      connectivityFilter.SetExtractionModeToAllRegions()
      connectivityFilter.ColorRegionsOn()
      connectivityFilter.Update()

      selectorCon = vtk.vtkThreshold()
      selectorCon.SetInputArrayToProcess(0, 0, vtk.vtkDataObject().FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes().SCALARS, "RegionId")
      selectorCon.SetInputConnection(connectivityFilter.GetOutputPort())
      selectorCon.AllScalarsOff()

      regions = list()

      for r in range(connectivityFilter.GetNumberOfExtractedRegions()):

        selectorCon.ThresholdBetween(r,r)
        geometryCon = vtk.vtkGeometryFilter()
        geometryCon.SetInputConnection(selectorCon.GetOutputPort())
        geometryCon.Update()
        cleanFilter = vtk.vtkCleanPolyData()
        cleanFilter.SetInputConnection(geometryCon.GetOutputPort())
        cleanFilter.Update()

        pd = vtk.vtkPolyData()
        pd.DeepCopy(cleanFilter.GetOutput())
        regions.append(pd)
      
      return regions

  @staticmethod
  def _transformPD(polydata, trf):
    transformPD = vtk.vtkTransformPolyDataFilter()
    transformPD.SetTransform(trf)
    transformPD.SetInputData(polydata)
    transformPD.Update()

    return transformPD.GetOutput()

  @staticmethod
  def _regLandmark(sourceLM, targetLM):

    def fiducialsToPoints(fiducials):
      points = vtk.vtkPoints()
      for fid in range(fiducials.GetNumberOfFiducials()):
        pos = [0]*3
        fiducials.GetNthFiducialPosition(fid, pos)
        points.InsertNextPoint(list(pos))
      return points
    
    targetPoints = fiducialsToPoints(targetLM)
    sourcePoints = fiducialsToPoints(sourceLM)

    landmarkTransform = vtk.vtkLandmarkTransform()
    landmarkTransform.SetSourceLandmarks(sourcePoints)
    landmarkTransform.SetTargetLandmarks(targetPoints)
    landmarkTransform.SetModeToRigidBody()
    landmarkTransform.Modified()
    landmarkTransform.Update()

    return landmarkTransform

  @staticmethod
  def _regICP(sourcePD, targetPD):

    icpTransform = vtk.vtkIterativeClosestPointTransform()
    icpTransform.SetSource(sourcePD)
    icpTransform.SetTarget(targetPD)
    icpTransform.GetLandmarkTransform().SetModeToRigidBody()
    icpTransform.StartByMatchingCentroidsOff()
    icpTransform.Modified()
    icpTransform.Update()

    return icpTransform

  @staticmethod
  def _calculateDistances(sourcePD, targetPD):
    dist = vtk.vtkDistancePolyDataFilter()
    dist.AddInputData(0, sourcePD)
    dist.AddInputData(1, targetPD)
    dist.SignedDistanceOff()
    dist.NegateDistanceOff()
    dist.ComputeSecondDistanceOff()
    dist.Update()
    
    return vtk.util.numpy_support.vtk_to_numpy(dist.GetOutput().GetPointData().GetArray("Distance"))

  @staticmethod
  def _regionMostCommonIds(regions, idArrayName, ids, minSize=0):
    for i in range(regions[0].GetPointData().GetNumberOfArrays()):
      if regions[0].GetPointData().GetArrayName(i) == idArrayName:
        break

    maxIntersects = -1
    r = -1
    for p, pd in enumerate(regions):
      ar = vtk.util.numpy_support.vtk_to_numpy(pd.GetPointData().GetArray(i))
  
      if len(np.intersect1d(ar,ids)) > 0 and len(ar) >= minSize:
        maxIntersects = len(np.intersect1d(ar,ids))
        r = p
    
    
    if r == -1:
      raise IndexError("No common region found.")
    
    return r

  @staticmethod
  def _getSurfaceArea(polydata):
    mass = vtk.vtkMassProperties()
    mass.SetInputData(polydata)
    mass.Update() 

    return mass.GetSurfaceArea()


  # Helper
  @staticmethod
  def _addTrf(trf, name=None):
    no = slicer.vtkMRMLLinearTransformNode()
    no.SetAndObserveTransformToParent(trf)
    t = slicer.mrmlScene.AddNode(no)
    if name:
      t.SetName(name)
    return t
  
  @staticmethod
  def _addArray(polydata, values, name):
    arr = vtk.util.numpy_support.numpy_to_vtk(values)
    arr.SetName(name)
    polydata.GetPointData().RemoveArray(name)
    return polydata.GetPointData().AddArray(arr)

  @staticmethod
  def _multiplyTransforms(trfList):
    def mult(trf1, trf2):
      dotMat = vtk.vtkMatrix4x4()
      vtk.vtkMatrix4x4().Multiply4x4(trf2.GetMatrix(), trf1.GetMatrix(), dotMat)
      trf = vtk.vtkTransform()
      trf.SetMatrix(dotMat)
      return trf
    
    outTrf = trfList[0]
    for t in range(len(trfList)-1):
      outTrf = mult(outTrf, trfList[t+1])
    
    return outTrf

  @staticmethod
  def getArrayNumber(polydata, name):
    for i in range(polydata.GetPointData().GetNumberOfArrays()):
      if polydata.GetPointData().GetArrayName(i) == name:
        return i

    raise NameError('Array name not found in passed polydata')

#
# SurfaceFragmentsRegistrationTest
#

class SurfaceFragmentsRegistrationTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear()
    self.logic = SurfaceFragmentsRegistrationLogic()

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    #self.test_SurfaceFragmentsRegistration1()
    self.test_SurfaceFragmentsRegistration2()

  def test_SurfaceFragmentsRegistration1(self):
    """Minimal working scenario with spheres.
    """

    self.delayDisplay("Starting test_SimilaritySubgroups1")

    # from SegmentStatistics import SegmentStatisticsLogic

    self.delayDisplay("Create models and fiducials")

    sphere1 = vtk.vtkSphereSource()
    sphere1.SetRadius(30)
    sphere1.SetCenter(0,0,0)
    sphere1.SetPhiResolution(50)
    sphere1.SetThetaResolution(50)
    sphere1.Update()

    sphere2 = vtk.vtkSphereSource()
    sphere2.SetRadius(30)
    sphere2.SetCenter(10,0,0)
    sphere2.SetPhiResolution(50)
    sphere2.SetThetaResolution(50)
    sphere2.Update()

    sphere3 = vtk.vtkSphereSource()
    sphere3.SetRadius(30)
    sphere3.SetCenter(50,50,50)
    sphere3.SetPhiResolution(50)
    sphere3.SetThetaResolution(50)
    sphere3.Update()

    boolean = vtk.vtkBooleanOperationPolyDataFilter()
    boolean.SetOperationToUnion()
    boolean.SetInputData(0, sphere1.GetOutput())
    boolean.SetInputData(1, sphere2.GetOutput())
    boolean.Update()

    smod = slicer.modules.models.logic().AddModel(boolean.GetOutput())
    tmod = slicer.modules.models.logic().AddModel(sphere3.GetOutput())

    sfid = slicer.vtkMRMLMarkupsFiducialNode()
    slicer.mrmlScene.AddNode(sfid)
    sfid.AddFiducial(0,30,0)
    sfid.AddFiducial(0,-30,0)
    sfid.AddFiducial(-30,0,0)

    tfid = slicer.vtkMRMLMarkupsFiducialNode()
    slicer.mrmlScene.AddNode(tfid)
    tfid.AddFiducial(50,80,50)
    tfid.AddFiducial(50,20,50)
    tfid.AddFiducial(20,50,50)

    self.delayDisplay("Run algorithm")

    parameterNode = self.logic.getParameterNode()
    self.logic.setDefaultParameters(parameterNode, overwrite=True)
    parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, smod.GetID())
    parameterNode.SetNodeReferenceID(PARAMETER_TARGETMODEL, tmod.GetID())
    parameterNode.SetNodeReferenceID(PARAMETER_SOURCELANDMARKS, sfid.GetID())
    parameterNode.SetNodeReferenceID(PARAMETER_TARGETLANDMARKS, tfid.GetID())
    parameterNode.SetParameter(PARAMETER_INITIALIZATIONCLUSTERRADIUS, str(20))
    parameterNode.SetParameter(PARAMETER_MINIMALCLUSTERAREA, str(10))
    parameterNode.SetParameter(PARAMETER_CUTOFFTHRESHOLD, str(0.5))
    parameterNode.SetParameter(PARAMETER_PREREGISTRATION, "true")
    parameterNode.SetParameter(PARAMETER_INITIALIZATIONITERATIONS, str(5))
    parameterNode.SetParameter(PARAMETER_OPENINGWIDTH, str(1))
    parameterNode.SetParameter(PARAMETER_MAXIMALITERATIONS, str(10))
    parameterNode.SetParameter(PARAMETER_ADDTRANSFORMATION, "true")
    parameterNode.SetParameter(PARAMETER_MARKDEVIATION, "true")
    parameterNode.SetParameter(PARAMETER_MARKFRAGMENTS, "true")

    self.logic.process(parameterNode)

    self.delayDisplay("Evaluate Result")

    arrays = list()
    for a in range(smod.GetPolyData().GetPointData().GetNumberOfArrays()):
      ar = smod.GetPolyData().GetPointData().GetArray(a)
      if ar.GetName().startswith(PREFIX_DEVIATION):
        arrays.append(vtk.util.numpy_support.vtk_to_numpy(ar))
    
    clusterDev = np.min(arrays, axis=0)

    meanMinDev = np.mean(clusterDev)
    self.assertLessEqual(meanMinDev, 0.5)

    self.delayDisplay('Test passed')

  def test_SurfaceFragmentsRegistration2(self):
    """Real world scenario with segmented hemipelvis.
    """

    self.delayDisplay("Starting test_SimilaritySubgroups2")

    # Get/create input data

    self.delayDisplay('Load Test Data Set')
    import SampleData
    registerSampleData()
    SampleData.downloadSamples('CTAAbdomenPanoramix')
    nodes = SampleData.downloadSamples('Surface Fragments Registration Models')
    sourceModel, targetModel = nodes
    sourceModel.GetDisplayNode().SetSliceIntersectionVisibility(1)
    targetModel.GetDisplayNode().SetSliceIntersectionVisibility(1)
    sourceModel.GetDisplayNode().SetSliceIntersectionThickness(2)

    self.delayDisplay("Run algorithm")
    parameterNode = self.logic.getParameterNode()
    self.logic.setDefaultParameters(parameterNode, overwrite=True)
    parameterNode.SetNodeReferenceID(PARAMETER_SOURCEMODEL, sourceModel.GetID())
    parameterNode.SetNodeReferenceID(PARAMETER_TARGETMODEL, targetModel.GetID())
    parameterNode.SetParameter(PARAMETER_OPENINGWIDTH, str(10))
    parameterNode.SetParameter(PARAMETER_CUTOFFTHRESHOLD, str(1.0))

    self.logic.process(parameterNode)

    self.delayDisplay("Evaluate Result")

    arrays = list()
    for a in range(sourceModel.GetPolyData().GetPointData().GetNumberOfArrays()):
      ar = sourceModel.GetPolyData().GetPointData().GetArray(a)
      if ar.GetName().startswith(PREFIX_DEVIATION):
        arrays.append(vtk.util.numpy_support.vtk_to_numpy(ar))
    
    clusterDev = np.min(arrays, axis=0)

    meanMinDev = np.mean(clusterDev)
    self.assertLessEqual(meanMinDev, 0.5)

    self.delayDisplay('Test passed')


DEFAULT_INITIALIZATIONCLUSTERRADIUS = 50
DEFAULT_INITIALIZATIONITERATIONS = 3
DEFAULT_MINIMALCLUSTERAREA = 500
DEFAULT_OPENINGWIDTH = 8
DEFAULT_CUTOFFTHRESHOLD = 1.5
DEFAULT_MAXIMALITERATIONS = 10
DEFAULT_PREREGISTRATION = False
DEFAULT_MARKDEVIATION = True
DEFAULT_ADDTRANSFORMATION = True
DEFAULT_MARKFRAGMENTS = False

PARAMETER_SOURCEMODEL = "SourceModel"
PARAMETER_TARGETMODEL = "TargetModel"
PARAMETER_SOURCELANDMARKS = "SourceLandmarks"
PARAMETER_TARGETLANDMARKS = "TargetLandmarks"
PARAMETER_INITIALIZATIONCLUSTERRADIUS = "InitializationClusterRadius"
PARAMETER_INITIALIZATIONITERATIONS = "InitializationIterations"
PARAMETER_MINIMALCLUSTERAREA = "MinimalClusterArea"
PARAMETER_OPENINGWIDTH = "OpeningWidth"
PARAMETER_CUTOFFTHRESHOLD = "CutoffThreshold"
PARAMETER_MAXIMALITERATIONS = "MaximalIterations"
PARAMETER_PREREGISTRATION = "Preregistration"
PARAMETER_MARKDEVIATION = "MarkDeviation"
PARAMETER_ADDTRANSFORMATION = "AddTransformation"
PARAMETER_MARKFRAGMENTS = "MarkFragments"

PREFIX_FRAGMENT = "sf_nr-"
PREFIX_DEVIATION = "sf_deviation-"
PREFIX_TRANSFORMATION = "sf_trf-"

LET_INTERMEDIATEDISTANCES = "intermediate-distances"
LET_INTERMEDIATEREGISTERED = "intermediate-registered"
LET_INTERMEDIATEIDS = "intermediate-ids"