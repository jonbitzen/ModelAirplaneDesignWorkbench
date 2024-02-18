import os
import FreeCADGui as Gui
import FreeCAD as App
from freecad.model_airplane_design import ICONPATH


class ModelAirplaneDesignWorkbench(Gui.Workbench):
    """
    class which gets initiated at startup of the gui
    """

    MenuText = "Model Airplane Design"
    ToolTip = "A workbench for designing model airplanes"
    Icon = os.path.join(ICONPATH, "icons8-airplane-64.png")
    toolbox = []

    def GetClassName(self):
        return "Gui::PythonWorkbench"

    def Initialize(self):
        """
        This function is called at the first activation of the workbench.
        here is the place to import all the commands
        """

        self.appendToolbar("Tools", self.toolbox)
        self.appendMenu("Tools", self.toolbox)

    def Activated(self):
        pass

    def Deactivated(self):
        pass


Gui.addWorkbench(ModelAirplaneDesignWorkbench())
