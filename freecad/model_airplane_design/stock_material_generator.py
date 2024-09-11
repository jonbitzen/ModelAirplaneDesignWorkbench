import FreeCAD as App
import numpy
import Part
from typing import List

def imp_frac(inches: float):
    return inches * 25.4

def create_cuboid(
        obj_name: str, 
        width: float = imp_frac(1/4), 
        height: float = imp_frac(1/4)
    ) -> App.DocumentObject:
    """
        Creates a cuboid stock object

        Parameters
        ----------
        obj_name: str
            name for the new object
        width: float
            width for the new object
        height: float
            height for the new object

        Return
        ------
        PartDesign.Body
            A body whose first feature is a cuboid solid
    """
    body = App.ActiveDocument.addObject("PartDesign::Body", obj_name)
    cuboid_stock = App.ActiveDocument.addObject("PartDesign::FeaturePython", "cuboid_base")
    StockCuboid(cuboid_stock, width, height)
    StockCuboidViewProvider(cuboid_stock.ViewObject)
    body.addObject(cuboid_stock)

def create_tube(
        obj_name: str,
        outer_diameter: float = 5.0,
        inner_diameter: float = 3.0
        ) -> App.DocumentObject:
    """
        Creates a tube stock object 

        Parameters
        ----------
        outer_diameter: float
            tube outer diameter
        inner_diameter: float
            tube inner diameter

        Return
        ------
        PartDesign.Body
            A body whose first feature is a tubular solid
    """
    body = App.ActiveDocument.addObject("PartDesign::Body", obj_name)
    tube = App.ActiveDocument.addObject("PartDesign::FeaturePython","tube_base")
    StockTube(tube, outer_diameter, inner_diameter)
    StockTubeViewProvider(tube.ViewObject)
    body.addObject(tube)

def create_cylinder(
        obj_name: str,
        diameter: float = imp_frac(1/4)
    ) -> App.DocumentObject:
    """
        Creates a cylindrical stock object

        Parameters
        ----------
        diameter: float
            cylinder diameter

        Return
        ------
        PartDesign.Body
            A body whose first feature is a cylindrical solid
    """
    body = App.ActiveDocument.addObject("PartDesign::Body", obj_name)
    tube = App.ActiveDocument.addObject("PartDesign::FeaturePython","cylinder_base")
    StockCylinder(tube, diameter)
    StockCylinderViewProvider(tube.ViewObject)
    body.addObject(tube)

class StockCuboid:
    def __init__(
            self,
            obj: App.DocumentObject,
            width: float,
            height: float
        ) -> None:
        
        obj.addProperty(
            "App::PropertyLength",
            "width",
            "Stock Cuboid",
            "Cuboid stock front profile width"
        ).width = width

        obj.addProperty(
            "App::PropertyLength",
            "height",
            "Stock Cuboid",
            "Cuboid stock front profile height"
        ).height = height

        obj.addProperty(
            "App::PropertyLength",
            "length",
            "Stock Cuboid",
            "Cuboid stock material length"
        ).length = 50.0

        self.makeAttachable(obj)
        obj.Proxy = self

    def makeAttachable(self, obj: App.DocumentObject):
        obj.addExtension('Part::AttachExtensionPython')
        obj.setEditorMode('Placement', 0)

    def execute(self, obj: App.DocumentObject):
        
        cuboid = Part.makeBox(obj.length, obj.width, obj.height)
        
        mtx = App.Matrix()
        mtx.move(-cuboid.BoundBox.Center)
        cuboid = cuboid.transformed(mtx, True)
        
        if not hasattr(obj, "positionBySupport"):
            self.makeAttachable(obj)
        obj.positionBySupport()
        cuboid.Placement = obj.Placement

        obj.Shape = cuboid

        # BaseFeature (shape property of type Part::PropertyPartShape) is provided
        # for uswith the PartDesign::FeaturePython and related classes, but it
        # might be empty if our object is the first object in the tree.  it's a
        # good idea to checkfor its existence in case we want to make type 
        # Part::FeaturePython, which won't have it
        if hasattr(obj, "BaseFeature") and obj.BaseFeature != None:
            if "Subtractive" in obj.TypeId:
                full_shape = obj.BaseFeature.Shape.cut(cuboid)
            else:
                full_shape = obj.BaseFeature.Shape.fuse(cuboid)
            full_shape.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.Shape = full_shape
        else:
            obj.Shape = cuboid

        # PartDesign::FeatureAdditivePython and PartDesign::FeatureSubtractivePython 
        # have this property but PartDesign::FeaturePython does not, it is the 
        # shape used for copying in pattern features for example in making a 
        # polar pattern
        if hasattr(obj,"AddSubShape"):
            cuboid.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.AddSubShape = cuboid

class StockCuboidViewProvider:
    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject):
        vobj.Proxy = self
        self.Object = vobj.Object

    def attach(self, vobj: App.Gui.ViewProviderDocumentObject):
        self.vobj = vobj

    def updateData(self, obj: App.DocumentObject, property: str):
        pass

    def getDisplayModes(self, vobj: App.Gui.ViewProviderDocumentObject) -> List[str]:
        modes=[]
        modes.append("Flat Lines")
        modes.append("Shaded")
        modes.append("Wireframe")
        return modes
    
    def getDefaultDisplayMode(self) -> str:
        return "Flat Lines"

    def setDisplayMode(self, mode: str) -> str:
        return mode
    
    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, property: str):
        pass

    def dumps(self) -> str:
        return None
    
    def loads(self, state):
        return None
    
class StockTube:
    def __init__(
            self,
            obj: App.DocumentObject,
            outer_diameter: float,
            inner_diameter: float
        ) -> None:

        obj.addProperty(
            "App::PropertyLength",
            "inner_diameter",
            "Tube",
            "inner_diameter"
        ).inner_diameter = inner_diameter

        obj.addProperty(
            "App::PropertyLength",
            "outer_diameter",
            "Tube",
            "outer_diameter"
        ).outer_diameter = outer_diameter

        obj.addProperty(
            "App::PropertyLength",
            "length",
            "Tube",
            "Length of the stock tube"
        ).length = 10

        self.makeAttachable(obj)
        obj.Proxy = self

    def makeAttachable(self, obj: App.DocumentObject):
        obj.addExtension('Part::AttachExtensionPython')
        obj.setEditorMode('Placement', 0) 

    def execute(self,  obj: App.DocumentObject):
        outer_cylinder = Part.makeCylinder(obj.outer_diameter/2, obj.length)
        inner_cylinder = Part.makeCylinder(obj.inner_diameter/2, obj.length)
        # just make cylinder
        if obj.inner_diameter == obj.outer_diameter: 
            tube_shape = outer_cylinder
        elif obj.inner_diameter < obj.outer_diameter:
            tube_shape = outer_cylinder.cut(inner_cylinder)
        else: #invert rather than error out
            tube_shape = inner_cylinder.cut(outer_cylinder)

        mtx = App.Matrix()
        mtx.move(-tube_shape.BoundBox.Center)
        mtx.rotateY(numpy.deg2rad(-90))
        tube_shape = tube_shape.transformed(mtx, True)

        if not hasattr(obj, "positionBySupport"):
            self.makeAttachable(obj)
        obj.positionBySupport()
        tube_shape.Placement = obj.Placement

        # TODO make the shape body centered so its easier to manage

        # BaseFeature (shape property of type Part::PropertyPartShape) is provided
        # for uswith the PartDesign::FeaturePython and related classes, but it
        # might be empty if our object is the first object in the tree.  it's a
        # good idea to checkfor its existence in case we want to make type 
        # Part::FeaturePython, which won't have it
        if hasattr(obj, "BaseFeature") and obj.BaseFeature != None:
            if "Subtractive" in obj.TypeId:
                full_shape = obj.BaseFeature.Shape.cut(tube_shape)
            else:
                full_shape = obj.BaseFeature.Shape.fuse(tube_shape)
            full_shape.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.Shape = full_shape
        else:
            obj.Shape = tube_shape

        # PartDesign::FeatureAdditivePython and PartDesign::FeatureSubtractivePython 
        # have this property but PartDesign::FeaturePython does not, it is the 
        # shape used for copying in pattern features for example in making a 
        # polar pattern
        if hasattr(obj,"AddSubShape"):
            tube_shape.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.AddSubShape = tube_shape

class StockTubeViewProvider:

    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject):
        '''Set this object to the proxy object of the actual view provider'''
        vobj.Proxy = self

    def attach(self, vobj: App.Gui.ViewProviderDocumentObject):
        self.vobj = vobj

    def updateData(self, obj: App.DocumentObject, property: str):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        pass

    def getDisplayModes(self,vobj: App.Gui.ViewProviderDocumentObject) -> List[str]:
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Flat Lines")
        modes.append("Shaded")
        modes.append("Wireframe")
        return modes

    def getDefaultDisplayMode(self) -> str:
        '''Return the name of the default display mode. It must be defined in getDisplayModes.'''
        return "Flat Lines"

    def setDisplayMode(self,mode) -> str:
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode

    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, property: str):
        '''Here we can do something when a single property got changed'''
        #App.Console.PrintMessage("Change property: " + str(prop) + "\n")
        pass

    def getIcon(self):
        '''Return the icon in XPM format which will appear in the tree view. This method is\
                optional and if not defined a default icon is shown.'''
        return """
            /* XPM */
            static const char * ViewProviderBox_xpm[] = {
            "16 16 6 1",
            "   c None",
            ".  c #141010",
            "+  c #615BD2",
            "@  c #C39D55",
            "#  c #000000",
            "$  c #57C355",
            "        ........",
            "   ......++..+..",
            "   .@@@@.++..++.",
            "   .@@@@.++..++.",
            "   .@@  .++++++.",
            "  ..@@  .++..++.",
            "###@@@@ .++..++.",
            "##$.@@$#.++++++.",
            "#$#$.$$$........",
            "#$$#######      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            " #$#$$$$$#      ",
            "  ##$$$$$#      ",
            "   #######      "};
            """

    def dumps(self) -> str:
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None

    def loads(self, state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None
    
class StockCylinder:
    def __init__(
            self,
            obj: App.DocumentObject,
            diameter: float
        ) -> None:

        obj.addProperty(
            "App::PropertyLength",
            "diameter",
            "Tube",
            "diameter"
        ).diameter = diameter

        obj.addProperty(
            "App::PropertyLength",
            "length",
            "Tube",
            "Length of the stock cylinder"
        ).length = 10

        self.makeAttachable(obj)
        obj.Proxy = self

    def makeAttachable(self, obj: App.DocumentObject):
        obj.addExtension('Part::AttachExtensionPython')
        obj.setEditorMode('Placement', 0) 

    def execute(self,  obj: App.DocumentObject):
        cylinder = Part.makeCylinder(obj.diameter/2, obj.length)
    
        mtx = App.Matrix()
        mtx.move(-cylinder.BoundBox.Center)
        mtx.rotateY(numpy.deg2rad(-90))
        cylinder = cylinder.transformed(mtx, True)

        if not hasattr(obj, "positionBySupport"):
            self.makeAttachable(obj)
        obj.positionBySupport()
        cylinder.Placement = obj.Placement

        # TODO make the shape body centered so its easier to manage

        # BaseFeature (shape property of type Part::PropertyPartShape) is provided
        # for uswith the PartDesign::FeaturePython and related classes, but it
        # might be empty if our object is the first object in the tree.  it's a
        # good idea to checkfor its existence in case we want to make type 
        # Part::FeaturePython, which won't have it
        if hasattr(obj, "BaseFeature") and obj.BaseFeature != None:
            if "Subtractive" in obj.TypeId:
                full_shape = obj.BaseFeature.Shape.cut(cylinder)
            else:
                full_shape = obj.BaseFeature.Shape.fuse(cylinder)
            full_shape.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.Shape = full_shape
        else:
            obj.Shape = cylinder

        # PartDesign::FeatureAdditivePython and PartDesign::FeatureSubtractivePython 
        # have this property but PartDesign::FeaturePython does not, it is the 
        # shape used for copying in pattern features for example in making a 
        # polar pattern
        if hasattr(obj,"AddSubShape"):
            cylinder.transformShape(obj.Placement.inverse().toMatrix(), True)
            obj.AddSubShape = cylinder

class StockCylinderViewProvider:

    def __init__(self, vobj: App.Gui.ViewProviderDocumentObject):
        '''Set this object to the proxy object of the actual view provider'''
        vobj.Proxy = self

    def attach(self, vobj: App.Gui.ViewProviderDocumentObject):
        self.vobj = vobj

    def updateData(self, obj: App.DocumentObject, property: str):
        '''If a property of the handled feature has changed we have the chance to handle this here'''
        pass

    def getDisplayModes(self,vobj: App.Gui.ViewProviderDocumentObject) -> List[str]:
        '''Return a list of display modes.'''
        modes=[]
        modes.append("Flat Lines")
        modes.append("Shaded")
        modes.append("Wireframe")
        return modes

    def getDefaultDisplayMode(self) -> str:
        '''Return the name of the default display mode. It must be defined in getDisplayModes.'''
        return "Flat Lines"

    def setDisplayMode(self,mode) -> str:
        '''Map the display mode defined in attach with those defined in getDisplayModes.\
                Since they have the same names nothing needs to be done. This method is optional'''
        return mode

    def onChanged(self, vobj: App.Gui.ViewProviderDocumentObject, property: str):
        '''Here we can do something when a single property got changed'''
        #App.Console.PrintMessage("Change property: " + str(prop) + "\n")
        pass

    def getIcon(self):
        '''Return the icon in XPM format which will appear in the tree view. This method is\
                optional and if not defined a default icon is shown.'''
        return """
            /* XPM */
            static const char * ViewProviderBox_xpm[] = {
            "16 16 6 1",
            "   c None",
            ".  c #141010",
            "+  c #615BD2",
            "@  c #C39D55",
            "#  c #000000",
            "$  c #57C355",
            "        ........",
            "   ......++..+..",
            "   .@@@@.++..++.",
            "   .@@@@.++..++.",
            "   .@@  .++++++.",
            "  ..@@  .++..++.",
            "###@@@@ .++..++.",
            "##$.@@$#.++++++.",
            "#$#$.$$$........",
            "#$$#######      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            "#$$#$$$$$#      ",
            " #$#$$$$$#      ",
            "  ##$$$$$#      ",
            "   #######      "};
            """

    def dumps(self) -> str:
        '''When saving the document this object gets stored using Python's json module.\
                Since we have some un-serializable parts here -- the Coin stuff -- we must define this method\
                to return a tuple of all serializable objects or None.'''
        return None

    def loads(self, state):
        '''When restoring the serialized object from document we have the chance to set some internals here.\
                Since no data were serialized nothing needs to be done here.'''
        return None
    