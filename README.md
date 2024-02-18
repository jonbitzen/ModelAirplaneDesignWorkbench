## FreeCAD Model Airplane Design Workbench

A workbench with tools to help generate model airplane structures in FreeCAD

## Development

Following are some tips to set up the project and have hot reloading and code assistance to improve development productivity


### Developing a FreeCAD Workbench with Hot Reloading

- set up a new workbench according to the workbench template available on Github
  - be sure to use python modules and packages that use "__init__.py" so that package paths get automatically imported

- to develop the workbench efficiently, softlink the top-level module folder into the FreeCAD user modules directory
  - `ln -s ${WB_DEV_PATH}/MyWorkBench ~/.local/share/FreeCAD/Mod` where "MyWorkBench" is the top-level folder of your workbench, containing all the artifacts that are ingested by FreeCAD to instantiate your workbench
  - when working on your workbench, import your top level workbench module to bring it into FreeCAD.  doing so will allow you to reload it as follows:

```
import my_wb_module

# experiment, make changes to your source, then reload them
from importlib import reload
reload(my_wb_module)
```

### Vscode Type Hinting
- install FreeCAD into the /opt folder as an exploded AppImage archive; we need to install FreeCAD as an exploded archive so that we can obtain Python library paths, and provide them to vscode for python code analysis
  - `./FreeCAD*.AppImage --appimage-extract`
    - this will usually create a default folder name like "squashfs", which you should rename to something more recognizable
- in the workbench folder, create a hidden folder ".vscode", and inside it create a file called "settings.json"
- in the file settings.json, create the following properties
```
{
"python.analysis.stubPath": "FREECAD_STUB_PATH",
"python.analysis.extraPaths": [
    "LIB_PATH1",
    "LIB_PATH2",
    ...
    "LIB_PATHN"
]
}
```
  - in `python.analysis.extraPaths`, fill in the entire python path set from the FreeCAD python console.  You can obtain this by starting FreeCAD, going to the python console, and entering the following:
```
import sys
print(sys.path)
```
  - you can use a text editor to reformat the path strings in a manner suitable for entry into the extraPaths variable above
    - format to get each path on its own line for viewing ease (replace ", " with ", \n" in a text editor)
    - replace single-quite string delimiters (`) with double-quotes (")
    - note that this will also contain the paths to installed workbenches that are usually in the user's local folder; if you add or remove workbenches from within FreeCAD, you may need to update the paths to include the changes

- Some of the FreeCAD work benches are provided by native libraries, and require stub files to enable code assistance.  To generate stubs, do the following:
  - create a folder `freecad-stub-gen`, then change directory inside it
  - checkout the FreeCAD source, `git clone https://github.com/FreeCAD/FreeCAD.git`
    - inside the FreeCAD folder, use git to checkout the branch or tagged version of FreeCAD that corresponds to the version you will be using to develop your workbench
  - checkout the repository freecad-stubs, `git clone https://github.com/ostr00000/freecad-stubs`, and change to the latest available tagged version
    - note that you may need to install Python 3.12 or later.  If your system doesn't have Python 3.12, you can install it from the deadsnakes PPA on Ubuntu
  - follow the directions in the freecad-stubs README to generate stubs for the native FreeCAD modules
    - Note, that you may need to add `freecad-stubs/lib` to your PYTHONPATH - this is a step that is omitted in the README directions in the repository
  - after generating the stubs, copy it to a sub-folder called `stubs` underneath your exploded AppImage FreeCAD installation; then provide the new stub directory (e.g., `/opt/freecad-0.21.2-appimage/stubs`) as the `FREECAD_STUB_PATH` variable noted above
    - this is a convenience, so that the generated stubs are available in a single place for other python workbench or macro projects
 
- in vscode, set the interpreter path to the python binary contained in your exploded AppImage archive; while not strictly necessary, this ensures the greatest consistency with the configured libraries for type hinting