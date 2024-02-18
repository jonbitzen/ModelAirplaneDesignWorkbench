from setuptools import setup
import os
# from freecad.workbench_starterkit.version import __version__
# name: this is the name of the distribution.
# Packages using the same name here cannot be installed together

version_path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 
                            "freecad", "model_airplane_design", "version.py")
with open(version_path) as fp:
    exec(fp.read())

setup(name='freecad.model_airplane_design',
      version=str(__version__),
      packages=['freecad',
                'freecad.model_airplane_design'],
      maintainer="Jonathan E. Railsback",
      maintainer_email="jonbitzen@hotmail.com",
      url="",
      description="workbench to design flying model airplanes",
      install_requires=[],
      include_package_data=True)
