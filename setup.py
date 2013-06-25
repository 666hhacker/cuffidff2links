from distutils.sysconfig import get_python_lib
from shutil import copy
import os

installPath = get_python_lib()
if os.path.exists(os.path.join(installPath,'cuffdiff2links.py')):
   os.remove(os.path.join(installPath,'cuffdiff2links.py'))
if os.path.exists(os.path.join(installPath,'cuffdiff2links.pyc')):
   os.remove(os.path.join(installPath,'cuffdiff2links.pyc'))
copy('cuffdiff2links.py',get_python_lib())
