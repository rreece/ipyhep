"""
    ipyhep - Python package for making ROOT plots

author: Ryan D. Reece  <ryan.reece@cern.ch>

TODO: description.

See also:
    http://www.hep.upenn.edu/~rreece/computing.html
    http://root.cern.ch/
    http://root.cern.ch/root/HowtoPyROOT.html
"""

#------------------------------------------------------------------------------
# module metadata
#------------------------------------------------------------------------------
__author__ = 'Ryan D. Reece'
__email__ = 'ryan.reece@cern.ch'
__copyright__ = 'Copyright 2015-2016 Ryan D. Reece'
__license__ = 'GPL http://www.gnu.org/licenses/gpl.html'

## load ROOT
import ROOT
#from . import rootlogon
#from . import rootnotes
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1001

# turn off warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

from rootpy import log
#log.setLevel(log.WARNING)
log.setLevel(log.ERROR)
from rootpy.plotting import set_style
set_style('ATLAS')
#set_style('default')
#set_style('cmstdr')

## modules in this package:
from . import file
from . import hist
from . import plot
#from . import poissonize
from . import sampleops
from . import stats
from . import style
from . import table
from . import tree

__all__ = ['file', 'hist', 'plot', 'sampleops', 'stats', 'style', 'table', 'tree']

