"""
TODO: write docstring.
"""

import glob
import math

import ROOT

#from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D
from rootpy.tree import TreeChain, Cut

import ipyhep


#------------------------------------------------------------------------------
# Functions for the user:
#
# -   project(x, *args, **kwargs):
#
#------------------------------------------------------------------------------

#______________________________________________________________________________
def project(x, *args, **kwargs):
    """
    x = '/path/to/file.root::mytree::varx:vary'
    """
    
    f = None
    t = None
    var = None
    h = None
    filepath = None
    treename = None

    assert x.count('::') == 2
    filepath, treename, var = x.split('::')
    
    assert filepath
    assert treename
    assert var

    if filepath.count('*'):
        paths = glob.glob(filepath)
        paths.sort()
        if paths:
            t = ROOT.TChain(treename)
            for p in paths:
                ## check that treename exists in path p
                tmp_f = ROOT.TFile.Open(p)
                tmp_keys = [k.GetName() for k in tmp_f.GetListOfKeys()]
                tmp_f.Close()
                if treename in tmp_keys:
                    t.Add(p)
    else:
        assert False
            
#### rootpy stuff gave bad Draw
##        assert len(paths) > 0
#        if len(paths) == 1:
#            filepath = paths[0]
#            ## get tree as rootpy tree
#            f = root_open(filepath)
#            t = getattr(f, treename)
#        if len(paths) > 1:
##            t = TreeChain(treename, paths)
#            t = ROOT.TChain(treename)
#            for p in paths:
#                t.Add(p)
#        else:
#            pass ## allow for empty samples
##        print 'Reading %i files.' % len(paths)

    if t:
        selection = kwargs.pop('selection', '')
        weight = kwargs.pop('weight', 1.0)
        bins = kwargs.pop('bins', None)
        includeover = kwargs.pop('includeover', False)
        h = _project(t, var, selection, weight, bins, includeover)

    if f:
        f.close()

    return h


## only the above functions are exposed
__all__ = ['project']


#------------------------------------------------------------------------------
# Private helper functions
#------------------------------------------------------------------------------

#______________________________________________________________________________
def _project(tree, var, selection='', weight=1.0, bins=None, includeover=False):

    h = None
    if var.count(':') == 0:
        ## Hist (1D)
        if bins:
            if isinstance(bins, tuple):
                assert len(bins) == 3
                h = Hist(*bins)
            elif isinstance(bins, list):
                h = Hist(bins)
            else:
                assert False
        else:
            assert False
    elif var.count(':') == 1:
        ## Hist2D
        ## like rootpy, we use a convention where var=x:y, unlike ROOT
        varx, vary = var.split(':')
        var = ':'.join([vary, varx])
        if bins:
            if isinstance(bins, tuple):
                assert len(bins) == 6
                h = Hist2D(*bins)
            elif isinstance(bins, list):
                ## TODO: support variable bins for Hist2D
                h = Hist2D(*bins)
                #assert False
            else:
                assert False
        else:
             assert False
    else:
        assert False

    assert h
#            kwargs['hist'] = h

    ## add the weight to the selection via a TCut
    weighted_selection = str(selection)
    if weight and weight != 1.0:
        weighted_selection = Cut(weighted_selection)
        weighted_selection = weighted_selection * weight
        weighted_selection = str(weighted_selection)

    tree.Draw('%s>>%s' % (var, h.GetName()), weighted_selection)

##    print tree.GetSelectedRows() ## debuging

#    for event in t:
#        x = getattr(event, var)
#        h.fill(x)

    if h:
        h.SetDirectory(0)

#        elist = ROOT.gDirectory.Get('elist')
#        elist.Print('all')

        if var.count(':') == 0 and includeover:
            ## include overflow for Hist (1D)
            nbins = h.GetNbinsX()
            c1 = h.GetBinContent(nbins)
            c2 = h.GetBinContent(nbins+1)
            e1 = h.GetBinError(nbins)
            e2 = h.GetBinError(nbins+1)
            h.SetBinContent(nbins, c1+c2)
            h.SetBinError(nbins, math.sqrt( (c1*e1*e1 + c2*e2*e2)/(c1+c2) ) if c1+c2!=0.0 else 0.0)
            h.SetBinContent(nbins+1, 0.0)
            h.SetBinError(nbins+1, 0.0)

    return h


