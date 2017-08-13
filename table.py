"""
TODO: write docstring.
"""

## std
import math
import copy

## ROOT
import ROOT

## mine
import ipyhep
from sigdigs import sigdigs
import tabletools


#------------------------------------------------------------------------------
# Functions for the user:
#
#------------------------------------------------------------------------------

#______________________________________________________________________________
def make_integrated_row(hists):
    row = list()
    for h in hists:
        count, err = ipyhep.hist.calc_integral_and_error(h)
        row.append( (count, err) )
    return row


#______________________________________________________________________________
def make_integrated_cutflow(hists, cuts):
    tab = list()
    for xcut in cuts:
        row = list()
        for h in hists:
            count, err = ipyhep.hist.calc_integral_and_error(h, xcut)
            row.append( (count, err) )
        tab.append(row)
    return tab


#______________________________________________________________________________
def make_rounded_str_table(tab, collables, rowlables):
    assert len(tab) == len(rowlables)
    tabout = list()
    rowout = [''] + collables
    tabout.append(rowout)
    for row, xcut in zip(tab, rowlables):
        assert len(row) == len(collables)
        rowout = list()
        if isinstance(xcut, float):
            rowout.append('%i' % round(xcut)) # NOTE: assumes xcuts falls on integers
        else:
            rowout.append(xcut)
        for count, err in row:
            s = make_rounded_str(count, err)
            rowout.append(s)
        tabout.append(rowout)
    print tabletools.ascii_table.make_str(tabout)
    return tabout


#______________________________________________________________________________
def make_rounded_str(x, ex):
    s = ''
    s = sigdigs.round_to_error(x, ex) 
    if s == '0':
        s = ''
    return s


#______________________________________________________________________________
def write_latex_table(tab, tex='table.tex'):
    _tab = copy.deepcopy(tab)
    return tabletools.latex_table.write(_tab, tex)


