"""
TODO: write docstring.
"""

import glob
import math
import os
import time

import ROOT

#from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D, HistStack, Graph, Canvas, Legend
from rootpy.plotting.utils import draw, get_limits

import ipyhep

## the normalization of your ntuples:
__ntuple_lumi = 1.0 # in inverse-pb
__default_lumi = 1000.0 # in inverse-pb
#__use_poissonize = True

## globals 
data = None
bkgs = None
sigs = None
treename = 'mytree'
datasearchpath = '/export/share/gauss/reece/TNT/data16_13TeV_n01/data16_13TeV.*.TNT.*/*.root*'
datadrivensearchpath = '/export/share/gauss/reece/TNT/data16_13TeV_n12/data16_13TeV.*.TNT.*/*.root*'
bkgsearchpath = '/export/share/gauss/reece/TNT/mc15_13TeV_n01/mc15_13TeV.%i.*.TNT.*/*.root*'
sigsearchpath = '/export/share/gauss/reece/TNT/mc15_13TeV_n01/mc15_13TeV.%i.*.TNT.*/*.root*'
lumi = __default_lumi
results = None

dir_of_this_file = os.path.dirname(os.path.abspath( __file__ ))
ROOT.gROOT.Macro( os.path.join(dir_of_this_file, 'helpers.C') )


#------------------------------------------------------------------------------
# Functions for the user:
#
# -   stack(x, *args, **kwargs):
# -   significance_scan_2d(x, bins=(16, 2000., 6000., 10, 100., 300.), xtitle='#font[12]{H}_{#font[132]{T}}  [GeV]', ytitle='#font[12]{E}_{#font[132]{T}}^{#font[132]{miss}}  [GeV]', selection='', weight='', bkgs=None, sig=None, lumi=None, bkgsearchpath=None, sigsearchpath=None, save=None)
#
#------------------------------------------------------------------------------

#______________________________________________________________________________
def set_samples(**kwargs):

    ## parse arguments
    _data = kwargs.pop('data', None)
    _bkgs = kwargs.pop('bkgs', None)
    _sigs = kwargs.pop('sigs', None)
    _treename = kwargs.pop('treename', None)
    _datasearchpath = kwargs.pop('datasearchpath', None)
    _datadrivensearchpath = kwargs.pop('datadrivensearchpath', None)
    _bkgsearchpath = kwargs.pop('bkgsearchpath', None)
    _sigsearchpath = kwargs.pop('sigsearchpath', None)
    _lumi = kwargs.pop('lumi', None)

    global data
    global bkgs
    global sigs
    global treename
    global datasearchpath
    global datadrivensearchpath
    global bkgsearchpath
    global sigsearchpath
    global lumi

    data = _data
    bkgs = _bkgs
    sigs = _sigs
    treename = _treename or treename
    datasearchpath = _datasearchpath or datasearchpath
    datadrivensearchpath = _datadrivensearchpath or datadrivensearchpath
    bkgsearchpath = _bkgsearchpath or bkgsearchpath
    sigsearchpath = _sigsearchpath or sigsearchpath
    if _lumi:
        lumi = float(_lumi)


#______________________________________________________________________________
def stack(x, *args, **kwargs):

    ## parse arguments
    _data = kwargs.pop('data', None)
    _bkgs = kwargs.pop('bkgs', None)
    _sigs = kwargs.pop('sigs', None)
    _treename = kwargs.pop('treename', None)
    _datasearchpath = kwargs.pop('datasearchpath', None)
    _datadrivensearchpath = kwargs.pop('datadrivensearchpath', None)
    _bkgsearchpath = kwargs.pop('bkgsearchpath', None)
    _sigsearchpath = kwargs.pop('sigsearchpath', None)
    _lumi = kwargs.pop('lumi', None)

    global data
    global bkgs
    global sigs
    global treename
    global datasearchpath
    global datadrivensearchpath
    global bkgsearchpath
    global sigsearchpath
    global lumi

    data = _data or data
    bkgs = _bkgs or bkgs
    sigs = _sigs or sigs
    treename = _treename or treename
    datasearchpath = _datasearchpath or datasearchpath
    datadrivensearchpath = _datadrivensearchpath or datadrivensearchpath
    bkgsearchpath = _bkgsearchpath or bkgsearchpath
    sigsearchpath = _sigsearchpath or sigsearchpath
    if _lumi:
        lumi = float(_lumi)

    xtitle = kwargs.pop('xtitle', '')
    ytitle = kwargs.pop('ytitle', '')
    logx   = bool(kwargs.pop('logx', False))
    logy   = bool(kwargs.pop('logy', False))
    blind  = kwargs.pop('blind', None)
    has_blinded_data = False

    ## save stuff to bookkeep and return
    stuff = dict()
    stuff['x'] = x

    ## get data histogram
    h_data = None
    if data:
        sp = datasearchpath # HACK: just data to True!
        newx = '%s::%s::%s' % (sp, treename, x)
        h_data = ipyhep.tree.project(newx, *args, **kwargs)
        if h_data:
            stuff['h_data'] = h_data

    ## blind the data?
    if h_data and not blind is None:
        if isinstance(blind, tuple):
            blind1, blind2 = blind
            nbins = h_data.GetNbinsX()
            for i_bin in xrange(1, nbins+2): # skip underflow (but not overflow)
                xval1 = h_data.GetXaxis().GetBinLowEdge(i_bin)
                xval2 = h_data.GetXaxis().GetBinUpEdge(i_bin)
                if xval1 >= blind1 and xval2 <= blind2:
                    h_data.SetBinContent(i_bin, 0.0)
                    h_data.SetBinError(i_bin, 0.0)
                    has_blinded_data = True
        else:
            nbins = h_data.GetNbinsX()
            for i_bin in xrange(1, nbins+2): # skip underflow (but not overflow)
                xval = h_data.GetXaxis().GetBinLowEdge(i_bin)
                if xval >= blind:
                    h_data.SetBinContent(i_bin, 0.0)
                    h_data.SetBinError(i_bin, 0.0)
                    has_blinded_data = True

    ## get background histograms
    h_bkgs = list()
    n_bkgs = list()
    if bkgs:
        for bkg in bkgs:
            if isinstance(bkg, list):
                h_subtotal = None
                for dsid in bkg:
                    assert isinstance(dsid, str)
                    h_bkg = None
                    if dsid.isdigit():
                        ## mc backgrounds
                        sp = bkgsearchpath % int(dsid)
                        newx = '%s::%s::%s' % (sp, treename, x)
                        h_bkg = ipyhep.tree.project(newx, *args, **kwargs)
                    else:
                        ## data-driven backgrounds
                        assert dsid == 'fakes' or dsid == 'efakes'
                        sp = datadrivensearchpath % dsid
                        newx = '%s::%s::%s' % (sp, treename, x)
                        h_bkg = ipyhep.tree.project(newx, *args, **kwargs)
                    if h_bkg:
                        if h_subtotal:
                            h_subtotal.Add(h_bkg)
                        else:
                            h_subtotal = h_bkg.Clone()
                if h_subtotal:
                    h_bkgs.append(h_subtotal)
                    dsid = bkg[0]
                    n_bkgs.append(dsid)
            else:
                dsid = bkg
                assert isinstance(dsid, str)
                h_bkg = None
                if dsid.isdigit():
                    ## mc backgrounds
                    sp = bkgsearchpath % int(dsid)
                    newx = '%s::%s::%s' % (sp, treename, x)
                    h_bkg = ipyhep.tree.project(newx, *args, **kwargs)
                else:
                    ## data-driven backgrounds
                    assert dsid == 'fakes' or dsid == 'efakes'
                    sp = datadrivensearchpath % dsid
                    newx = '%s::%s::%s' % (sp, treename, x)
                    h_bkg = ipyhep.tree.project(newx, *args, **kwargs)
                if h_bkg:
                    h_bkgs.append(h_bkg)
                    n_bkgs.append(dsid)
        if h_bkgs:
            stuff['h_bkgs'] = h_bkgs

    ## get signal histograms
    h_sigs = list()
    n_sigs = list()
    if sigs:
        for dsid in sigs:
            sp = sigsearchpath % int(dsid)
            newx = '%s::%s::%s' % (sp, treename, x)
            h_sig = ipyhep.tree.project(newx, *args, **kwargs)
            if h_sig:
                h_sigs.append(h_sig)
                n_sigs.append(dsid)
        if h_sigs:
            stuff['h_sigs'] = h_sigs

    assert h_sigs

    ## style data
    if h_data:
        h_data.title = 'Data'
        h_data.linecolor = ipyhep.style.black
        h_data.linewidth = 2
        h_data.markercolor = ipyhep.style.black
        h_data.markerstyle = 20
        h_data.markersize = 1.2
        h_data.fillstyle = ipyhep.style.fill_hollow
        h_data.drawstyle = 'PE'
        h_data.legendstyle = 'LP'

    ## scale and style background histograms
    if h_bkgs:
        assert len(h_bkgs) == len(n_bkgs), '%s\n%s' % (h_bkgs, n_bkgs)

        for h, dsid in zip(h_bkgs, n_bkgs):
            sf = ipyhep.sampleops.get_sf(dsid)
            if dsid.isdigit():
                sf *= lumi/__ntuple_lumi
            h.Scale(sf)

            h.title = ipyhep.sampleops.get_label(dsid)
            h.linecolor = ipyhep.style.black
            h.linewidth = 1
            h.markercolor = ipyhep.sampleops.get_color(dsid)
            h.fillcolor = ipyhep.sampleops.get_color(dsid)
            h.fillstyle = ipyhep.style.fill_solid
            h.legendstyle = 'F'

    ## calculate stat error on total background
    h_bkg_total = None
    if h_bkgs:
        for h_bkg in h_bkgs:
            if h_bkg_total:
                h_bkg_total.Add(h_bkg)
            else:
                h_bkg_total = h_bkg.Clone()
        stuff['h_bkg_total'] = h_bkg_total

    ## style h_bkg_total
    if h_bkg_total:
        h_bkg_total.title = 'stat. uncert.'
        h_bkg_total.linecolor = ipyhep.style.black
        h_bkg_total.linewidth = 1
        h_bkg_total.markerstyle = 0
        h_bkg_total.fillcolor = ipyhep.style.dark_gray
        h_bkg_total.fillstyle = ipyhep.style.fill_lines
        h_bkg_total.drawstyle = 'E2'
        h_bkg_total.legendstyle = 'LF'

    ## scale and style signal histograms
    if h_sigs:
        assert len(h_sigs) == len(n_sigs)
        for h, dsid in zip(h_sigs, n_sigs):
            sf = ipyhep.sampleops.get_sf(dsid)
            sf *= lumi/__ntuple_lumi
            h.Scale(sf)

            h.title = ipyhep.sampleops.get_label(dsid)
            h.linecolor = ipyhep.sampleops.get_color(dsid)
            h.linewidth = 3
            h.fillstyle = ipyhep.style.fill_hollow
            h.markerstyle = 0
            h.drawstyle = 'HIST'
            h.legendstyle = 'L'

    ## build list of all_hists
    all_hists = list()
    main_hists = list()
    if h_data:
        all_hists.append(h_data)
        main_hists.append(h_data)
    if h_bkgs:
        all_hists.extend(h_bkgs)
        main_hists.extend(h_bkgs)
    if h_bkg_total:
        all_hists.append(h_bkg_total)
        main_hists.append(h_bkg_total)
    if h_sigs:
        all_hists.extend(h_sigs)

    ## get statistics
    if all_hists:
        stats_list = list()
        for h in all_hists:
            stats_list.extend( get_stats(h) )
        html = convert_table_to_html( convert_stats_to_table( stats_list ) )
        stuff['html'] = html

    ## renormalize for bin widths
    bins = kwargs.pop('bins', None)
    if bins and isinstance(bins, list):
        for h in all_hists:
            renormalize_for_bin_widths(h, bins)

    ## stack background histograms
    if h_bkgs:
        assert len(h_bkgs) == len(n_bkgs), '%s\n%s' % (h_bkgs, n_bkgs)

        h_bkgs.reverse()
        n_bkgs.reverse()

        hstack = HistStack()
        for h in h_bkgs:
            hstack.Add(h)
        hstack.title = 'stack sum'
        hstack.drawstyle='HIST'
        stuff['stack'] = hstack

        h_bkgs.reverse()
        n_bkgs.reverse()

#    ## convert data to TGraphAsymmErrors
#    g_data = None
#    if h_data:
#        if __use_poissonize:
#            g_data = poissonize.GetPoissonizedGraph(h_data)
#        else:
#            g_data = ROOT.TGraphAsymmErrors()
#            i_g = 0
#            nbins = h_data.GetNbinsX()
#            for i_bin in xrange(1, nbins+1): # skip underflow/overflow
#                c = h_data.GetBinContent(i_bin)
#                e = h_data.GetBinError(i_bin)
#                if c != 0.0:
#                    g_data.SetPoint(i_g, h_data.GetBinCenter(i_bin), c)
#                    g_ratio.SetPointError(i_g, 
#                            h_data.GetBinWidth(i_bin)/2.,
#                            h_data.GetBinWidth(i_bin)/2.,
#                            e,
#                            e)
#                i_g += 1


    ## build list of objects to draw
    objects = list()
    if h_bkgs:
        objects.append(stuff['stack'])
        objects.append(stuff['h_bkg_total'])
    if h_sigs:
        objects.extend(h_sigs)
    if h_data:
        objects.append(h_data)

    ## set xlimits and ylimits
    ypadding = 0.21
    logy_crop_value = 7e-3
    xmin, xmax, ymin, ymax = 0.0, 1.0, 0.0, 1.0
    if objects:
        xmin, xmax, ymin, ymax = get_limits(objects, logx=logx, logy=logy, ypadding=ypadding, logy_crop_value=logy_crop_value)
    if logy:
        ymin = 7e-3
    else:
        ymin = 0.0
    xlimits = (xmin, xmax)
    ylimits = (ymin, ymax)
    stuff['xlimits'] = xlimits
    stuff['ylimits'] = ylimits

    ## remove xtitle for do_ratio
    _xtitle = xtitle
    if h_data and h_bkg_total and kwargs.get('do_ratio'):
        _xtitle = ''

    ## make canvas
    canvas = Canvas(800, 600)
    stuff['canvas'] = canvas

    ## draw the objects
    if objects:
        canvas.cd()
        draw(objects, pad=canvas, xtitle=_xtitle, ytitle=ytitle, xlimits=xlimits, ylimits=ylimits)

    ## set log x/y, for some reason doesn't work before draw
    if logx or logy:
        if logx:
            canvas.SetLogx()
        if logy:
            canvas.SetLogy()
        canvas.Update()

    ## draw blind_line
    if has_blinded_data:
        if isinstance(blind, tuple):
            blind_list = list(blind)
        else:
            blind_list = [blind]
        blind_lines = list()
        for bl in blind_list:
            line_y1 = ymin
            line_y2 = ymax
            blind_line = ROOT.TLine(bl, line_y1, bl, line_y2)  
            blind_line.SetLineColor(ROOT.kGray+2)
            blind_line.SetLineStyle(7)
            blind_line.SetLineWidth(2)
            blind_line.Draw()
            blind_lines.append(blind_line)
        stuff['blind_lines'] = blind_lines
        canvas.Update()

    ## legend
    lefty = True
    if h_bkg_total:
        lefty = is_left_sided(h_bkg_total)
    elif h_data:
        lefty = is_left_sided(h_data)
    elif h_sigs:
        lefty = is_left_sided(h_sigs[0])

    if main_hists:
        header = '%.1f fb^{-1}, 13 TeV' % (lumi/1000.0)
        if lefty:
            legend = Legend(main_hists,
                            pad=canvas,
                            header=header,
                            textsize=16,
                            topmargin=0.03,
                            leftmargin=0.60,
                            rightmargin=0.02,
                            entrysep=0.01,
                            entryheight=0.04)
        else:
            legend = Legend(main_hists,
                            pad=canvas,
                            header=header,
                            textsize=16,
                            topmargin=0.03,
                            leftmargin=0.03,
                            rightmargin=0.59,
                            entrysep=0.01,
                            entryheight=0.04)
        legend.Draw()
        stuff['legend'] = legend

    if h_sigs:
#        header = 'ATLAS Internal'
        header = ''
        if lefty:
            legend2 = Legend(h_sigs,
                            pad=canvas,
                            header=header,
                            textsize=16,
                            topmargin=0.03,
                            leftmargin=0.37,
                            rightmargin=0.23,
                            entrysep=0.01,
                            entryheight=0.04)
        else:
            legend2 = Legend(h_sigs,
                            pad=canvas,
                            header=header,
                            textsize=16,
                            topmargin=0.03,
                            leftmargin=0.20,
                            rightmargin=0.40,
                            entrysep=0.01,
                            entryheight=0.04)
        legend2.Draw()
        stuff['legend2'] = legend2

    ## do_ratio
    if h_data and h_bkg_total and kwargs.get('do_ratio'):

        ## top canvas
        top_canvas = stuff.pop('canvas')
        stuff['top_canvas'] = top_canvas

        ## make SM/SM with error band: h_ratio_band
        i_sfratio = int(kwargs.get('sfratio', -1))
        if i_sfratio < 0: # ratio plot of Data/Model
            h_ratio_band = h_bkg_total.Clone()
            nbins = h_ratio_band.GetNbinsX()
            for i_bin in xrange(nbins+2):
                h_ratio_band.SetBinContent(i_bin, 1.0)
                c = h_bkg_total.GetBinContent(i_bin)
                e = h_bkg_total.GetBinError(i_bin) / c if c > 0.0 else 0.0
                h_ratio_band.SetBinError(i_bin, e)
            stuff['h_ratio_band'] = h_ratio_band
        else: # ratio plot of Scale Factor for ith background
            hi = h_bkgs[i_sfratio]
            h_ratio_band = hi.Clone()
            nbins = h_ratio_band.GetNbinsX()
            for i_bin in xrange(nbins+2):
                h_ratio_band.SetBinContent(i_bin, 1.0)
                c = hi.GetBinContent(i_bin)
                e = hi.GetBinError(i_bin) / c if c > 0.0 else 0.0
                h_ratio_band.SetBinError(i_bin, e)
            stuff['h_ratio_band'] = h_ratio_band

        ## make data/(SM) h_ratio
        if i_sfratio < 0:
            h_ratio = h_data.Clone()
            h_ratio.Divide(h_data, h_bkg_total, 1.0, 1.0)
            stuff['h_ratio'] = h_ratio
        else:
            ## SF1 = 1.0 + (data - MCtot) / MC1
            sfname = kwargs.get('sfname')
            sffile = kwargs.get('sffile')
            if not sfname:
                sfname = 'h_sf'
            hi = h_bkgs[i_sfratio]
            h_numer = h_data.Clone()
            h_numer.Add(h_bkg_total, -1.0)
            ## do the division
            h_ratio = h_data.Clone(sfname)
            h_ratio.Divide(h_numer, hi, 1.0, 1.0)
            ## add the 1.0
            nbins = h_ratio.GetNbinsX()
            for i_bin in xrange(nbins+2):
                c = h_ratio.GetBinContent(i_bin)
                h_ratio.SetBinContent(i_bin, c+1.0)
                h_ratio_band.SetBinContent(i_bin, c+1.0)
            ## ignore bins with no data for SF
            for i_bin in xrange(nbins+2):
                c = h_data.GetBinContent(i_bin)
                if c <= 0:
                    h_ratio.SetBinContent(i_bin, 0.0)
                    h_ratio.SetBinError(i_bin, 0.0)
                    h_ratio_band.SetBinError(i_bin, 0.0)
            stuff['h_ratio'] = h_ratio
            if sffile:
                f_out = ipyhep.file.write(h_ratio, sffile)
#                f_out.Close()

        ## convert ratio to a TGraphErrors so that Draw('E0')
        ## shows error bars for points off the pad
        g_ratio = ROOT.TGraphErrors()
        i_g = 0
        for i_bin in xrange(1, nbins+1): # skip underflow/overflow
            ratio_content = h_ratio.GetBinContent(i_bin)
            if ratio_content != 0.0:
                g_ratio.SetPoint(i_g, h_ratio.GetBinCenter(i_bin), ratio_content)
                g_ratio.SetPointError(i_g,
                        h_ratio.GetBinWidth(i_bin)/2.,
                        h_ratio.GetBinError(i_bin) )
                i_g += 1
            else:
                h_ratio.SetBinError(i_bin, 0.0)
        stuff['g_ratio'] = g_ratio

        ## style ratio
        h_ratio_band.title = 'bkg uncert.'
        if i_sfratio < 0:
            h_ratio_band.linecolor = ipyhep.style.yellow
        else:
            h_ratio_band.linecolor = ipyhep.style.light_gray
        h_ratio_band.linewidth = 0
        h_ratio_band.markerstyle = 0
        if i_sfratio < 0:
            h_ratio_band.fillcolor = ipyhep.style.yellow
        else:
            h_ratio_band.linecolor = ipyhep.style.light_gray
        h_ratio_band.fillstyle = ipyhep.style.fill_solid
        h_ratio_band.drawstyle = 'E2'
        h_ratio_band.legendstyle = 'F'
        h_ratio.title = 'ratio'
        h_ratio.linecolor = ipyhep.style.black
        h_ratio.linewidth = 2
        h_ratio.markercolor = ipyhep.style.black
        h_ratio.markerstyle = 20
        h_ratio.markersize = 1.2
        h_ratio.fillstyle = ipyhep.style.fill_hollow
        h_ratio.drawstyle = 'PE'
        h_ratio.legendstyle = 'LP'

        ## bottom canvas
        bottom_canvas = Canvas(800, 600)
        bottom_canvas.cd()
        stuff['bottom_canvas'] = bottom_canvas

        ## set ratio ylimits
        ratio_min = kwargs.get('ratio_min', -0.2)
        ratio_max = kwargs.get('ratio_max',  2.2)
        ratio_ylimits = (ratio_min, ratio_max)

        ## draw ratio band
        if i_sfratio < 0:
            _ytitle = 'Data / Model'
        else:
            hi = h_bkgs[i_sfratio]
            _ytitle = 'SF(%s)' % hi.title
        draw([h_ratio_band], pad=bottom_canvas, xtitle=xtitle, ytitle=_ytitle, xlimits=xlimits, ylimits=ratio_ylimits)

        ## set log x/y, for some reason doesn't work before draw?
        if logx:
            bottom_canvas.SetLogx()
            bottom_canvas.Update()

        ### make horiz lines in ratio plot every 0.5:
        line_ys = [y / 10.0 for y in range(10*int(round(ratio_min)), 10*int(round(ratio_max))+5, 5)]
        line_x1 = canvas.GetUxmin()
        line_x2 = canvas.GetUxmax()
        line_xwidth = abs(line_x2-line_x1)
        lines = []
        for line_y in line_ys:
            line = ROOT.TLine(line_x1+0.02*line_xwidth, line_y, line_x2-0.02*line_xwidth, line_y)  
            line.SetLineWidth(1)
            line.SetLineStyle(7)
            if line_y == 1.0:
                line.SetLineColor(ROOT.kGray+2)
            else:
                line.SetLineColor(ROOT.kGray+0)
            line.Draw()
            lines.append(line)
        stuff['lines'] = lines

        ## draw blind_line
        if has_blinded_data:
            if isinstance(blind, tuple):
                blind_list = list(blind)
            else:
                blind_list = [blind]
            blind_lines = list()
            for bl in blind_list:
                line_y1 = ymin
                line_y2 = ymax
                blind_line = ROOT.TLine(bl, line_y1, bl, line_y2)  
                blind_line.SetLineColor(ROOT.kGray+2)
                blind_line.SetLineStyle(7)
                blind_line.SetLineWidth(2)
                blind_line.Draw()
                blind_lines.append(blind_line)
            stuff['blind_lines2'] = blind_lines
            canvas.Update()

        ## draw ratio
        g_ratio.Draw('PE0')
#        h_ratio.GetYaxis().SetRangeUser(ratio_min, ratio_max)
#        h_ratio.Draw('PE,SAME')

        ## shared canvas
        shared_canvas = Canvas(800, 800)
        shared_plot = plot_shared_axis(
            top_canvas,
            bottom_canvas,
            canvas = shared_canvas,
            split=0.35,
            axissep=0.01)
        stuff['canvas'] = shared_canvas
        canvas = shared_canvas

    ## save figures
    save = kwargs.get('save')

    if save is None: # NOTE: save can be False to skip saving
        save = ['pdf', 'png']

    if save:
        ipyhep.file.save_figures(canvas, x, save)

    global results
    results = stuff
    return stuff


#______________________________________________________________________________
def calc_fake_factors(name, numers_data, numers_bkg, denoms_data, denoms_bkg, labels, xtitle, ytitle, save=None, ylimits=None):

    assert numers_data, numers_bkg

    if not isinstance(numers_data, list):
        numers_data = [numers_data]
    if not isinstance(numers_bkg, list):
        numers_bkg = [numers_bkg]
    if not isinstance(denoms_data, list):
        denoms_data = [denoms_data]
    if not isinstance(denoms_bkg, list):
        denoms_bkg = [denoms_bkg]

    if save is None:
        save = ['pdf', 'png']

    stuff = dict()
    objects = list()

    colors = [
        ipyhep.style.red,
        ipyhep.style.blue,
        ipyhep.style.orange,
        ipyhep.style.violet,
        ]

    i_ff = 0
    for numer_data, numer_bkg, denom_data, denom_bkg in zip(numers_data, numers_bkg, denoms_data, denoms_bkg):

        for h in (numer_data, numer_bkg, denom_data, denom_bkg):
            if h and not h.GetSumw2N():
                h.Sumw2()

        h_numer = numer_data.Clone()
        h_numer.SetDirectory(0)
        if numer_bkg:
            set_negative_bins_zero(numer_bkg)
            h_numer.Add(numer_bkg, -1.0)
            
        h_denom = denom_data.Clone()
        h_denom.SetDirectory(0)
        if denom_bkg:
            set_negative_bins_zero(denom_bkg)
            h_denom.Add(denom_bkg, -1.0)
            
        h_ff = h_numer.Clone()
        h_ff.SetDirectory(0)
        divide_option = ''
        h_ff.Divide(h_numer, h_denom, 1.0, 1.0, divide_option)

        h_ff.linecolor = colors[i_ff]
        h_ff.markercolor = colors[i_ff]
        h_ff.markerstyle = 20
        h_ff.markersize = 1.2
        h_ff.fillstyle = ipyhep.style.fill_hollow
        h_ff.drawstyle = 'PE'
        h_ff.legendstyle = 'LP'
        if labels:
            h_ff.title = '%s' % labels[i_ff]
        else:
            h_ff.title = '%i' % i_ff
        h_ff.name = 'h_ff_%s' % h_ff.title

        stuff['h_ff_%s' % h_ff.title] = h_ff
        objects.append(h_ff)
        i_ff += 1

    ## make canvas
    canvas = Canvas(800, 600)
    stuff['canvas'] = canvas

    xmin, xmax, ymin, ymax = get_limits(objects)
    xlimits = (xmin, xmax)
#    ylimits = (ymin, ymax)
    ylimits = ylimits or (0.0, 1.0)

    draw(objects, pad=canvas, xtitle=xtitle, ytitle=ytitle, xlimits=xlimits, ylimits=ylimits)

    ## legend
    header = ''
    legend = Legend(objects,
                    pad=canvas,
                    header=header,
                    textsize=16,
                    topmargin=0.02,
                    leftmargin=0.60,
                    rightmargin=0.02,
                    entrysep=0.01,
                    entryheight=0.04)

    legend.Draw()
    stuff['legend'] = legend

    ## save figures
    if save:
        ipyhep.file.save_figures(canvas, 'h_ff', save)

    ## save fake factors
    timestamp = time.strftime('%Y_%m_%d_%Hh%M')
    outfile = 'fakefactors-%s.root' % timestamp
    f_out = None
    for h_ff in objects:
        f_out = ipyhep.file.write(h_ff, outfile)
#    f_out.Close()

#    ipyhep.file.close_all_files()

    return stuff

#______________________________________________________________________________
def set_negative_bins_zero(h, warn=False):
    for bin in xrange(0, h.GetNbinsX()+2):
        if h.GetBinContent(bin) < 0:
            h.SetBinContent(bin, 0.0)
            h.SetBinError(bin, 0.0)  ## NOTE
            if warn:
                print 'WARNING (fakeestimate.py): setting negative bin to zero in'
                print '  %s' % h.GetName()
    return h


#______________________________________________________________________________
#def significance_scan_2d(x, bins=(16, 2000., 6000., 10, 100., 300.), xtitle='#font[12]{H}_{#font[132]{T}}  [GeV]', ytitle='#font[12]{E}_{#font[132]{T}}^{#font[132]{miss}}  [GeV]', selection='', weight='', bkgs=None, sig=None, lumi=None, datadrivensearchpath=None, bkgsearchpath=None, sigsearchpath=None, save=None):
def significance_scan_2d(x, bins=(16, 2000., 6000., 10, 100., 300.), xtitle='#font[12]{H}_{#font[132]{T}}  [GeV]', ytitle='#font[12]{E}_{#font[132]{T}}^{#font[132]{miss}}  [GeV]', selection='', weight=''):
    """
    The signficance is calculated using Cowan's analytic formula for discovery
    significance of a signal with background and uncertainty on the background.
    See: https://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf
    """

    global data
    global bkgs
    global sigs
    global treename
    global datasearchpath
    global datadrivensearchpath
    global bkgsearchpath
    global sigsearchpath
    global lumi

    assert x
    assert len(bins) == 6
    assert bkgs
    assert sigs
    assert lumi
    assert bkgsearchpath
    assert sigsearchpath

    save = ['pdf', 'png']

    ## only use first sig
    assert len(sigs) == 1
    sig = sigs[0]

    ## save stuff to bookkeep and return
    stuff = dict()
    stuff['x'] = x

    ## get background histograms
    h_bkgs = list()
    n_bkgs = list()
    if bkgs:
        for bkg in bkgs:
            if isinstance(bkg, list):
                h_subtotal = None
                for dsid in bkg:
                    assert isinstance(dsid, str)
                    h_bkg = None
                    if dsid.isdigit():
                        ## mc backgrounds
                        sp = bkgsearchpath % int(dsid)
                        newx = '%s::%s::%s' % (sp, treename, x)
                        h_bkg = ipyhep.tree.project(newx, bins=bins, xtitle=xtitle, ytitle=ytitle, selection=selection, weight=weight)
                    else:
                        ## data-driven backgrounds
                        assert dsid == 'fakes' or dsid == 'efakes'
                        sp = datadrivensearchpath
                        newx = '%s::%s::%s' % (sp, treename, x)
                        h_bkg = ipyhep.tree.project(newx, bins=bins, xtitle=xtitle, ytitle=ytitle, selection=selection, weight=weight)
                    if h_bkg:
                        if h_subtotal:
                            h_subtotal.Add(h_bkg)
                        else:
                            h_subtotal = h_bkg.Clone()
                if h_subtotal:
                    h_bkgs.append(h_subtotal)
                    dsid = bkg[0] ## lists of combined backgrounds use the first dsid
                    n_bkgs.append(dsid)
            else:
                dsid = bkg
                assert isinstance(dsid, str)
                h_bkg = None
                if dsid.isdigit():
                    ## mc backgrounds
                    sp = bkgsearchpath % int(dsid)
                    newx = '%s::%s::%s' % (sp, treename, x)
                    h_bkg = ipyhep.tree.project(newx, bins=bins, xtitle=xtitle, ytitle=ytitle, selection=selection, weight=weight)
                else:
                    ## data-driven backgrounds
                    assert dsid == 'fakes' or dsid == 'efakes'
                    sp = datadrivensearchpath
                    newx = '%s::%s::%s' % (sp, treename, x)
                    h_bkg = ipyhep.tree.project(newx, bins=bins, xtitle=xtitle, ytitle=ytitle, selection=selection, weight=weight)
                if h_bkg:
                    h_bkgs.append(h_bkg)
                    n_bkgs.append(dsid)
    assert h_bkgs

    ## get signal histograms
    h_sig = None
    if sig:
        dsid = sig
        sp = sigsearchpath % int(dsid)
        newx = '%s::%s::%s' % (sp, treename, x)
        h_sig = ipyhep.tree.project(newx, bins=bins, xtitle=xtitle, ytitle=ytitle, selection=selection, weight=weight)
    assert h_sig

    ## scale background histograms
    if h_bkgs:
        assert len(h_bkgs) == len(n_bkgs), '%s\n%s' % (h_bkgs, n_bkgs)
        for h, dsid in zip(h_bkgs, n_bkgs):
            sf = ipyhep.sampleops.get_sf(dsid)
            if dsid.isdigit():
                sf *= lumi/__ntuple_lumi
            h.Scale(sf)

    ## scale signal histogram
    if h_sig:
        sf = ipyhep.sampleops.get_sf(sig)
        sf *= lumi/__ntuple_lumi
        h_sig.Scale(sf)

    ## set systematic errors on backgrounds
    if h_bkgs:
        assert len(h_bkgs) == len(n_bkgs), '%s\n%s' % (h_bkgs, n_bkgs)

        # HACK: systematics updated 2017-06-08 for the SUSY diphoton analysis
        syst_fracs = dict()
        syst_fracs['407013'] = 0.50 # #gamma#gamma
        syst_fracs['361039'] = 0.50 # #gammaj+jj
        syst_fracs['fakes']  = 0.50 # fakes
        syst_fracs['efakes'] = 0.20 # efakes
        syst_fracs['301890'] = 0.20 # W#gamma
        syst_fracs['301899'] = 0.20 # Z#gamma
        syst_fracs['407022'] = 0.27 # W#gamma#gamma
        syst_fracs['407025'] = 0.45 # Z#gamma#gamma
        syst_fracs['407028'] = 0.45 # Z#gamma#gamma

        for h, dsid in zip(h_bkgs, n_bkgs):
            nbinsx = h.GetNbinsX()
            nbinsy = h.GetNbinsY()
            for i_x in xrange(nbinsx):
                for i_y in xrange(nbinsy):
                    c = h.GetBinContent(i_x, i_y)
                    e = h.GetBinError(i_x, i_y)
                    syst = syst_fracs.get(dsid, 0.0)
                    assert syst
                    h.SetBinError(i_x, i_y, math.sqrt(e*e + c*c*syst*syst))

    ## total background
    h_bkg_total = None
    if h_bkgs:
        for h_bkg in h_bkgs:
            if h_bkg_total:
                h_bkg_total.Add(h_bkg)
            else:
                h_bkg_total = h_bkg.Clone()
    assert h_bkg_total

    ## significance scan
    h_signif = Hist2D(*bins)
    h_signif.SetTitle(';%s;%s' % (xtitle, ytitle))
    nbinsx = h_signif.GetNbinsX()
    nbinsy = h_signif.GetNbinsY()
    for i_x in xrange(1, nbinsx+1):
        for i_y in xrange(1, nbinsy+1):
            sigma = ROOT.Double(0)
            s =  h_sig.Integral(i_x, nbinsx+1, i_y, nbinsy+1)
            b =  h_bkg_total.IntegralAndError(i_x, nbinsx+1, i_y, nbinsy+1, sigma)
            sigma = float(sigma)
#            if b < 0.10: ## HACK: force background to be at least 0.10 events with 100% uncert.
#                b = 0.10
#                sigma = 0.10
            if b < 0.20 and sigma < 0.20:
                sigma = 0.20
            za = _calculate_significance(s, b, sigma)
            xval = h_signif.GetXaxis().GetBinLowEdge(i_x)
            yval = h_signif.GetYaxis().GetBinLowEdge(i_y)
            h_signif.Fill(xval, yval, za)

    ## make canvas
    canvas = Canvas(800, 600)
    canvas.SetRightMargin(0.17)
    canvas.cd()
    stuff['canvas'] = canvas

    ## draw
    ROOT.gStyle.SetPaintTextFormat("3.1f")
    h_signif.Draw('COLZ')
    h_signif2 = h_signif.Clone()
    h_signif2.SetMarkerSize(1.2)
    h_signif2.SetMarkerColor(ipyhep.style.black)
    h_signif2.Draw('TEXT SAME')
    stuff['h_signif'] = h_signif
    stuff['h_signif2'] = h_signif2

    canvas.Update()

    ## save figures
    if save:
        ipyhep.file.save_figures(canvas, x, save)

    return stuff
        

## only the above functions are exposed
__all__ = ['results', 'significance_scan_2d', 'stack']


#------------------------------------------------------------------------------
# Private helper functions
#------------------------------------------------------------------------------

#______________________________________________________________________________
def _calculate_significance(s, b, sigma):
    bt = 0.0
    sigmat = 0.0
    if isinstance(b, list):
        assert len(b) == len(sigma)
        for bi, sigmai in zip(b, sigma):
            bt += bi
            sigmat += bi*sigmai*sigmai
        sigmat = sigmat / bt if bt else 0.0
        sigmat = math.sqrt(sigmat)
    else:
        bt = b
        sigmat = sigma
    l1 = 0.0
    l2 = 0.0
    sq = 0.0
    za = 0.0
    if s and bt:
        try:
#            za = math.sqrt( 2.0*( (s+bt)*math.log(( (s+bt)*(bt+sigmat*sigmat) )/( bt*bt + (s+bt)*sigmat*sigmat )) - ((bt*bt)/(sigmat*sigmat))*math.log(1.0+(sigmat*sigmat*s)/(bt*(bt+sigmat*sigmat))) ) )
            l1 = ( (s+bt)*(bt+sigmat*sigmat) )/( bt*bt + (s+bt)*sigmat*sigmat )
            l2 = 1.0+(sigmat*sigmat*s)/(bt*(bt+sigmat*sigmat))
            if l1 > 0.0 and l2 > 0.0:
                sq = 2.0*( (s+bt)*math.log(l1) - ((bt*bt)/(sigmat*sigmat))*math.log(l2) )
                if sq > 0.0:
                    za = math.sqrt( sq )
                else:
                    print 'WARNING: sq = %5g' % sq
            else:
                print 'WARNING: l1 = %5g, l2 = %5g' % (l1, l2)
        except:
            print s, bt, sigmat
            assert False
    return za


#_______________________________________________________________________________
def renormalize_for_bin_widths(h, bins):
    ## renormalize based on bin width
    ## TODO: support variable bins for Hist2D
    assert len(bins) > 2
    nbinedges = len(bins)
    nbins = nbinedges - 1
    binwidths = list()
    for i_bin in xrange(1, nbinedges):
        x1 = bins[i_bin-1]
        x2 = bins[i_bin]
        bw = abs(x2 - x1)
        binwidths.append(bw)
    assert len(binwidths) + 1 == len(bins)
    min_bw = min(binwidths)
    for i_bw, bw in enumerate(binwidths):
        i_bin = i_bw + 1
        binunits = float(bw)/min_bw
        assert binunits.is_integer()
        if binunits > 1:
            binweight = 1.0/binunits
            h.SetBinContent(i_bin, h.GetBinContent(i_bin)*binweight)
            h.SetBinError(i_bin, h.GetBinError(i_bin)*binweight)
    return h


#_______________________________________________________________________________
def is_left_sided(hist):
#    mean  = hist.GetMean() # gives nan sometimes?
    mode  = hist.GetXaxis().GetBinCenter( hist.GetMaximumBin() )
    x_min = hist.GetXaxis().GetXmin()
    x_max = hist.GetXaxis().GetXmax()
    mid   = (x_max + x_min)/2.0
    return mode < mid


#_______________________________________________________________________________
def get_stats(h):
    los = list()
    if isinstance(h, ROOT.TH1) and not isinstance(h, ROOT.TH2):
        nbins       = h.GetNbinsX()
        entries     = h.GetEntries()
        err = ROOT.Double(0)
        integral    = h.IntegralAndError(0, nbins+1, err)
        mean        = h.GetMean()
        rms         = h.GetRMS()
        under       = h.GetBinContent(0)
        over        = h.GetBinContent(nbins+1)
        stats = dict()
        name = h.GetTitle()
        stats['name']    = name
        stats['entries'] = '%i' % entries
        stats['int']     = ('%i' % round(integral)) if integral >= 1000. else ('%.3g' % integral)
        stats['err']     = '%.2g' % err
        stats['mean']    = '%.3g' % mean
        stats['rms']     = '%.3g' % rms
        stats['under']   = '%.3g' % under
        stats['over']    = '%.3g' % over
        los.append( stats )
    elif isinstance(h, ROOT.TH2):
        nbins_x     = h.GetNbinsX()
        nbins_y     = h.GetNbinsY()
        entries     = h.GetEntries()
        err = ROOT.Double(0)
        integral    = h.IntegralAndError(0, nbins_x+1, 0, nbins_y+1, err)
        mean_x      = h.GetMean(1)
        rms_x       = h.GetRMS(1)
        mean_y      = h.GetMean(2)
        rms_y       = h.GetRMS(2)
        stats = dict()
        name = h.GetTitle()
        stats['name']    = name
        stats['entries'] = '%i'   % entries
        stats['int']     = ('%i' % round(integral)) if integral >= 1000. else ('%.3g' % integral)
        stats['err']     = '%.2g' % err
        stats['mean_x']  = '%.3g' % mean_x
        stats['rms_x']   = '%.3g' % rms_x
        stats['mean_y']  = '%.3g' % mean_y
        stats['rms_y']   = '%.3g' % rms_y
        los.append( stats )
    elif isinstance(h, ROOT.TGraph) \
            or isinstance(h, ROOT.TGraphErrors) \
            or isinstance(h, ROOT.TGraphAsymmErrors):
        print 'WARNING: get_stats( %s ) not implemented.' % type(h)
    elif isinstance(h, ROOT.THStack):
#        stack_stats = get_stats( h.GetStack().Last() )
#        assert len(stack_stats) == 1, type(h.GetStack().Last())
#        stack_stats[0]['name'] = 'stack sum'
#        los.extend( stack_stats )
        stack_hists_stats = list()
        for hist in h.GetHists():
            stack_hists_stats.extend( get_stats(hist) )
        stack_hists_stats.reverse()
        los.extend(stack_hists_stats)
    else:
        print 'WARNING: get_stats( %s ) not implemented.' % type(h)
    return los


#__________________________________________________________________________
def convert_stats_to_table(los):
    ## HACK: need to come up with a way to determine which stats to expect,
    ## and how to organize the table(s)
    if los[0].has_key('rms_x'): # TH2
        top_row = ['name', 'int', 'err', 'entries', 'mean_x', 'rms_x', 'mean_y', 'rms_y']
    else:
        top_row = ['name', 'int', 'err', 'entries', 'mean', 'rms', 'under', 'over']
    tab = [top_row]
    for stats in los:
        row = []
        for x in top_row:
            row.append( stats.get(x, '') )
        tab.append(row)
    return tab


#__________________________________________________________________________
def convert_table_to_html(tab):
    html = ['    <table>\n']
    is_first = True
    for row in tab:
        html += ['        <tr>']
        for i_col, col in enumerate(row):
            row[i_col] = str(col)
        if is_first:
            for col in row:
                html += ['<th>%s</th>' % col]
            is_first = False
        else:
            for col in row:
                html += ['<td>%s</td>' % col]
        html += ['</tr>\n']
    html += ['    </table>\n']
    html = ''.join(html)
    return html


#______________________________________________________________________________
def plot_shared_axis(top_canvas, bottom_canvas, canvas, split=0.5, axissep=0.0):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    canvas.cd()
    top_pad = ROOT.TPad('top_pad', '',  0, split, 1, 1, 0, 0, 0)
    top_pad.Draw()

    bottom_pad = ROOT.TPad('bottom_pad', '',  0, 0, 1, split, 0, 0, 0)
    bottom_pad.Draw()

    top_pad.cd()
    top_canvas.DrawClonePad()
    bottom_pad.cd()
    bottom_canvas.DrawClonePad()

    top_pad.SetTopMargin(canvas.GetTopMargin())#*1.0/(1.0-split))
    top_pad.SetBottomMargin(axissep)
    top_pad.SetRightMargin(canvas.GetRightMargin())
    top_pad.SetLeftMargin(canvas.GetLeftMargin());
    top_pad.SetFillStyle(0) # transparent
    top_pad.SetBorderSize(0)

    bottom_pad.SetTopMargin(axissep)
    bottom_pad.SetBottomMargin(0.33) #(canvas.GetBottomMargin())#*1.0/split)
    bottom_pad.SetRightMargin(canvas.GetRightMargin())
    bottom_pad.SetLeftMargin(canvas.GetLeftMargin());
    bottom_pad.SetFillStyle(0) # transparent
    bottom_pad.SetBorderSize(0)

    pads = [top_pad, bottom_pad]
    for i_pad, pad in enumerate(pads):
        prims = [ p.GetName() for p in pad.GetListOfPrimitives() ]
        for name in prims:
            h = pad.GetPrimitive(name)
            if isinstance(h, ROOT.TH1) or isinstance(h, ROOT.THStack) or isinstance(h, ROOT.TGraph) or isinstance(h, ROOT.TGraphErrors) or isinstance(h, ROOT.TGraphAsymmErrors):
                if isinstance(h, ROOT.TGraph) or isinstance(h, ROOT.THStack) or isinstance(h, ROOT.TGraphErrors) or isinstance(h, ROOT.TGraphAsymmErrors):
                    h = h.GetHistogram()
                h.GetXaxis().SetLabelFont(43)
                h.GetXaxis().SetLabelSize(30)
                h.GetXaxis().SetTitleFont(43)
                h.GetXaxis().SetTitleSize(32)
                h.GetYaxis().SetLabelFont(43)
                h.GetYaxis().SetLabelSize(30)
                h.GetYaxis().SetTitleFont(43)
                h.GetYaxis().SetTitleSize(32)
                h.SetTitleOffset(3.2, 'X')
                h.SetTitleOffset(1.6, 'Y')
                if i_pad == 0:
                    h.SetLabelSize(0.0, 'X')
                    h.GetXaxis().SetTitle("")
                    h.GetXaxis().SetNdivisions(507)
                    h.GetYaxis().SetNdivisions(506)
                else:
                    h.GetXaxis().SetNdivisions(507)
                    h.GetYaxis().SetNdivisions(403)

    # crash in ROOT::delete_TPad ()
    # see: http://root.cern.ch/phpBB3/viewtopic.php?f=14&t=13392
    ROOT.SetOwnership(top_canvas, False)
    ROOT.SetOwnership(bottom_canvas, False)
    ROOT.SetOwnership(canvas, False)

    canvas.Update()

    return {'canvas':canvas, 'top_pad':top_pad, 'bottom_pad':bottom_pad, 'top_canvas':top_canvas, 'bottom_canvas':bottom_canvas}


