
import ROOT

#______________________________________________________________________________
def calc_integral_and_error(h, xmin=None, xmax=None):
    nbins = h.GetNbinsX()
    error = ROOT.Double(0)
    if xmin is None:
        bin1 = 0 
    else:
        bin1 = h.FindBin(xmin)
    if xmax is None:
        bin2 = nbins+1
    else:
        bin2 = h.FindBin(xmax)-1
    integral = h.IntegralAndError(bin1, bin2, error)
    error = float(error)
    return integral, error


