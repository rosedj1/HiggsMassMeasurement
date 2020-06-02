def Root_Hist_GetLastBinRightEdge(hist):
    """
    Returns the right-most bin edge of a ROOT.TH1 object. 
    """
    n_bins = hist.GetNbinsX()
    edge = hist.GetBinLowEdge(n_bins) + hist.GetBinWidth(n_bins)
    return edge