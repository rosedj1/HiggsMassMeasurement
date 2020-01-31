charge_dict = {
    '+':'-13', 
    '-':'13',
}

def combine_cut_list(cut_list):
    total_cuts = ''
    for cut in cut_list:
        total_cuts += cut + ' && '
    cut_str = total_cuts.rstrip(' && ')
    return cut_str

def calc_num_bins(bin_min, bin_max, bin_width):
    return int(round( (bin_max-bin_min)/bin_width ))

def make_str_title_friendly(cut_str, keep_whitespace=False):
    title = cut_str.replace('&&','_')
    title = title.replace('<','_lt_')
    title = title.replace('>','_gt_')
    if not (keep_whitespace): 
        title = title.replace(' ','')
    title = title.replace('.','p')
    title = title.replace('(','')
    title = title.replace(')','')
    title = title.replace('-','neg')
    title = title.replace('+','pos')
    title = title.replace('||','or')
    title = title.replace('-13','muPlus')
    title = title.replace('+13','muMinus')
    return title

def calc_ymin_for_legend(n_graphs, text_height=0.042):
    """
    Allows for dynamic expansion of bottom axis of TLegend. 
    There is probably an in-built function, but I want to write my own!
    
    Parameters
    ----------
    n_obj : int
        Number of objects (like graphs, histos) that will show up in your legend.
    text_height : float between 0 and 1
        Essentially sets the font. Was 0.033333 but you'll have to play with it. 
    """
    delta_y = text_height*n_graphs
    return 0.9 - delta_y # When you do TCanvas(), it puts the top axis at y = 0.9.

def print_header_message(msg):
    n = len(msg)
    octothorpes = (n+12)*'#'  # I'm probably the first person ever to name a variable 'octothorpes'.
    buff = 5*'#'
    print octothorpes 
    print buff, msg, buff
    print octothorpes
