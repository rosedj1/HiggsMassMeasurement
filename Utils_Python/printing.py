import json

def announce(msg, pad_char="-", n_center_pad_chars=3):
    """Print a message to stdout with custom padding characters.
    
    Example:
    msg = "Neat!"
    output:
    -------------
    --- Neat! ---
    -------------
    """
    banner = pad_char * (len(msg) + n_center_pad_chars*2 + 2)
    pad_char_short = pad_char * n_center_pad_chars
    middle = f"{pad_char_short} {msg} {pad_char_short}"
    print(banner)
    print(middle)
    print(banner)

def print_size(obj, name):
    """Print info about size of `obj` with `name`."""
    print(f"{name} has allocated memory: {obj.__sizeof__()}")
    
def info():
    """Do Suzanne's fancy colored printing here."""
    pass

def print_skipevent_msg(reason, evt_num, run=None, lumi=None, event=None):
    """Print info on why this event was skipped.
    
    Args:
        reason (str): Your explanation for why evt was skipped.
        evt_num (int): The entry in TTree (row).
        run, lumi, event
    """
    print(f"Skipping event index {evt_num}: {reason}")
    evtid = f"{run}:{lumi}:{event}".replace("None", "")
    if evtid != "::":
        # Provided at least one of run, lumi, event.
        print(f"  with Run:Lumi:Event = {evtid}")

def print_periodic_evtnum(evt_num, n_tot, print_every=500000):
    """Print event info: 'Event `print_every`/`n_tot`.' """
    if (evt_num % print_every) == 0:
        print(f"Processing event {evt_num}/{n_tot}.")

def pretty_print_dict(d):
    """Pretty print a dictionary."""
    print(json.dumps(d, indent=4))