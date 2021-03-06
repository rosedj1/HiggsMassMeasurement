def print_header_message(msg, pad_char="-", n_center_pad_chars=3):
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

def info():
    """Do Suzanne's fancy colored printing here."""
    pass