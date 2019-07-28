from TaurenMD import log

def is_selection_special_string(selection):
    """
    str -> bool
    
    Confirms if a string is a special selection string.
    
    Special selection strings have the following formats:
    - :INT
    - INT:
    - INT:INT
    - INT::INT
    - INT:INT:INT
    """
    try:
        sels = selection.split(":")
    except AttributeError:
        return False
    
    if len(sels) > 3:
        return False
    else:
        return all((n.isdigit() or n=='') for n in sels)


def is_selection_string_list(selection):
    
    try:
        valid = [
            bool(selection),
            isinstance(selection, str),
            all(n.isdigit() for n in selection.split())
            ]
    except AttributeError:  # if selection is None
        return False
    else:
        return all(valid)


def is_selection_list(selection):
    try:
        valid = [
            bool(selection),
            isinstance(selection, list),
            all(isinstance(n, int) for n in selection),
            ]
    except TypeError:  # if selection is None
        return False
    else:
        return all(valid)


def get_slice_from_special_string(selection):
    
    ls = selection.split(":")
    
    slices = [None, None, None]
    
    def add_to_slice(n):
        try:
            return int(n)
        except ValueError:
            return None
    
    for i, n in enumerate(ls):
        slices[i] = add_to_slice(n)
    
    return slice(*slices)
