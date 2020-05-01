from itertools import chain

def weave_lists(list1, list2):
    """
    Weave the elements of two lists together into a single list, in a predictable way.
    The initial lists can be of different lengths, in which case the remaining elements 
    of the longer list are append to the end of the merged list.

    Example:
    list1 = [1, 2, 3, 4]
    list2 = [5, 6]

    merged_list = [1, 5, 2, 6, 3, 4]

    """
    merged_list = list(chain.from_iterable(zip(list1, list2)))

    len1 = len(list1)
    len2 = len(list2)
    if (len1 != len2):
        # Missing remaining elements from longer list.
        min_len = min(len1, len2)

        if (min_len == len1):
            # list1 is the short boy.
            merged_list.extend(list2[min_len:])
        else:
            merged_list.extend(list1[min_len:])

    return merged_list 