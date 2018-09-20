from pdblib import parse_pdb_file, parse_pdb_line
from .pdblist import PDBList
from math import isclose

def equal(lhs, rhs):
    """
    An equal function considering floating points
    """
    try:
        return isclose(lhs, rhs)
    except:
        return lhs == rhs

def entry_match(lhs, rhs, matching_keys = 
                ["name", "resName", "chainID", "resSeq", "segment"]):
    """
    Return True if lhs and rhs match for all keys in matching_keys
    """
    for key in matching_keys:
        if not equal(lhs[key], rhs[key]):
            return False
    return True

def entry_compare(lhs, rhs, cmp_keys = ["x", "y", "z", "element"]):
    """
    Return True if lhs and rhs mismatch for any key in cmp_keys
    """
    return not entry_match(lhs, rhs, cmp_keys)

def diff_pdb(query, subject, match=entry_match, compare=entry_compare,
             filter=None, exclude=None):
    """
    calculate the difference between pdb files

    query   : in, short pdb filename
    subject : in, long pdb filename
    match   : in, a callable to determine whether two entries match,
              match(query_entry, subject_entry) = True if they match
    compare : in, a callable to determine whether two matching entries are 
              inequivalent, return True if not identical

    filter  : a function that changes a query item, return an entry, comes first
    exclude : a function that removes a query item, return True if excluding

    matchlist, the entries in subject file that match the query file
    diff1list, diff2lis, the entries where two files are different
    notfoundlist, the entries in the query file not found in the subject file
    """
    # pdb1 < pdb2
    pdb1 = parse_pdb_file(query)
    pdb2 = parse_pdb_file(subject)

    if filter is not None:
        pdb1 = PDBList([filter(entry) for entry in pdb1])

    matchlist = {
            "query"   : PDBList(),
            "subject" : PDBList(),
                }
    difflist = {
            "query"   : PDBList(),
            "subject" : PDBList(),
               }
    notfoundlist = PDBList()

    for entry in pdb1:
        if exclude is not None and exclude(entry):
            continue
        found = False
        for idx2 in range(len(pdb2)):
            if match(entry, pdb2[idx2]):
                matchlist["query"].append(entry)
                matchlist["subject"].append(pdb2[idx2])
                if compare(entry, pdb2[idx2]):
                    difflist["query"].append(entry)
                    difflist["subject"].append(pdb2[idx2])
                del pdb2[idx2]
                found = True
                break
        if not found:
            notfoundlist.append(entry)

    return matchlist, difflist, notfoundlist
