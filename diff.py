from pdblib import parse_pdb_file, parse_pdb_line
from .pdblist import PDBList
#from sys import stderr

def equal(lhs, rhs, tol=1e-6):
    """
    An equal function considering floating points
    """
    try:
        return abs(lhs - rhs) < tol
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

def entry_compare(lhs, rhs, cmp_keys = ["x", "y", "z"]):
    """
    Return True if lhs and rhs mismatch for any key in cmp_keys
    """
    return not entry_match(lhs, rhs, cmp_keys)

def diff_pdb(query, subject, match=entry_match, compare=entry_compare):
    """
    calculate the difference between pdb files

    query   : in, short pdb filename
    subject : in, long pdb filename
    match   : in, a callable to determine whether two entries match,
              match(query_entry, subject_entry) = True if they match
    compare : in, a callable to determine whether two matching entries are 
              inequivalent, return True if not identical

    matchlist, the entries in subject file that match the query file
    diff1list, diff2lis, the entries where two files are different
    notfoundlist, the entries in the query file not found in the subject file
    """
    # pdb1 < pdb2
    pdb1 = parse_pdb_file(query)
    pdb2 = parse_pdb_file(subject)

    #lpdb2 = len(pdb2)
    matchlist = PDBList()
    diff1list = PDBList()
    diff2list = PDBList()
    notfoundlist = PDBList()
    #idx2 = 0
    for entry in pdb1:
        found = False
        #for idx2 in range(idx2, lpdb2):
        for idx2 in range(len(pdb2)):
            if match(entry, pdb2[idx2]):
                matchlist.append(pdb2[idx2])
                if compare(entry, pdb2[idx2]):
                    diff1list.append(entry)
                    diff2list.append(pdb2[idx2])
                del pdb2[idx2]
                found = True
                break
        #if idx2 == lpdb2:
        if not found:
            #print("Warning, entry not found:", entry, file=sys.stderr)
            #matchlist.append(parse_pdb_line(""))
            notfoundlist.append(entry)
            #return False, [], []

    return matchlist, diff1list, diff2list, notfoundlist
