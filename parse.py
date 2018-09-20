from collections import OrderedDict
from .pdblist import PDBList

def parse_pdb_line(l):
    """
    Parse a line in a pdb file and convert it to a OrderedDict
    """
    l = l.rstrip()
    return OrderedDict(
        [
            # TODO support more than 100000 atoms
            ("Record"    ,       l[ 0: 6].strip()),
            ("serial"    ,   int(l[ 6:11])       ),
            ("name"      ,       l[12:16].strip()),
            ("altLoc"    ,       l[16:17].strip()),
            ("resName"   ,       l[17:21].strip()),
            ("chainID"   ,       l[21:22].strip()),
            # TODO support more than 10000 residues
            ("resSeq"    ,   int(l[22:26])       ),
            ("iCode"     ,       l[26:27].strip()),
            ("x"         , float(l[30:38])       ),
            ("y"         , float(l[38:46])       ),
            ("z"         , float(l[46:54])       ),
            ("occupancy" , float(l[54:60])       ),
            ("tempFactor", float(l[60:66])       ),
            ("segment"   ,       l[72:76].strip()),
            ("element"   ,       l[76:78].strip()),
            ("charge"    ,       l[78:80].strip()),
        ]
    )

def parse_pdb_file(filename):
    """
    Given the filename, parse a pdb file and convert it to a list of OrderedDict
    """
    with open(filename) as f:
        rtn = PDBList()
        for l in f:
            if l[:6] in ["ATOM  ", "HETATM"]:
                rtn.append(parse_pdb_line(l))
    return rtn

def write_pdb_line(entry):
    """
    Convert a OrderedDict to a pdb line
    """
    name = entry["name"]
    if len(name) < 4:
        name = ' ' + name
    mod_entry = entry.copy()
    mod_entry["name"] = name
    return '{: <6s}{: >5d} {: <4s}{:1s}{: <4s}{:1s}{: >4d}{:1s}   {: >8.3f}{: >8.3f}{: >8.3f}{: >6.2f}{: >6.2f}      {: <4s}{: >2s}{: >2s}'.format(*mod_entry.values()).rstrip() + '\n'

def write_pdb_lines(pdbdata):
    """
    Convert a list of OrderedDict to a pdb file string
    """
    return ''.join([write_pdb_line(entry) for entry in pdbdata])

def write_pdb_file(pdbdata, filename):
    """
    Convert a list of OrderedDict to a pdb file
    """
    with open(filename, "w") as f:
        f.write(write_pdb_lines(pdbdata))

if __name__ == "__main__":
    #parse_pdb_line("ATOM   1914  SOD SOD S1127       0.016  -3.389  -0.040  1.00 58.57      S   NA")
    source = "/home/zhiwang/my_proj/charmm-gui_e86p_A51_drude/step2_drude.pdb"
    source = "/home/zhiwang/my_proj/charmm-gui_e86p_A51/step5_assembly.namd.pdb"
    #source = "orig.pdb"
    pdb = (parse_pdb_file(source))
    #print(pdb[6].values())
    with open("dd.pdb", "w") as f:
        f.write(write_pdb_file(pdb))
    with open(source) as f:
        with open('ref.pdb', 'w') as fout:
            for l in f:
                fout.write(l.rstrip() + '\n')
