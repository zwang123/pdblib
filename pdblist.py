class PDBList(list):
    """
    store as list of dict, but can function as a dict of list
    """
    def __getitem__(self, key):
        #if type(key) is int:
        try:
            if abs(float(key) - int(key)) > 1e-6:
                raise ValueError
            return super().__getitem__(int(key))
        except ValueError:
            try:
                return [x[key] for x in self]
            except:
                return [[x[key] for x in y ] for y in self]
