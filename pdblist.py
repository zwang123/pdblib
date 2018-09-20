class PDBList(list):
    """
    store as list of dict, but can function as a dict of list
    """
    def __getitem__(self, key):
        if type(key) is str:
            #if key is str
            try:
                return [x[key] for x in self]
            except:
                return [[x[key] for x in y ] for y in self]
        else:
            try:
                iter(key)
            except TypeError:
                #if key is not iterable or str
                return super().__getitem__(key)
        #if key is iterable but not str
        return [[x[k] for k in key] for x in self]

    def find_serial(self, serial):
        for i, x in enumerate(self):
            if x["serial"] == serial:
                return i
