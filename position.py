class Position:
    
    def __init__(self, pos, seq_dict):
        self.pos = pos
        self.seq_dict = seq_dict
        self.count_dict = self.initial_count()
        self.control_aa = self.find_control_aa()
        self.aa_dict = self.create_aa_dict()
        
    def initial_count(self):
        count_dict = {}
        for aa in self.aa_dict.values():
            if aa not in count_dict:
                count_dict[aa] = 1
            else:
                count_dict[aa] += 1
        return count_dict

    def unique(self):
        return len(self.count_dict)

    def total(self):
        return len(self.aa_dict)

    def find_control_aa(self):
        max_key = next(iter(self.count_dict))
        for key in self.count_dict:
            if self.count_dict[key] > self.count_dict[max_key]:
                max_key = key
        return max_key

    def control_count(self):
        return self.count_dict[self.control_aa]

    def proportion(self):
        total_aa = self.total()
        control_aa_count = self.control_count()
        return control_aa_count / total_aa

    def create_aa_dict(self):
        aa_dict = {}
        for species, sequences in self.seq_dict.items():
            aa_dict[species] = sequences[self.pos]
        return aa_dict
    


        
