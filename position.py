class Position:
    
    def __init__(self, pos, seq_dict):
        self.pos = pos #amino acid posiiton on multiple sequence alignment
        self.seq_dict = seq_dict #dictionary from Biopython (key = species, value = amino acid sequence)
        self.aa_dict = self.create_aa_dict()
        self.count_dict = self.initial_count() #dictionary that keeps track of amino acids and presence
        self.control_aa = self.find_control_aa() #most common amino acid at position
        self.sample_size = self.get_sample_size() 
        
    def initial_count(self): #keeps count of amino acids
        count_dict = {}
        for aa in self.aa_dict.values():
            if aa not in count_dict:
                count_dict[aa] = 1
            else:
                count_dict[aa] += 1
        return count_dict

    def unique(self): #creates count of unique amino acids, based on how many amino acids represented in count_dict
        return len(self.count_dict)

    def total(self): #total number of amino acids at position, based on 
        return len(self.aa_dict)

    def find_control_aa(self): #finds the most common amino acid at position ("control")
        max_key = next(iter(self.count_dict))
        for key in self.count_dict:
            if self.count_dict[key] > self.count_dict[max_key]:
                max_key = key
        return max_key

    def control_count(self): #finds how many times the control amino acid is present at position
        return self.count_dict[self.control_aa]

    def proportion(self): #calculates the proportion of amino acids at position being the control
        total_aa = self.total()
        control_aa_count = self.control_count()
        return control_aa_count / total_aa

    def create_aa_dict(self): 
        aa_dict = {}
        for species, sequences in self.seq_dict.items():
            aa_dict[species] = sequences[self.pos]
        return aa_dict

    def get_sample_size(self):
        print(self.seq_dict)
        sample_size = len(self.seq_dict)
    


        
