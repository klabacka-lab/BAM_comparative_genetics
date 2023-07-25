class Position:
    
    def __init__(self, pos, aa_dict):
        self.pos = pos
        self.aa_dict = aa_dict
        self.count_dict = {}
        self.control_aa = ""
        
    def initial_count(self):
        for aa in self.aa_dict.values():
            if aa not in self.count_dict:
                self.count_dict[aa] = 1
            else:
                self.count_dict[aa] += 1

    def unique(self):
        return len(self.count_dict)

    def total(self):
        return len(self.aa_dict)

    def control(self):
        max_key = next(iter(self.count_dict))
        for key in self.count_dict:
            if self.count_dict[key] > self.count_dict[max_key]:
                max_key = key
        self.control_aa = max_key

    def control_count(self):
        return self.count_dict(self.control_aa)

    def proportion(self):
        total_aa = self.total()
        control_aa_count = self.control_count()
        return control_aa_count / total_aa

    def dict_create(self, filename):
        with open(filename, "r") as fasta_file:
            read_file = fasta_file.read()
            lines = read_file.split(">")
            for line in lines:
                temp = line.split("\n")
                aa_data = temp[1]
                self.aa_dict[temp[0]] = aa_data[self.pos]
    


        
