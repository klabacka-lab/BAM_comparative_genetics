class Position:

    """
    Position is a class that looks at the amino acids from a sequence alignment across different speices at a particular position.

    Attributes:
    pos (int) = number that indicates sequence position on alignment that object will focus on
    seq (dict) = dictionary with species as key and respective amino acid sequence as value
    aa_dict (dict) = dictionary with species as key and amino acid at position as value
    control_aa (str) = most common amino acid found at position
    sample_size (int) = number of amino acids at position (excludes '-')

    """

    def __init__(self, pos, seq_dict):

        """
        Initializes position object

        Parameters:
        pos (int) = number that indicates sequence position on alignment that object will focus on
        seq(dict) = dictionary with species as key and respective amino acid sequence as value
        """

        self.pos = pos
        self.seq_dict = seq_dict
        self.aa_dict = self.create_aa_dict()
        self.count_dict = self.create_count_dict()
        self.control_aa = self.find_control_aa()
        self.sample_size = self.get_sample_size()
        self.proportion = self.get_proportion()
        self.total_aa = self.find_total_aa()



    def create_count_dict(self):

        """
        Creates a dictionary with amino acids as key and how many times it is represented as the value

        Returns:
        count_dict (dict) = dictionary of amino acids and how many times it is represented at the position

        """

        count_dict = {}
        for aa in self.aa_dict.values():
            if aa not in count_dict:
                count_dict[aa] = 1
            else:
                count_dict[aa] += 1
        return count_dict

    def count_unique_aa(self):

        """
        Returns the number of unique amino acids at the position designated by object

        Returns:
            unique_count (int) = number of unique amino acids at the position
        """

        unique_count = len(self.count_dict)
        return unique_count

    def find_total_aa(self):

        """
        Returns the total number of amino acids at the position designated by the object

        Returns:
        total_aa (int) = total number of amino acids at the position
        """

        total_aa = len(self.aa_dict)
        return total_aa

    def find_control_aa(self):

        """
        Finds the most common amino acid at the position designated by the object

        Returns:
        max_key (str) = Most common amino acid as determined by count found in count_dict
        """

        max_key = next(iter(self.count_dict))
        for key in self.count_dict:
            if self.count_dict[key] > self.count_dict[max_key]:
                max_key = key
        return max_key

    def count_control(self):

        """
        Returns the amount of times the control amino acid (most common amino acid) is represented at the position designated by the object

        Returns:
        control_count (int) = Number of times control amino acid is found at position in alignment
        """

        control_count = self.count_dict[self.control_aa]
        return control_count

    def get_proportion(self):

        """
        Calculate's the proportion of the most common amino acid to the total amount of amino acids found at the position designated by the object

        Returns:
            aa_proportion (float) = proportion of most common amino acid to total amino acids
        """

        total_aa = self.find_total_aa()
        control_aa_count = self.count_control()
        aa_proportion = control_aa_count / total_aa
        return aa_proportion

    def create_aa_dict(self):

        """
        Creates a dictionary with the species as a key and the amino acid at the particular position, designated by the object, of the species' amino acid sequence

        Returns:
            aa_dict (dict) = dictionary with species name as the key and the amino acid at position as value
        """

        aa_dict = {}
        for species, sequences in self.seq_dict.items():
            aa_dict[species] = sequences[self.pos]
        return aa_dict

    def get_sample_size(self):

        """
        Calculate's sample size at position, separating gaps "-" from the
        amount of amino acids

        Returns:
            sample_size (int) = total number of present amino acids
        """

        sample_size = sum(1 for value in self.aa_dict.values() if value != '-')
        return(sample_size)




