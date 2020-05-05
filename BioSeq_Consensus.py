HELP_DOC = """
BIOLOGICAL SEQUENCE CONSENSUS
(version 1.0)
by Angelo Chan

A library of functions useful for deriving consensus sequences for DNA or Amino
Acids.

Different thresholds can be specified for determining consensus, including
minimum count required, maximum error permitted, and minimum percentage 
required.

Inputs can be provided as a list of strings (sequences), a list of dictionaries 
(count of occurences of pseudo-occurences), or a standard matrix (nucleotides/
amino acid positions in alphabetical order).

The produced consensus sequence can either be strict (one output per position),
or standard notation. DNA consensus sequences also have a redundant nucleotide
option.
"""



# Configurations ###############################################################



# Imported Modules #############################################################



# Enums ########################################################################

class OUTPUT_MODE:
    SINGULAR = 1
    STANDARD = 2
    AMBIGUOUS = 3
    
class INPUT_MODE:
    SINGULAR = 1
    STANDARD = 2
    AMBIGUOUS = 3
    
class SEQUENCE_TYPE:
    DNA = 1
    AMINO = 2
    RNA = 3



# Strings ######################################################################

STR__invalid_seq_type = "\nERROR: Invalid sequence type specified."
STR__invalid_consesus_format = "\nERROR: Invalid consensus format specified."
STR__different_lengths = "\nERROR: Sequences are of differing lengths."

STR__no_ambiguous_AA_input = "\nERROR: Ambiguous input option has not been "\
"implemented for amino acids."
STR__no_ambiguous_AA_output = "\nERROR: Ambiguous output option has not been "\
"implemented for amino acids."

STR__invalid_seq_char = "\nERROR: Invalid sequence character: {c}"
STR__invalid_consensus_seq = "\nERROR: Invalid consensus sequence:\n\t{s}"

STR__invalud_cutoff = "\nERROR: Invalid cutoff specified."



# Lists ########################################################################

LIST__dna = ["DNA", "Dna", "dna", "D", "d"]
LIST__amino = [] # Populated later



LIST__n = ["A", "C", "G", "T"]
LIST__n_r = ["A", "C", "G", "T", "R", "Y", "S", "W", "M", "K", "N", "X"]
LIST__n_2 = ["R", "Y", "S", "W", "M", "K"]
LIST__n_3 = ["B", "D", "H", "V"]
LIST__n_4 = ["N", "X"]
LIST__a = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q",
        "R", "S", "T", "V", "W", "Y"]
LIST__a_20 = ["X"]
LIST__all_a3 = ["Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys",
        "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp",
        "Tyr"]
LIST__n_x = LIST__n + ["N"]
LIST__a_x = LIST__a + ["X"]
LIST__n_x_ = LIST__n + ["-"]
LIST__a_x_ = LIST__a + ["-"]
LIST__fillers = ["_", "-", " ", "."]
###########################
# A # Ala # Alanine       #
# C # Cys # Cysteine      #
# D # Asp # Aspartate     #
# E # Glu # Glutamate     #
# F # Phe # Phenylalanine #
# G # Gly # Glycine       #
# H # His # Histidine     #
# I # Ile # Isoleucine    #
# K # Lys # Lysine        #
# L # Leu # Leucine       #
# M # Met # Methionine    #
# N # Asn # Asparagine    #
# P # Pro # Proline       #
# Q # Gln # Glutamine     #
# R # Arg # Arginine      #
# S # Ser # Serine        #
# T # Thr # Threonine     #
# V # Val # Valine        #
# W # Trp # Tryptophan    #
# Y # Tyr # Tyrosine      #
###########################



# Basic matches
LIST__n_match__A = ["A"] # [A]denin
LIST__n_match__C = ["C"] # [C]ytosine
LIST__n_match__G = ["G"] # [G]uanine
LIST__n_match__T = ["T"] # [T]hymine
LIST__n_match__R = ["A", "G"] # Pu[r]ine
LIST__n_match__Y = ["C", "T"] # P[y]rimidine
LIST__n_match__S = ["C", "G"] # [S]trong
LIST__n_match__W = ["A", "T"] # [W]eak
LIST__n_match__M = ["A", "C"] # A[m]ino
LIST__n_match__K = ["G", "T"] # [K]eto
LIST__n_match__B = ["C", "G", "T"] # B+1 (Not A)
LIST__n_match__D = ["A", "G", "T"] # D+1 (Not C)
LIST__n_match__H = ["A", "C", "T"] # H+1 (Not G)
LIST__n_match__V = ["A", "C", "G"] # T+2 (Not T)
LIST__n_match__N = ["A", "C", "G", "T"]
LIST__n_match__X = LIST__n_match__N



# Dictionaries #################################################################

DICT__matches = {
    "A": LIST__n_match__A, "C": LIST__n_match__C, "G": LIST__n_match__G,
    "T": LIST__n_match__T, "R": LIST__n_match__R, "Y": LIST__n_match__Y,
    "S": LIST__n_match__S, "W": LIST__n_match__W, "M": LIST__n_match__M,
    "K": LIST__n_match__K, "B": LIST__n_match__B, "D": LIST__n_match__D,
    "H": LIST__n_match__H, "V": LIST__n_match__V, "N": LIST__n_match__N,
    "X": LIST__n_match__X}

DICT__matches_S = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "M": "AC", "K": "GT",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGT", "X": "ACGT"}

DICT__matches_inverse_redundant = {
    ("A", "C"): "M", ("A", "G"): "R", ("A", "T"): "W",
    ("C", "G"): "S", ("C", "T"): "Y", ("G", "T"): "K",
    ("A", "C", "G"): "V", ("A", "C", "T"): "H",
    ("A", "G", "T"): "D", ("C", "G", "T"): "B"}



TEMPLATE__counts_DNA = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
TEMPLATE__counts_AA = {"A": 0, "C": 0, "D": 0, "E": 0, "F": 0, "G": 0, "H": 0,
        "I": 0, "K": 0, "L": 0, "M": 0, "N": 0, "P": 0, "Q": 0, "R": 0, "S": 0,
        "T": 0, "V": 0, "W": 0, "Y": 0, "X": 0}



# Resolve Lists ################################################################

for s in ["Amino Acids", "Amino Acid", "Amino", "A", "Protein"]:
    LIST__amino.append(s.upper())
    LIST__amino.append(s.capitalize())
    LIST__amino.append(s.lower())



# Resolve Dictionaries #########################################################



# Functions ####################################################################

def Consensus_From_List(sequences, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD, input_mode=INPUT_MODE.SINGULAR):
    """
    """
    counts_dicts = Sequence_List_To_Counts_Dict(sequences, seq_type, input_mode)
    consensus = Consensus_From_Dicts(counts_dicts, cutoff, min_coverage,
            seq_type, output_mode)
    return consensus

def Consensus_From_Dicts(counts_dicts, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD):
    """
    """
    matrix = Counts_Dict_To_Standard_Matrix(counts_dicts, seq_type)
    consensus = Consensus_From_Standard_Matrix(matrix, cutoff, min_coverage,
            seq_type, output_mode)
    return consensus

def Consensus_From_Standard_Matrix(matrix, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD):
    """
    """
    # Setup
    sb = ""
    if seq_type == SEQUENCE_TYPE.DNA or seq_type in LIST__dna:
        range_ = range(len(LIST__n_x_))
        list_ = LIST__n_x_
    elif seq_type == SEQUENCE_TYPE.AMINO or seq_type in LIST__amino:
        if output_mode == OUTPUT_MODE.AMBIGUOUS:
            raise Exception(STR__no_ambiguous_AA_output)
        range_ = range(len(LIST__a_x_))
        list_ = LIST__a_x_
    else:
        raise Exception(STR__invalid_seq_type)
    # Main loop
    for sublist in matrix:
        # min_coverage
        coverage = sum(sublist)
        if coverage < min_coverage: sb += "-"
        else:
            best_count = 0
            candidates = []
            # For every possible nucleotide/residue
            for i in range_:
                char = list_[i]
                count = sublist[i]
                # Best possible
                if cutoff == 0:
                    if count == best_count: candidates.append(char)
                    elif count > best_count:
                        best_count = count
                        candidates = [char]
                # Complete match
                elif cutoff == 1:
                    if count == coverage: candidates.append(char)
                # Minimum count
                elif cutoff > 1:
                    if count >= cutoff: candidates.append(char)
                # Percentage
                elif cutoff > 0:
                    percentage = (float(count))/coverage
                    if percentage > cutoff: candidates.append(char)
                # Mismatches
                elif cutoff < 0:
                    if count >= cutoff + coverage: candidates.append(char)
                # None of the above
                else:
                    raise Exception(STR__invalud_cutoff )
            # Processing
            length = len(candidates)
            # No matches
            if len(candidates) == 0: sb += "-"
            # 1 match
            elif len(candidates) == 1: sb += candidates[0]
            # Universal match
            elif "-" in candidates: sb += "-"
            # Multiple matches
            else:
                if output_mode == OUTPUT_MODE.SINGULAR: sb += "-"
                elif output_mode == OUTPUT_MODE.STANDARD:
                    sb += ("[" + "".join(candidates) + "]")
                elif output_mode == OUTPUT_MODE.AMBIGUOUS:
                    key = tuple(candidates)
                    c = DICT__matches_inverse_redundant.get(key, "-")
                    sb += c
                else:
                    raise Exception(STR__invalid_consesus_format)
    # Return
    return sb

def Sequence_List_To_Counts_Dict(sequences, seq_type,
        input_mode=INPUT_MODE.SINGULAR, normalize=True, existing_dicts=[]):
    """
    """
    if seq_type == SEQUENCE_TYPE.DNA or seq_type in LIST__dna:
        if input_mode == INPUT_MODE.SINGULAR:
            result = Sequence_List_To_Counts_Dict__DNA_Singular(sequences,
                    normalize, existing_dicts)
        elif input_mode == INPUT_MODE.STANDARD:
            result = Sequence_List_To_Counts_Dict__DNA_Standard(sequences,
                    normalize, existing_dicts)
        elif input_mode == INPUT_MODE.AMBIGUOUS:
            result = Sequence_List_To_Counts_Dict__DNA_Ambiguous(sequences,
                    normalize, existing_dicts)
        else:
            raise Exception(STR__invalid_consesus_format)
    elif seq_type == SEQUENCE_TYPE.AMINO or seq_type in LIST__amino:
        if input_mode == INPUT_MODE.SINGULAR:
            result = Sequence_List_To_Counts_Dict__AA_Singular(sequences,
                    normalize, existing_dicts)
        elif input_mode == INPUT_MODE.STANDARD:
            result = Sequence_List_To_Counts_Dict__AA_Standard(sequences,
                    normalize, existing_dicts)
        elif input_mode == INPUT_MODE.AMBIGUOUS:
            raise Exception(STR__no_ambiguous_AA_input)
        else:
            raise Exception(STR__invalid_consensus_format)
    else:
        raise Exception(STR__invalid_seq_type)
    return result

def Counts_Dict_To_Standard_Matrix(counts, seq_type):
    """
    """
    # Determine matrix size based on sequence type
    if seq_type == SEQUENCE_TYPE.DNA or seq_type in LIST__dna:
        order = LIST__n_x
    elif seq_type == SEQUENCE_TYPE.AMINO or seq_type in LIST__amino:
        order = LIST__a_x
    else:
        raise Exception(STR__invalid_seq_type)
    
    # Prepare
    result = []
    
    # Populate using counts dictionary
    for dictionary in counts:
        temp = []
        for k in order:
            count = dictionary.get(k, 0)
            temp.append(count)
        result.append(temp)
    
    # Return
    return result

def Parse_Standard_Consensus_Sequence(string):
    """
    """
    result = []
    # Setup
    length = len(string)
    range_ = range(length)
    flag_brackets_curve = False
    flag_brackets_square = False
    sb_prev = ""
    sb_multi = ""
    sb_number = ""
    # Main loop
    for i in range_:
        c = string[i]
        # Number
        if c.isdigit():
            if (not flag_brackets_curve) or flag_brackets_square:
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            sb_number += c
        # Open brackets (curved)
        elif c == "(":
            if flag_brackets_curve or flag_brackets_square:
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            flag_brackets_curve = True
        # Close brackets (curved)
        elif c == ")":
            if (not flag_brackets_curve) or flag_brackets_square:
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            try:
                number = int(sb_number) - 1
            except:
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            result += [sb_prev] * number
            sb_number = ""
            flag_brackets_curve = False
        # Open brackets (square)
        elif c == "[":
            if flag_brackets_curve or flag_brackets_square:
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            flag_brackets_square = True
        # Close brackets (square)
        elif c == "]":
            if flag_brackets_curve or (not flag_brackets_square):
                raise Exception(STR__invalid_consensus_seq.format(s=string))
            result.append(sb_multi)
            sb_multi = ""
            flag_brackets_square = False
        # Other
        else:
            # In curve brackets
            if flag_brackets_curve and (not flag_brackets_square):
                sb_number += c
            # In square brackets
            elif (not flag_brackets_curve) and flag_brackets_square:
                sb_multi += c
            # Other
            else:
                result.append(c)
                sb_prev = c
    # Follow up
    if flag_brackets_curve or flag_brackets_square:
        raise Exception(STR__invalid_consensus_seq.format(s=string))
    # Return
    return result

def Parse_Ambiguous_Sequence(string):
    """
    """
    result = []
    # Main loop
    for n in string:
        matches = DICT__matches_S.get(n, "")
        if matches: result.append(matches)
        else: Exception(STR__invalid_consensus_seq.format(s=string))
    # Return
    return result 



# Sub-Functions ################################################################

def Sequence_List_To_Counts_Dict__DNA_Singular(sequences, normalize=True,
        existing_dicts=[]):
    """
    Subfunction of Sequence_List_To_Counts_Dict.
    """
    # Initialize
    length = len(sequences[0])
    range_ = range(length)
    if existing_dicts:
        result = []
        for d in existing_dicts: result.append(dict(d))
    else:
        result = []
        for i in range_:
            result.append(dict(TEMPLATE__counts_DNA))
    # Main loop
    for sequence in sequences:
        # Check length
        check = len(sequence)
        if check != length: raise Exception(STR__different_lengths)
        # Process
        for i in range_:
            c = sequence[i]
            if c in LIST__n_x:
                result[i][c] += 1
            elif c not in LIST__fillers:
                raise Exception(STR__invalid_seq_char.format(c=c))
    # Return
    return result

def Sequence_List_To_Counts_Dict__DNA_Standard(sequences, normalize=True,
        existing_dicts=[]):
    """
    Subfunction of Sequence_List_To_Counts_Dict.
    """
    # Initialize
    list_ = Parse_Standard_Consensus_Sequence(sequences[0])
    length_check = len(list_)
    range_ = range(length_check)
    if existing_dicts:
        result = []
        for d in existing_dicts: result.append(dict(d))
    else:
        result = []
        for i in range_:
            result.append(dict(TEMPLATE__counts_DNA))
    # Main loop
    for sequence in sequences:
        # Check length
        list_ = Parse_Standard_Consensus_Sequence(sequence)
        length = len(list_)
        if length != length_check: raise Exception(STR__different_lengths)
        # Process
        for i in range_:
            nucleotides = list_[i]
            score = 12/(len(nucleotides)) # 12 is the LCM of [1,2,3,4]
            for n in nucleotides: # Go through every N at that position
                # Nucleotides listed
                if n in LIST__n:
                    result[i][n] += score
                # "N" or "X"
                elif n in LIST__n_4:
                    for n2 in LIST__n: result[i][n2] += 3 # 12/4
                # Error
                elif n not in LIST__fillers:
                    raise Exception(STR__invalid_seq_char.format(c=n))
    # Normalize
    if normalize:
        for dictionary in result:
            for key in dictionary:
                dictionary[key] = dictionary[key]/232792560.0
    # Return
    return result

def Sequence_List_To_Counts_Dict__DNA_Ambiguous(sequences, normalize=True,
        existing_dicts=[]):
    """
    Subfunction of Sequence_List_To_Counts_Dict.
    """
    # Initialize
    list_ = Parse_Ambiguous_Sequence(sequences[0])
    length_check = len(list_)
    range_ = range(length_check)
    if existing_dicts:
        result = []
        for d in existing_dicts: result.append(dict(d))
    else:
        result = []
        for i in range_:
            result.append(dict(TEMPLATE__counts_DNA))
    # Main loop
    for sequence in sequences:
        # Check length
        list_ = Parse_Ambiguous_Sequence(sequence)
        length = len(list_)
        if length != length_check: raise Exception(STR__different_lengths)
        # Process
        for i in range_:
            nucleotides = list_[i]
            score = 12/(len(nucleotides)) # 12 is the LCM of [1,2,3,4]
            for n in nucleotides: # Go through every N at that position
                # Nucleotides listed
                if n in LIST__n:
                    result[i][n] += score
                # "N" or "X"
                elif nucleotides in LIST__n_4:
                    for n2 in LIST__n: result[i][n2] += 3 # 12/4
                # Error
                elif n not in LIST__fillers:
                    raise Exception(STR__invalid_seq_char.format(c=n))
    # Normalize
    if normalize:
        for dictionary in result:
            for key in dictionary:
                dictionary[key] = dictionary[key]/232792560.0
    # Return
    return result

def Sequence_List_To_Counts_Dict__AA_Singular(sequences, normalize=True,
        existing_dicts=[]):
    """
    Subfunction of Sequence_List_To_Counts_Dict.
    """
    # Initialize
    length = len(sequences[0])
    range_ = range(length)
    if existing_dicts:
        result = []
        for d in existing_dicts: result.append(dict(d))
    else:
        result = []
        for i in range_:
            result.append(dict(TEMPLATE__counts_AA))
    # Main loop
    for sequence in sequences:
        # Check length
        check = len(sequence)
        if check != length: raise Exception(STR__different_lengths)
        # Process
        for i in range_:
            c = sequence[i]
            if c in LIST__a_x:
                result[i][c] += 1
            elif c not in LIST__fillers:
                raise Exception(STR__invalid_seq_char.format(c=c))
    # Return
    return result

def Sequence_List_To_Counts_Dict__AA_Standard(sequences, normalize=True,
        existing_dicts=[]):
    """
    Subfunction of Sequence_List_To_Counts_Dict.
    """
    # Initialize
    list_ = Parse_Standard_Consensus_Sequence(sequences[0])
    length_check = len(list_)
    range_ = range(length_check)
    if existing_dicts:
        result = []
        for d in existing_dicts: result.append(dict(d))
    else:
        result = []
        for i in range_:
            result.append(dict(TEMPLATE__counts_AA))
    # Main loop
    for sequence in sequences:
        # Check length
        list_ = Parse_Standard_Consensus_Sequence(sequence)
        length = len(list_)
        if length != length_check: raise Exception(STR__different_lengths)
        # Process
        for i in range_:
            aminos = list_[i]
            score = 232792560/(len(aminos)) # 232792560 is the LCM of [1,...,20]
            for a in aminos: # Go through every AA at that position
                # Nucleotides listed
                if a in LIST__a:
                    result[i][a] += score
                # "N" or "X"
                elif aminos in LIST__a_20: # 11639628 is 232792560/20
                    for a2 in LIST__a: result[i][a2] += 11639628
                # Error
                elif a not in LIST__fillers:
                    raise Exception(STR__invalid_seq_char.format(c=a))
    # Normalize
    if normalize:
        for dictionary in result:
            for key in dictionary:
                dictionary[key] = dictionary[key]/232792560.0
    # Return
    return result



# Classes ######################################################################

class Consensus_Builder:
    """
    """
    def __init__(self, seq_type, input_mode=INPUT_MODE.SINGULAR):
        """
        """
        # Stored args
        self._type = seq_type
        self._mode = input_mode
        # Stored data
        self._counts_dicts = []
        self._normalized = []
        # Variables checking
        if not (self._type == SEQUENCE_TYPE.DNA or self._type in LIST__dna or
                self._type == SEQUENCE_TYPE.AMINO or self._type in LIST__amino):
            raise Exception()
    
    def Add_Seq(self, seq, input_mode=None):
        """
        """
        if not input_mode: input_mode = self._mode
        self._counts_dicts = Sequence_List_To_Counts_Dict([seq], self._type,
                input_mode, False, self._counts_dicts)
    
    def Find_Consensus(self, cutoff, min_coverage,
            output_mode=OUTPUT_MODE.STANDARD):
        """
        """
        self.Calculate_Normalized()
        consensus = Consensus_From_Dicts(self._counts_dicts, cutoff,
                min_coverage, self._type, output_mode)
        return consensus
    
    def Calculate_Normalized(self):
        """
        """
        # Determine factor
        if self._type == SEQUENCE_TYPE.DNA or self._type in LIST__dna:
            factor = 12.0
        else: factor = 232792560.0
        # Copy
        copy = []
        for d in self._counts_dicts: copy.append(dict(d))
        # Normalize
        for dictionary in copy:
            for key in dictionary:
                dictionary[key] = dictionary[key]/factor
        # Update
        self._normalized = copy



# Controlled Print Statements ##################################################



# Main Loop ####################################################################


