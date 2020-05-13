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
    
EXAMPLE USAGE:

    >>> cb = Consensus_Builder("DNA", 1)
    >>> cb.Add_Seq("ATGC")
    >>> cb.Add_Seq("ATAT")
    >>> cb.Add_Seq("A-A-")
    >>> consensus_majority = cb.Find_Consensus(0.5, 2, 2)
    >>> consensus_absolute = cb.Find_Consensus(1, 2, 2)
    >>> consensus_majority_coverage3 = cb.Find_Consensus(0.5, 3, 2)

EXAMPLE USAGE:

    >>> consensus_1 = Consensus_From_List(["ATGC", "ATAT", "A-GC"], 0.5, 1,
            "DNA", 1, 1)
    >>> consensus_2 = Consensus_From_List(["ATGC", "ATAT", "A-GC"], 0.5, 3,
            "DNA", 1, 1)
    >>> consensus_3 = Consensus_From_List(["ATGC", "ATAT", "A-GC"], 1, 1,
            "DNA", 1, 1)
    >>> consensus_4 = Consensus_From_List(["ATGC", "ATAT", "A-CC"], 0, 1,
            "DNA", 1, 1)
    >>> consensus_5 = Consensus_From_List(["ATGC", "ATAT", "A-CC"], 0, 1,
            "DNA", 2, 1)
    >>> consensus_6 = Consensus_From_List(["ATGC", "ATAT", "[AT]-GC"], 2, 1,
            "DNA", 1, 2)
    >>> consensus_7 = Consensus_From_List(["ATGC", "ATAT", "[AT]-GC"], 2, 3,
            "DNA", 1, 2)
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
LIST__n_r = ["A", "C", "G", "T", "R", "Y", "S", "W", "M", "K",
        "B", "D", "H", "V", "N", "X"]
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

for s in ["Amino Acids", "Amino Acid", "Amino", "A", "Protein", "P"]:
    LIST__amino.append(s.upper())
    LIST__amino.append(s.capitalize())
    LIST__amino.append(s.lower())



# Resolve Dictionaries #########################################################



# Functions ####################################################################

def Consensus_From_List(sequences, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD, input_mode=INPUT_MODE.SINGULAR):
    """
    Derive a consensus sequence from a list of residue count dictionaries.
    
    @sequences
            (list<str>)
            A list of biological sequences. All sequences need to be the same
            length. Unknown residues ("N" and "X" for DNA, "X" for amino acids)
            are accepted, as are gaps. ("_", "-", " ", ".") For DNA sequences,
            ambiguous nucleotides (RYSWMKBDHVNX) are also accepted.
            All sequences must be the same length. If an existing list of
            dictionaries is supplied, all sequences must also be the same length
            as said list.
        
    @cutoff
            (int/float)
            The cutoff used to determine the qualification for a residue to be
            regarded as the consensus residue. Different cutoff methods will be
            used depending on the number specified:
                
                MISMATCH LIMIT
                (C < 0)
                    Only regard a residue as a/the consensus residue if there
                    are |@cutoff| or fewer mismatches across all counts.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                BEST AVAILABLE
                (C = 0)
                    Regard the residue(s) with the highest number of counts as
                    the consensus residue.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                PERCENTAGE
                (0 < C < 1)
                    Regard the residue(s) which account for a @cutoff fraction
                    or more of the all counts at that position. If the singular
                    output mode is specified, residues with more than one such
                    consensus residue will be regarded as a blank. (No consensus
                    residue at that position.)
                
                PERFECT MATCH
                (C = 1)
                    Only regard a residue as the consensus residue if all
                    residues at that position are that residue.
                
                MINIMUM COUNT
                (C > 1)
                    Only regard a residue as a/the consensus residue if there
                    are @cutoff or more counts of that residue at that position.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.
    
    @min_coverage
            (int)
            The minimum coverage required at a position in order for that
            position to have a consensus residue.
    
    @seq_type
            (int/str)
            An integer or string denoting the type of sequences that will be
            inputted.
            Accepted values for specifying DNA:
                SEQUENCE_TYPE.DNA
                1
                Any capitalization of "DNA" or it's abbreviations.
            Accepted values for specifying Amino Acid residues:
                SEQUENCE_TYPE.AMINO
                2
                Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                "Protein", or their abbreviations.
    
    @output_mode
            (int)
            An integer which specifies the format of the output consensus
            sequence:
                1:  SINGULAR
                    Only one residue at each position.
                2:  STANDARD
                    Accepts multiple residues possible for each position.
                    (Ex. "A[AT]GC")
                3:  AMBIGUOUS
                    (DNA only)
                    Only one character at each position, but that character can
                    represent multiple nucleotides.)
    
    @input_mode
            (int)
            An integer which specifies the format of the input sequences:
                1:  SINGULAR
                    Only one residue at each position.
                2:  STANDARD
                    Accepts multiple residues possible for each position.
                    (Ex. "A[AT]GC") Also acceptes repeats of the same residue
                    being denoted in short hand. (Ex. "AT(3)GC(2)" = "ATTTGCC")
                3:  AMBIGUOUS
                    (DNA only)
                    Only one character at each position, but that character
                    can represent multiple nucleotides.
    
    Consensus_From_Dicts(list<str>, int/float, int, int/str, int, int) -> str
    """
    counts_dicts = Sequence_List_To_Counts_Dict(sequences, seq_type, input_mode)
    consensus = Consensus_From_Dicts(counts_dicts, cutoff, min_coverage,
            seq_type, output_mode)
    return consensus

def Consensus_From_Dicts(counts_dicts, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD):
    """
    Derive a consensus sequence from a list of residue count dictionaries.

    @counts_dicts
            (list<dict<str:int/float>>)
            A list of dictionaries containing the residue counts at the
            corresponding positions.
        
    @cutoff
            (int/float)
            The cutoff used to determine the qualification for a residue to be
            regarded as the consensus residue. Different cutoff methods will be
            used depending on the number specified:
                
                MISMATCH LIMIT
                (C < 0)
                    Only regard a residue as a/the consensus residue if there
                    are |@cutoff| or fewer mismatches across all counts.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                BEST AVAILABLE
                (C = 0)
                    Regard the residue(s) with the highest number of counts as
                    the consensus residue.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                PERCENTAGE
                (0 < C < 1)
                    Regard the residue(s) which account for a @cutoff fraction
                    or more of the all counts at that position. If the singular
                    output mode is specified, residues with more than one such
                    consensus residue will be regarded as a blank. (No consensus
                    residue at that position.)
                
                PERFECT MATCH
                (C = 1)
                    Only regard a residue as the consensus residue if all
                    residues at that position are that residue.
                
                MINIMUM COUNT
                (C > 1)
                    Only regard a residue as a/the consensus residue if there
                    are @cutoff or more counts of that residue at that position.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.
    
    @min_coverage
            (int)
            The minimum coverage required at a position in order for that
            position to have a consensus residue.
    
    @seq_type
            (int/str)
            An integer or string denoting the type of sequences that will be
            inputted.
            Accepted values for specifying DNA:
                SEQUENCE_TYPE.DNA
                1
                Any capitalization of "DNA" or it's abbreviations.
            Accepted values for specifying Amino Acid residues:
                SEQUENCE_TYPE.AMINO
                2
                Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                "Protein", or their abbreviations.
    
    @output_mode
            (int)
            An integer which specifies the format of the output consensus
            sequence:
                1:  SINGULAR
                    Only one residue at each position.
                2:  STANDARD
                    Accepts multiple residues possible for each position.
                    (Ex. "A[AT]GC")
                3:  AMBIGUOUS
                    (DNA only)
                    Only one character at each position, but that character can
                    represent multiple nucleotides.)
    
    Consensus_From_Dicts(list<dict<str:int/float>>, int/float, int, int/str,
            int) -> str
    """
    matrix = Counts_Dict_To_Standard_Matrix(counts_dicts, seq_type)
    consensus = Consensus_From_Standard_Matrix(matrix, cutoff, min_coverage,
            seq_type, output_mode)
    return consensus

def Consensus_From_Standard_Matrix(matrix, cutoff, min_coverage, seq_type,
        output_mode=OUTPUT_MODE.STANDARD):
    """
    Derive a consensus sequence from a residue counts matrix.
    
    @matrix
            (list<list<int>>)
            A list of residue count lists. Each sublist contains the residue
            counts for the corresponding position. The integers in the sublist
            correspond to the number of ACGT or ACDEFGHIKLMNPQRSTVWY residues at
            that position. (ACGT for DNA, ACDEFGHIKLMNPQRSTVWY for amino acids)
            The last integer in each sublist indicates the number of residues
            which were ambiguous.
        
    @cutoff
            (int/float)
            The cutoff used to determine the qualification for a residue to be
            regarded as the consensus residue. Different cutoff methods will be
            used depending on the number specified:
                
                MISMATCH LIMIT
                (C < 0)
                    Only regard a residue as a/the consensus residue if there
                    are |@cutoff| or fewer mismatches across all counts.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                BEST AVAILABLE
                (C = 0)
                    Regard the residue(s) with the highest number of counts as
                    the consensus residue.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.)
                
                PERCENTAGE
                (0 < C < 1)
                    Regard the residue(s) which account for a @cutoff fraction
                    or more of the all counts at that position. If the singular
                    output mode is specified, residues with more than one such
                    consensus residue will be regarded as a blank. (No consensus
                    residue at that position.)
                
                PERFECT MATCH
                (C = 1)
                    Only regard a residue as the consensus residue if all
                    residues at that position are that residue.
                
                MINIMUM COUNT
                (C > 1)
                    Only regard a residue as a/the consensus residue if there
                    are @cutoff or more counts of that residue at that position.
                    If the singular output mode is specified, residues with more
                    than one such consensus residue will be regarded as a blank.
                    (No consensus residue at that position.
    
    @min_coverage
            (int)
            The minimum coverage required at a position in order for that
            position to have a consensus residue.
    
    @seq_type
            (int/str)
            An integer or string denoting the type of sequences that will be
            inputted.
            Accepted values for specifying DNA:
                SEQUENCE_TYPE.DNA
                1
                Any capitalization of "DNA" or it's abbreviations.
            Accepted values for specifying Amino Acid residues:
                SEQUENCE_TYPE.AMINO
                2
                Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                "Protein", or their abbreviations.
    
    @output_mode
            (int)
            An integer which specifies the format of the output consensus
            sequence:
                1:  SINGULAR
                    Only one residue at each position.
                2:  STANDARD
                    Accepts multiple residues possible for each position.
                    (Ex. "A[AT]GC")
                3:  AMBIGUOUS
                    (DNA only)
                    Only one character at each position, but that character can
                    represent multiple nucleotides.)
    
    Consensus_From_Standard_Matrix(list<list<int>>, int/float, int, int/str,
            int) -> str
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
    Take a list of sequences and return a count of the number of each residue at
    each position in the form of a list of dictionaries, where each dictionary
    corresponds to its sequence position.
    
    For sequences which denote more than one possible residue per location,
    (in the same input sequence) all counts are multiplied (12 for nucleotides,
    232792560 for amino acids) so that all partial counts are whole numbers
    while remaining fractionally proportional, and then optionally renormalized.
    (Ex. For DNA, if there's only one nucleotide at a given position, it gives
    12 "counts" for that position. If there's two possible nucleotides at a
    given position, it gives 6 "counts for that position. Three possible
    nucleotides, 4. Four possible nucleotides, 3. If normalization is enabled,
    then once the counting is finished, all counts are divided by 12. Using
    scaled whole counts and then normalizing at the end, rather than using
    fractional counts throughout the whole process, is done to minimize
    inaccuracy from Python's floating point limitations.)
    
    If an existing list of dictionaries is supplied, that list of dictionaries
    will be updated with the new counts. If no existing list of dictionaries is
    supplied, an empty one will be created.
    
    @sequences
            (list<str>)
            A list of biological sequences. All sequences need to be the same
            length. Unknown residues ("N" and "X" for DNA, "X" for amino acids)
            are accepted, as are gaps. ("_", "-", " ", ".") For DNA sequences,
            ambiguous nucleotides (RYSWMKBDHVNX) are also accepted.
            All sequences must be the same length. If an existing list of
            dictionaries is supplied, all sequences must also be the same length
            as said list.
    
    @seq_type
            (int/str)
            An integer or string denoting the type of sequences that will be
            inputted.
            Accepted values for specifying DNA:
                SEQUENCE_TYPE.DNA
                1
                Any capitalization of "DNA" or it's abbreviations.
            Accepted values for specifying Amino Acid residues:
                SEQUENCE_TYPE.AMINO
                2
                Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                "Protein", or their abbreviations.
    
    @input_mode
            (int)
            An integer which specifies the format of the input sequences:
                1:  SINGULAR
                    Only one residue at each position.
                2:  STANDARD
                    Accepts multiple residues possible for each position.
                    (Ex. "A[AT]GC") Also acceptes repeats of the same residue
                    being denoted in short hand. (Ex. "AT(3)GC(2)" = "ATTTGCC")
                3:  AMBIGUOUS
                    (DNA only)
                    Only one character at each position, but that character
                    can represent multiple nucleotides.
    @normalize
            (bool)
            Whether or not to renormalize the counts.
    
    @existing_dicts
            (list<dict<str:int/float>>)
            An existing list of residue count dictionaries can be supplied, and
            the new counts will be added to said dictionary. Otherwise, a
            fresh list of dictionaries (0 counts) will be created and used
            instead.
    
    Sequence_List_To_Counts_Dict(list<str>, int/str, int, bool,
            list<dict<str:int/float>>) -> list<dict<str:int/float>>
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
    Return a count matrix derived from count dictionaries.
    
    The input is a list of dictionaries. Each dictionary contains the residue
    counts corresponding to that position. Only standard Nucleotide or standard
    Amino Acid residues (ACGT or ACDEFGHIKLMNPQRSTVWY) will be counted.
    
    Ambiguous residues ("N" for nucleotides, "X" for amino acids) will also be
    accepted.
    
    The input dictionaries do not need to be "complete". (i.e., residues with a
    count of 0 don't need to be explicitly listed in the dictionaries)
    
    The returned matrix is a list of lists. Each list contains the residue
    counts corresponding to that position, and the counts are ordered
    alphabetically, (Ex. for nucleotide sequence matrices, the first integer is
    the number of Adenosines and the second integer is the number of Cytosines)
    with the final integer corresponding to the number of unknown residues.

    @counts
            (list<dict<str:int/float>>)
            A list of dictionaries containing the residue counts at the
            corresponding positions.
    
    @seq_type
            (int/str)
            An integer or string denoting the type of sequences that will be
            inputted.
            Accepted values for specifying DNA:
                SEQUENCE_TYPE.DNA
                1
                Any capitalization of "DNA" or it's abbreviations.
            Accepted values for specifying Amino Acid residues:
                SEQUENCE_TYPE.AMINO
                2
                Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                "Protein", or their abbreviations.
    
    Counts_Dict_To_Standard_Matrix(list<dict<str:int/float>>, int/str) ->
            list<list<int>>
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
    Return a standard-format consensus sequence as a list of strings. Each
    string contains the possible basic nucleotides at the corresponding
    position.
    
    Multiple possible residues at a position are denoted by square brackets,
    while repeats of the same residue are denoted by curly brackets.

    E.x:
        "A[TG]C(3)A"
    ->
        ["A", "TG", "C", "C", "C", "A"]
    
    @string
            (str)
            The consensus sequence.
    
    Parse_Standard_Consensus_Sequence(str) -> list<str>
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
    Return a DNA sequence containing ambiguous nucleotides (RYSWMKBDHVNX) as
    a list of strings. Each string contains the possible basic nucleotides at
    the corresponding position.
    
    @string
            (str)
            The consensus sequence.
    
    Parse_Ambiguous_Sequence(str) -> list<str>
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
                result[i][c] += 12
            elif c not in LIST__fillers:
                raise Exception(STR__invalid_seq_char.format(c=c))
    # Normalize
    if normalize:
        for dictionary in result:
            for key in dictionary:
                dictionary[key] = dictionary[key]/12
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
                dictionary[key] = dictionary[key]/12.0
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
                dictionary[key] = dictionary[key]/12.0
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
                result[i][c] += 232792560
            elif c not in LIST__fillers:
                raise Exception(STR__invalid_seq_char.format(c=c))
    # Normalize
    if normalize:
        for dictionary in result:
            for key in dictionary:
                dictionary[key] = dictionary[key]/232792560
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
    A Class used to derive the consensus sequence from a series of nucleotide
    or amino acid sequences.
    
    Unlike functions which require all the input upfront, Consensus_Builder is
    designed to be intialized and then have the counts updated as sequences are
    added one at a time.
    
    Added sequences need to all be of the same length.
    
    EXAMPLE USAGE:
    
        >>> cb = Consensus_Builder("DNA", 1)
        >>> cb.Add_Seq("ATGC")
        >>> cb.Add_Seq("ATAT")
        >>> cb.Add_Seq("A-A-")
        >>> consensus_majority = cb.Find_Consensus(0.5, 2, 2)
        >>> consensus_absolute = cb.Find_Consensus(1, 2, 2)
        >>> consensus_majority_coverage3 = cb.Find_Consensus(0.5, 3, 2)
    """
    def __init__(self, seq_type, input_mode=INPUT_MODE.SINGULAR):
        """
        Initialize the Consensus_Builder object.
        
        @seq_type
                (int/str)
                An integer or string denoting the type of sequences that will be
                inputted.
                Accepted values for specifying DNA:
                    SEQUENCE_TYPE.DNA
                    1
                    Any capitalization of "DNA" or it's abbreviations.
                Accepted values for specifying Amino Acid residues:
                    SEQUENCE_TYPE.AMINO
                    2
                    Any capitalization of "Amino Acids", "Amino Acid", "Amino",
                    "Protein", or their abbreviations.
        @input_mode
                (int)
                (DEFAULT: 1)
                An integer which specifies the input sequence format and how
                they are to be counted:
                    1:  SINGULAR
                        Only one residue at each position. Unknown residues
                        ("N" for nucleotides, "X" for amino acids) count towards
                        coverage but not towards individual counts.
                    2:  STANDARD
                        Accepts multiple residues at each position in the
                        standard consensus sequence format. (i.e. "A[AT]GC")
                        Multiple residues at a position count as partial counts
                        for that position.
                    3:  AMBIGUOUS
                        (DNA only)
                        Only one character at each position, but that character
                        can represent an ambiguous nucleotide. An ambiguous
                        nucleotide at any particular position count as partial
                        counts for that position.
        
        __init__(str/int, int) -> None
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
        Update the counts dictionaries with another sequence.
        
        An input mode may be specified. (See Class documentation) If no input
        mode is specified, the input mode specified during initialization will
        be used.
        
        @seq
                (str)
                The sequence being added to the counts.
        @input_mode
                (int)
                An integer which specifies the input sequence format and how
                they are to be counted:
                    1:  SINGULAR
                        Only one residue at each position. Unknown residues
                        ("N" for nucleotides, "X" for amino acids) count towards
                        coverage but not towards individual counts.
                    2:  STANDARD
                        Accepts multiple residues at each position in the
                        standard consensus sequence format. (Ex. "A[AT]GC")
                        Multiple residues at a position count as partial counts
                        for that position.
                        Also accepts repeats. (Ex. "AN(3)T" = "ANNNT")
                    3:  AMBIGUOUS
                        (DNA only)
                        Only one character at each position, but that character
                        can represent an ambiguous nucleotide. An ambiguous
                        nucleotide at any particular position count as partial
                        counts for that position.
        
        Add_Seq(str) -> None
        Add_Seq(str, int) -> None
        """
        if not input_mode: input_mode = self._mode
        self._counts_dicts = Sequence_List_To_Counts_Dict([seq], self._type,
                input_mode, False, self._counts_dicts)
    
    def Find_Consensus(self, cutoff, min_coverage,
            output_mode=OUTPUT_MODE.STANDARD):
        """
        Calculate the consensus sequence according to the given parameters.
        
        @cutoff
                (int/float)
                The cutoff used to determine the qualification for a residue to
                be regarded as the consensus residue. Different cutoff methods
                will be used depending on the number specified:
                
                    MISMATCH LIMIT
                    (C < 0)
                        Only regard a residue as a/the consensus residue if
                        there are |@cutoff| or fewer mismatches across all
                        counts.
                        If the singular output mode is specified, residues with
                        more than one such consensus residue will be regarded as
                        a blank. (No consensus residue at that position.)
                    
                    BEST AVAILABLE
                    (C = 0)
                        Regard the residue(s) with the highest number of counts
                        as the consensus residue.
                        If the singular output mode is specified, residues with
                        more than one such consensus residue will be regarded as
                        a blank. (No consensus residue at that position.)
                    
                    PERCENTAGE
                    (0 < C < 1)
                        Regard the residue(s) which account for a @cutoff
                        fraction or more of the all counts at that position.
                        If the singular output mode is specified, residues with
                        more than one such consensus residue will be regarded as
                        a blank. (No consensus residue at that position.)
                    
                    PERFECT MATCH
                    (C = 1)
                        Only regard a residue as the consensus residue if all
                        residues at that position are that residue.
                
                    MINIMUM COUNT
                    (C > 1)
                        Only regard a residue as a/the consensus residue if
                        there are @cutoff or more counts of that residue at that
                        position.
                        If the singular output mode is specified, residues with
                        more than one such consensus residue will be regarded as
                        a blank. (No consensus residue at that position.)
        
        @min_coverage
                (int)
                The minimum coverage required at a position in order for that
                position to have a consensus residue.
        
        @output_mode
                (int)
                An integer which specifies the format of the output consensus
                sequence:
                    1:  SINGULAR
                        Only one residue at each position.
                    2:  STANDARD
                        Accepts multiple residues possible for each position.
                        (Ex. "A[AT]GC")
                    3:  AMBIGUOUS
                        (DNA only)
                        Only one character at each position, but that character
                        can represent multiple nucleotides.
        """
        self.Calculate_Normalized()
        consensus = Consensus_From_Dicts(self._counts_dicts, cutoff,
                min_coverage, self._type, output_mode)
        return consensus
    
    def Calculate_Normalized(self):
        """
        THIS METHOD SHOULD ONLY BE CALLED INTERNALLY
        
        Calculate a "normalized" version of the counts dictionary.
        
        Whole nucleotide counts are multiplied by 12, while whole amino acid
        residue counts are multiplied by 232792560. Partial counts are
        multiplied by a fraction of this, depending on the number of
        possibilities.
        
        Normalization divides all counts by either 12 (for DNA) or 232792560
        (for Amino Acids) to get the "true" count.
        
        By adding whole numbers and normalizing at the end, instead of adding
        fractions to the count, the margin of error from Python floating point
        limitations is minimized. 
        """
        # Copy
        copy = []
        for d in self._counts_dicts: copy.append(dict(d))
        # If 
        if self._mode != INPUT_MODE.SINGULAR:
            # Determine factor
            if self._type == SEQUENCE_TYPE.DNA or self._type in LIST__dna:
                factor = 12.0
            else: factor = 232792560.0
            # Normalize
            for dictionary in copy:
                for key in dictionary:
                    dictionary[key] = dictionary[key]/factor
        # Update
        self._normalized = copy



# Controlled Print Statements ##################################################



# Main Loop ####################################################################


