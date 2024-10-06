""" This module helps to work with sequences of nucleic acids """


import dna_rna_module
import filter_fastq_module

dna = {"a", "t", "g", "c"}
rna = {"a", "u", "g", "c"}


def run_dna_rna_tools(*args: tuple) -> (str, int, float, bool, list):
    """Function transcribes, reverts,
    returns complement and reverse complement,
    finds gc percentage, determines palindromic condition
    and counts triplets of a list of string or a single string.
    Only one action for one launch is used.



    Args: *args: tuple of strings
    The last string should always be a function name



    Returns: (str, list): string, integer, float or bool
            or list of strings, integers, floats or bools depending on function


    """
    function = args[-1]
    strands = args[:-1]
    normal_output = []
    for strnd in strands:
        if (set(strnd.lower()).issubset(dna) or
           set(strnd.lower()).issubset(rna)):
            if function == "transcribe":
                normal_output.append(dna_rna_module.transcribe(strnd))
            if function == "reverse":
                normal_output.append(dna_rna_module.reverse(strnd))
            if function == "complement":
                normal_output.append(dna_rna_module.complement(strnd))
            if function == "reverse_complement":
                normal_output.append(dna_rna_module.reverse_complement(strnd))
            if function == "gc_content":
                normal_output.append(dna_rna_module.gc_content(strnd))
            if function == "is_palindrom":
                normal_output.append(dna_rna_module.is_palindrom(strnd))
            if function == "triplets_count":
                normal_output.append(dna_rna_module.triplets_count(strnd))
    if len(normal_output) == 1:
        normal_output = normal_output[0]
    return normal_output


def filter_fastq(
    seqs: dict(type=str),
    gc_bounds: (tuple, int, float) = (0, 100),
    length_bounds: (tuple, int, float) = (0, 2 ** 32),
    quality_threshold: int = 0,
) -> dict:
    """Function filters fastq sequences by their length, GC % and quality.


    Args:
        seqs (dict(type = str)): dict of string items for analysis
        gc_bounds (tuple, int, float) = (0, 100): the GC% for filtering;
            tuple, integer or float; by default is (0, 100).
            If int or float, it is assumed that this is the upper bound.
        length_bounds(tuple, int, float) = (0, 2**32):
            the length interval for filtering;
            everything is similar to gc_bounds, but by default is (0, 2**32).
        quality_threshold(int) = 0: integer,
            the threshold value of the average read quality,
            filtering is 0 by default (phred33 scale).


    Returns: dict: dictionary of strings


    """
    output_seqs = {}
    gc_bounds = filter_fastq_module.var_converse(gc_bounds)
    length_bounds = filter_fastq_module.var_converse(length_bounds)
    for key in seqs.keys():
        length = len(seqs[key][0])
        gc_content = dna_rna_module.gc_content(seqs[key][0])
        quality = filter_fastq_module.quality_threshold(seqs[key][1])
        if (
            gc_bounds[0] <= gc_content <= gc_bounds[1]
            and length_bounds[0] <= length <= length_bounds[1]
            and quality >= quality_threshold
        ):
            output_seqs[key] = seqs[key]
    return output_seqs
