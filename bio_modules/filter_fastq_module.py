""" Module with auxiliary functions for bio_tools.py, filter_fastq()"""


def quality_threshold(seq: str) -> float:
    total_quality = 0
    for symbol in seq:
        symbol_quality = ord(symbol) - 33
        total_quality += symbol_quality
    mean_quality = total_quality / len(seq)
    return mean_quality


def var_converse(var: (tuple, int, float)) -> tuple:
    if isinstance(var, (int, float)):
        var = (0, var)
    return var
