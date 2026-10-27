"""Make PROPKA 3.5.1 run on Python 3.14.

PROPKA's parameter-file reader does `self.__annotations__.get(name)` to decide how to parse
each line. Under PEP 649 (Python 3.14) a class's `__annotations__` is served by a descriptor
on `type` rather than sitting in the class `__dict__`, so an *instance* can no longer see it
and every PROPKA run dies before reading its own parameters. Assigning the dict back onto the
class does not help: 3.14 redirects that write to `__annotations_cache__`.

So the one method that needs the annotations gets them passed in instead. Everything else in
the library is untouched, and the installed package is not modified - the patch is applied at
import time by whatever imports this module.

Import this before `propka.run`, then call `check()` once to confirm the patch took.
"""
import inspect

from propka.parameters import Parameters
import propka.parameters as _pp

#: Raw, UNevaluated annotations, captured once from the class. Not eval_str=True: PROPKA
#: compares against constants that are themselves a mix of typing aliases and plain strings
#: (_T_MATRIX is the string 'InteractionMatrix'), so evaluating turns the strings into
#: classes and every identity comparison then fails silently into parse_parameter.
_ANN = inspect.get_annotations(Parameters)


def _parse_line(self, line):
    """PROPKA's own parse_line, with the annotation lookup fixed."""
    comment_pos = line.find("#")
    if comment_pos != -1:
        line = line[:comment_pos]
    words = line.split()
    if not words:
        return
    ann = _ANN.get(words[0])
    if ann is _pp._T_NUMBER_DICTIONARY:
        self.parse_to_number_dictionary(words)
    elif ann is _pp._T_STRING_LIST:
        self.parse_to_string_list(words)
    elif ann is _pp._T_STRING:
        self.parse_string(words)
    elif ann is _pp._T_LIST_DICTIONARY:
        self.parse_to_list_dictionary(words)
    elif ann is _pp._T_MATRIX or ann is _pp._T_PAIR_WISE_MATRIX:
        self.parse_to_matrix(words)
    elif ann is _pp._T_STRING_DICTIONARY:
        self.parse_to_string_dictionary(words)
    else:
        self.parse_parameter(words)


Parameters.parse_line = _parse_line


def check():
    """Prove the patch works on PROPKA's own parameter file, not on a toy input.

    If the parameters failed to load, the dictionaries below come back empty and every pKa
    PROPKA then reports would be silently wrong rather than absent - which is the failure
    mode worth guarding against.
    """
    from propka.input import read_parameter_file
    p = read_parameter_file("propka.cfg", Parameters())
    n_model = len(getattr(p, "model_pkas", {}) or {})
    if n_model < 5:
        raise SystemExit("propka parameters did not load (%d model pKa values); "
                         "the compatibility patch is not working" % n_model)
    return n_model
