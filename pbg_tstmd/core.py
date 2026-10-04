"""Explicit core builder so the workspace's own Process registers regardless of
how the editable install was made (see viva-superpowers discovery conventions)."""
from process_bigraph import allocate_core

from .processes import TstmdEcmProcess


def build_core(core=None):
    if core is None:
        core = allocate_core()
    core.register_link("TstmdEcmProcess", TstmdEcmProcess)
    return core
