"""Simulation and testing of candidate solutions for ischaemic heart disease.

The package simulates *scenarios*, not patients. Every therapeutic effect is
traceable to an entry in ``data/solutions/trial_anchors.csv``, and benefit is only
credited when the evidence tier and the endpoint class permit it. Nothing here is
a clinical prediction, a treatment recommendation, or evidence of a cure.
"""

from ihd_solutions.evidence import CreditMode, EndpointClass, EvidenceTier

__all__ = ["CreditMode", "EndpointClass", "EvidenceTier", "__version__"]
__version__ = "1.0.0"
