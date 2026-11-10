"""plab -- an in-silico protein design sandbox for leukemia surface targets.

Scope and honesty note
----------------------
This package generates candidate binder sequences against a real,
experimentally-determined target surface and scores them with published
biophysical formulas plus a protein language model. It is a *design and
triage* tool, of the kind used to decide what to build first in a lab.

It does not predict clinical efficacy. No score produced here means a
molecule will help a patient. Validating any candidate requires expressing
the protein, measuring binding, testing on primary cells, then animal models
and trials. The composite score is reported relative to a clinically
characterised reference binder so that the number has a defined meaning.
"""

__version__ = "0.1.0"
