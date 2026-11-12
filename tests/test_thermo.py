import json
import math
import os
import random

import pytest
from Bio.SeqUtils import MeltingTemp as mt

from primerJinn.multiplex import thermo

HERE = os.path.dirname(os.path.abspath(__file__))
NEB = os.path.join(HERE, "..", "docs", "research", "reference_code", "neb_calib.json")


def test_reference_self_tests():
    thermo._self_test(verbose=False)


def test_neb_scale_reproduces_neb_calculator():
    data = json.load(open(NEB))
    for conc, values in data["neb"].items():
        nM = float(conc) * 1000
        for seq, neb in zip(data["seqs"], values):
            assert abs(thermo.tm_neb(seq, nM) - neb) < 0.01, (seq, conc)


def test_physical_matches_biopython():
    rng = random.Random(1)
    for _ in range(100):
        n = rng.randint(18, 30)
        s = "".join(rng.choice("ACGT") for _ in range(n))
        c = thermo.Conditions(k_mM=50, mg_mM=2, dntp_mM=0.8, primer_nM=200)
        ours = thermo.analyze(s, c).tm_C
        bio = mt.Tm_NN(s, nn_table=mt.DNA_NN4, Na=0, K=50, Tris=0, Mg=2, dNTPs=0.8,
                       saltcorr=7, dnac1=200, dnac2=0)
        assert abs(ours - bio) < 0.01


def test_theta_and_structure_cost():
    c = thermo.Conditions(primer_nM=100)
    r = thermo.analyze("TACTTCCAGTGCTCAGCGTA", c)
    assert abs(thermo.theta_from_dg(r.dG(r.tm_C), r.tm_C, 100) - 0.5) < 1e-6
    # a -1 kcal/mol hairpin costs about 1.13 kcal/mol at 60 C (fact-check value)
    assert abs(thermo.structure_cost(-1.0, 60.0) - 1.13) < 0.01
    assert abs(thermo.structure_cost(0.0, 60.0) - 0.4589) < 0.001
    assert thermo.structure_cost(5.0, 60.0) < 0.01
    assert math.isfinite(thermo.structure_cost(-50.0, 60.0))
