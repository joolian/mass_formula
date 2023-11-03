import mass_formula as ms


def test_ppm_error():
    assert ms.ppm_error(500, 1) == 2000


def test_ppm_to_mass():
    assert  ms.ppm_to_mass(500, 2000) == 1
