import pytest


def pytest_configure(config):
    config.addinivalue_line("markers", "bitlevel: check that results are reproducible at the bit level")
    config.addinivalue_line("markers", "requires_hyperion_binaries: require the Hyperion Fortran binaries to run")


def pytest_addoption(parser):
    parser.addoption('--generate-reference', help="generate reference results for bit-level tests", type=str)
    parser.addoption('--enable-bit-level-tests', help="enable bit-level tests", action="store_true")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--enable-bit-level-tests"):
        return
    skip_bit_level = pytest.mark.skip(reason="need --enable-bit-level-tests option to run")
    for item in items:
        if "bitlevel" in item.keywords:
            item.add_marker(skip_bit_level)
