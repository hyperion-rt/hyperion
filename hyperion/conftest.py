import os


def pytest_addoption(parser):
    parser.addoption('--generate-reference', help="generate reference results for bit-level tests", type="string")


def pytest_runtest_setup(item):
    if 'generate_reference' in item.keywords:
        if item.config.getvalue("generate_reference"):
            if item.config.getvalue("generate_reference").startswith('-'):
                raise Exception("Need to specify output directory for generating reference files")
            if not item.config.getvalue("generate_reference").startswith('/'):
                raise Exception("Need to specify output directory for generating reference files as an absolute path")
            item.funcargs['generate'] = item.config.getvalue("generate_reference")
