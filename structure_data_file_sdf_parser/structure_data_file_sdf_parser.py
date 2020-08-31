#!/usr/bin/env python3.7

# pylint: disable=C0111     # docstrings are always outdated and wrong
# pylint: disable=W0511     # todo is encouraged
# pylint: disable=R0902     # too many instance attributes
# pylint: disable=C0302     # too many lines in module
# pylint: disable=C0103     # single letter var names
# pylint: disable=R0911     # too many return statements
# pylint: disable=R0912     # too many branches
# pylint: disable=R0915     # too many statements
# pylint: disable=R0913     # too many arguments
# pylint: disable=R1702     # too many nested blocks
# pylint: disable=R0914     # too many local variables
# pylint: disable=R0903     # too few public methods
# pylint: disable=E1101     # no member for base
# pylint: disable=W0201     # attribute defined outside __init__
## pylint: disable=W0703     # catching too general exception

import os
import sys
import pprint
import tempfile
from openbabel import pybel
import click
from pathlib import Path
from icecream import ic
from kcl.configops import click_read_config
from kcl.configops import click_write_config_entry
from kcl.inputops import enumerate_input
from contextlib import contextmanager


ic.configureOutput(includeContext=True)
# import IPython; IPython.embed()
# import pdb; pdb.set_trace()
# from pudb import set_trace; set_trace(paused=False)

APP_NAME = 'structure_data_file_sdf_parser'
# https://stackoverflow.com/questions/14921929/python-progam-to-read-sdf-chemistry-file

# https://stackoverflow.com/questions/1430446/create-a-temporary-fifo-named-pipe-in-python
@contextmanager
def temp_fifo(verbose=False):
    """Context Manager for creating named pipes with temporary names."""
    tmpdir = tempfile.mkdtemp()
    filename = os.path.join(tmpdir, 'fifo')  # Temporary filename
    if verbose:
        ic(filename)
    os.mkfifo(filename)  # Create FIFO
    try:
        yield filename
    finally:
        os.unlink(filename)  # Remove file
        os.rmdir(tmpdir)  # Remove directory


def molecule_dict_generator(path, verbose=False):
    if path.endswith('.gz'):
        with temp_fifo(verbose=verbose) as fifo_file:
            for mol in pybel.readfile('sdf', fifo_file):
                yield dict(mol.data)
    else:
        for mol in pybel.readfile('sdf', path):
            yield dict(mol.data)


# DONT CHANGE FUNC NAME
@click.command()
@click.argument("paths", type=str, nargs=-1)
@click.option('--add', is_flag=True)
@click.option('--verbose', is_flag=True)
@click.option('--debug', is_flag=True)
@click.option('--ipython', is_flag=True)
@click.option("--null", is_flag=True)
#@click.group()
def cli(paths,
        add,
        verbose,
        debug,
        ipython,
        null):

    config, config_mtime = click_read_config(click_instance=click,
                                             app_name=APP_NAME,
                                             verbose=verbose)
    if verbose:
        ic(config, config_mtime)

    for index, path in enumerate_input(iterator=paths,
                                       null=null,
                                       debug=debug,
                                       verbose=verbose):
        if verbose:
            ic(index, path)

        for mol_data in molecule_dict_generator(path):
            pprint.pprint(mol_data, indent=1)
            if ipython:
                import IPython; IPython.embed()
                break

        if add:
            section = "test_section"
            key = "test_key"
            value = "test_value"
            config, config_mtime = click_write_config_entry(click_instance=click,
                                                            app_name=APP_NAME,
                                                            section=section,
                                                            key=key,
                                                            value=value,
                                                            verbose=verbose)
            if verbose:
                ic(config)

