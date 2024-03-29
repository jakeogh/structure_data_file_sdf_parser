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

import pprint

import click
from configtool import click_read_config
from configtool import click_write_config_entry
from enumerate_input import enumerate_input
from openbabel import pybel

# import IPython; IPython.embed()
# import pdb; pdb.set_trace()
# from pudb import set_trace; set_trace(paused=False)

APP_NAME = 'structure_data_file_sdf_parser'
# https://stackoverflow.com/questions/14921929/python-progam-to-read-sdf-chemistry-file

# can read gzipped files
def molecule_dict_generator(path, verbose=False):
    for mol in pybel.readfile('sdf', path):
        yield dict(mol.data)


@click.command()
@click.argument("paths", type=str, nargs=-1)
@click.option('--add', is_flag=True)
@click.option('--verbose', is_flag=True)
@click.option('--debug', is_flag=True)
@click.option('--ipython', is_flag=True)
@click.option("--null", is_flag=True)
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

        for mol_data in molecule_dict_generator(path, verbose=verbose):
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

