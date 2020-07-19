#!/usr/bin/env python3

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
import click
from pathlib import Path
from shutil import get_terminal_size
from icecream import ic
from kcl.configops import click_read_config
from kcl.configops import click_write_config_entry
from kcl.inputops import input_iterator


ic.configureOutput(includeContext=True)
ic.lineWrapWidth, _ = get_terminal_size((80, 20))
# import IPython; IPython.embed()
# import pdb; pdb.set_trace()
# from pudb import set_trace; set_trace(paused=False)

APP_NAME = 'structuredatafile'
# https://stackoverflow.com/questions/14921929/python-progam-to-read-sdf-chemistry-file

# DONT CHANGE FUNC NAME
@click.command()
@click.argument("paths", type=str, nargs=-1)
@click.argument("sysskel",
                type=click.Path(exists=False,
                                dir_okay=True,
                                file_okay=False,
                                path_type=str,
                                allow_dash=False), nargs=1, required=True)
@click.option('--add', is_flag=True)
@click.option('--verbose', is_flag=True)
@click.option('--debug', is_flag=True)
@click.option('--ipython', is_flag=True)
@click.option("--null", is_flag=True)
#@click.group()
def cli(paths,
        sysskel,
        add,
        verbose,
        debug,
        ipython,
        null):

    byte = b'\n'
    if null:
        byte = b'\x00'

    config, config_mtime = click_read_config(click_instance=click,
                                             app_name=APP_NAME,
                                             verbose=verbose)
    if verbose:
        ic(config, config_mtime)

    for index, path in enumerate(input_iterator(strings=paths,
                                                null=null,
                                                verbose=verbose)):
        if verbose:
            ic(index, path)

        with open(path, 'rb') as fh:
            path_bytes_data = fh.read()

        if ipython:
            import IPython; IPython.embed()

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

