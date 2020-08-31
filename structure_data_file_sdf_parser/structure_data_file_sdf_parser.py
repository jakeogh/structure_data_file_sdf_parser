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
import gzip
from openbabel import pybel
import click
from pathlib import Path
from icecream import ic
from kcl.configops import click_read_config
from kcl.configops import click_write_config_entry
from kcl.inputops import enumerate_input
from kcl.commandops import run_command

ic.configureOutput(includeContext=True)
# import IPython; IPython.embed()
# import pdb; pdb.set_trace()
# from pudb import set_trace; set_trace(paused=False)

APP_NAME = 'structure_data_file_sdf_parser'
# https://stackoverflow.com/questions/14921929/python-progam-to-read-sdf-chemistry-file

def readfile(format, filename, opt=None):
    """Iterate over the molecules in a file.

    Required parameters:
       format - see the informats variable for a list of available
                input formats
       filename

    Optional parameters:
       opt    - a dictionary of format-specific options
                For format options with no parameters, specify the
                value as None.

    You can access the first molecule in a file using the next() method
    of the iterator (or the next() keyword in Python 3):
        mol = readfile("smi", "myfile.smi").next() # Python 2
        mol = next(readfile("smi", "myfile.smi"))  # Python 3

    You can make a list of the molecules in a file using:
        mols = list(readfile("smi", "myfile.smi"))

    You can iterate over the molecules in a file as shown in the
    following code snippet:
    >>> atomtotal = 0
    >>> for mol in readfile("sdf", "head.sdf"):
    ...     atomtotal += len(mol.atoms)
    ...
    >>> print atomtotal
    43
    """
    if opt is None:
        opt = {}
    obconversion = ob.OBConversion()
    formatok = obconversion.SetInFormat(format)
    for k, v in opt.items():
        if v is None:
            obconversion.AddOption(k, obconversion.INOPTIONS)
        else:
            obconversion.AddOption(k, obconversion.INOPTIONS, str(v))
    if not formatok:
        raise ValueError("%s is not a recognised Open Babel format" % format)
    #if not os.path.isfile(filename):
    #    raise IOError("No such file: '%s'" % filename)

    def filereader():
        obmol = ob.OBMol()
        notatend = obconversion.ReadFile(obmol, filename)
        while notatend:
            yield Molecule(obmol)
            obmol = ob.OBMol()
            notatend = obconversion.Read(obmol)
    return filereader()


pybel.readfile = readfile




def molecule_dict_generator(path, verbose=False):
    if path.endswith('.gz'):
        tmpdir = tempfile.mkdtemp()

        # https://github.com/cybernoid/archivemount/issues/16
        command = "archivemount {path} {mountpoint} -o auto_unmount -o readonly".format(path=path, mountpoint=tmpdir)
        output = run_command(command=command,
                             shell=True,
                             verbose=verbose,
                             expected_exit_code=0)
        ic(output)
        import IPython; IPython.embed()
        #with gzip.open(path) as gfh:
        #    if verbose:
        #        ic(gfh)
        #    with temp_fifo(verbose=verbose) as fifo_file:
        #        with open(fifo_file, 'wb') as ffh:
        #            sdf_chunk = gfh.read(4096*4)
        #            ic(len(sdf_chunk))
        #            ffh.write(sdf_chunk)
        #            for mol in pybel.readfile('sdf', fifo_file):
        #                sdf_chunk = gfh.read(4096*4)
        #                ic(len(sdf_chunk))
        #                ffh.write(sdf_chunk)
        #                yield dict(mol.data)
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

