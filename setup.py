# -*- coding: utf-8 -*-

import sys
import fastentrypoints
from setuptools import find_packages, setup
if not sys.version_info[0] == 3:
    sys.exit("Python 3 is required. Use: \'python3 setup.py install\'")

dependencies = ["icecream", "click", "colorama"]

config = {
    "version": "0.1",
    "name": "structure_data_file_sdf_parser",
    "url": "https://github.com/jakeogh/structure_data_file_sdf_parser",
    "license": "ISC",
    "author": "Justin Keogh",
    "author_email": "github.com@v6y.net",
    "description": "Parse (chemical) Strcuture Data Files (SDF)",
    "long_description": __doc__,
    "packages": find_packages(exclude=['tests']),
    "include_package_data": True,
    "zip_safe": False,
    "platforms": "any",
    "install_requires": dependencies,
    "entry_points": {
        "console_scripts": [
            "structure_data_file_sdf_parser=structure_data_file_sdf_parser.structure_data_file_sdf_parser:cli",
        ],
    },
}

setup(**config)
