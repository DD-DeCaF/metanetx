# Copyright (c) 2019, Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Data classes for MetaNetX data."""

import logging
import re
import os
from gzip import GzipFile
from io import BytesIO, TextIOWrapper

import requests


logger = logging.getLogger(__name__)


class Compartment:
    def __init__(self, mnx_id, name, xref):
        self.mnx_id = mnx_id
        self.name = name
        self.xref = xref


class Reaction:
    def __init__(self, mnx_id, equation, ec):
        self.mnx_id = mnx_id
        self.equation = equation
        self.ec = ec

    def parse_equation(self):
        """
        Returns a list of all metabolites used in the equation, with the
        following keys:

            metabolite_id: The metanetx id of the metabolite
            compartment_id: The metanetx id of the compartment
            coefficient: The stoichiometry coefficient (negative if substrate,
                positive if product)
        """
        equation = []
        substrates, products = self.equation.split(" = ")
        metabolite_regex = r"(\d+) (\w+)@(\w+)"
        for match in re.findall(metabolite_regex, substrates):
            coefficient, metabolite_id, compartment_id = match
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": int(coefficient) * -1,
                }
            )
        for match in re.findall(metabolite_regex, products):
            coefficient, metabolite_id, compartment_id = match
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": int(coefficient),
                }
            )
        return equation


class Metabolite:
    def __init__(self, mnx_id, description):
        self.mnx_id = mnx_id
        self.description = description


# MetaNetX data will be read into memory into the following dicts, keyed by ID.
compartments = {}
reactions = {}
metabolites = {}


def load_metanetx_data():
    for line in _retrieve("comp_prop.tsv"):
        if line.startswith("#"):
            continue
        mnx_id, name, xref = line.rstrip("\n").split("\t")
        compartments[mnx_id] = Compartment(mnx_id, name, xref)
    logger.info(f"Loaded {len(compartments)} compartments")

    for line in _retrieve("reac_prop.tsv"):
        if line.startswith("#"):
            continue
        mnx_id, equation, _, _, ec, _ = line.rstrip("\n").split("\t")
        reactions[mnx_id] = Reaction(mnx_id, equation, ec)
    logger.info(f"Loaded {len(reactions)} reactions")

    for line in _retrieve("chem_prop.tsv"):
        if line.startswith("#"):
            continue
        mnx_id, description, _, _, _, _, _, _, _ = line.rstrip("\n").split("\t")
        metabolites[mnx_id] = Metabolite(mnx_id, description)
    logger.info(f"Loaded {len(metabolites)} metabolites")


def _retrieve(filename):
    if os.environ.get("LOCAL_METANETX_DATA"):
        logger.debug(f"Reading data/{filename}")
        with open(f"data/{filename}") as file_:
            for line in file_:
                yield line
    else:
        logger.debug(f"Downloading {filename} from external storage")
        r = requests.get(
            f"https://storage.googleapis.com/dd-decaf/metanetx/{filename}.gz"
        )
        with TextIOWrapper(GzipFile(fileobj=BytesIO(r.content))) as file_:
            for line in file_:
                yield line
