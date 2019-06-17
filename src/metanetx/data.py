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

import json
import logging
import os
import re
from collections import defaultdict
from gzip import GzipFile
from io import BytesIO, TextIOWrapper

import requests


logger = logging.getLogger(__name__)


class Compartment:
    def __init__(self, mnx_id, name, xref):
        self.mnx_id = mnx_id
        self.name = name
        self.xref = xref

    @property
    def annotation(self):
        return compartment_xrefs.get(self.mnx_id, {})


class Reaction:
    metabolite_regex = re.compile(r"([\d|\.]+) (\w+)@(\w+)")

    def __init__(self, mnx_id, name, equation_string, ec):
        self.mnx_id = mnx_id
        self.name = name
        self.equation_string = equation_string
        self.equation_parsed = Reaction.parse_equation(self.equation_string)
        self.ec = ec

    @property
    def annotation(self):
        return reaction_xrefs.get(self.mnx_id, {})

    @staticmethod
    def parse_equation(equation_string):
        """
        Parse the reactions equation string.

        Returns a list of all metabolites used in the equation, with the
        following keys:

            metabolite_id: The metanetx id of the metabolite
            compartment_id: The metanetx id of the compartment
            coefficient: The stoichiometry coefficient (negative if substrate,
                positive if product)
        """
        equation = []
        substrates, products = equation_string.split(" = ")
        for match in re.findall(Reaction.metabolite_regex, substrates):
            coefficient, metabolite_id, compartment_id = match
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": float(coefficient) * -1,
                }
            )
        for match in re.findall(Reaction.metabolite_regex, products):
            coefficient, metabolite_id, compartment_id = match
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": float(coefficient),
                }
            )
        return equation


class Metabolite:
    def __init__(self, mnx_id, name):
        self.mnx_id = mnx_id
        self.name = name

    @property
    def annotation(self):
        return metabolite_xrefs.get(self.mnx_id, {})


# MetaNetX data will be read into memory into the following dicts, keyed by ID.
compartments = {}
reactions = {}
metabolites = {}

# Cross-references
compartment_xrefs = {}
reaction_xrefs = {}
metabolite_xrefs = {}

# Reaction names are collected separately and merged
reaction_names = {}


def load_metanetx_data():
    with _retrieve("reaction_names.json") as file_:
        reaction_names = json.load(file_)
    logger.info(f"Loaded {len(reaction_names)} reaction name mappings")

    for line in _iterate_tsv(_retrieve("comp_prop.tsv")):
        mnx_id, name, xref = line.rstrip("\n").split("\t")
        compartments[mnx_id] = Compartment(mnx_id, name, xref)
    logger.info(f"Loaded {len(compartments)} compartments")

    for line in _iterate_tsv(_retrieve("comp_xref.tsv")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in compartment_xrefs:
                compartment_xrefs[mnx_id] = defaultdict(list)
            compartment_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(compartment_xrefs)} compartment cross-references")

    for line in _iterate_tsv(_retrieve("reac_prop.tsv")):
        mnx_id, equation, _, _, ec, _ = line.rstrip("\n").split("\t")
        name = reaction_names.get(mnx_id)
        reactions[mnx_id] = Reaction(mnx_id, name, equation, ec)
    logger.info(f"Loaded {len(reactions)} reactions")

    for line in _iterate_tsv(_retrieve("reac_xref.tsv")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in reaction_xrefs:
                reaction_xrefs[mnx_id] = defaultdict(list)
            reaction_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(reaction_xrefs)} reaction cross-references")

    for line in _iterate_tsv(_retrieve("chem_prop.tsv")):
        mnx_id, name, _, _, _, _, _, _, _ = line.rstrip("\n").split("\t")
        metabolites[mnx_id] = Metabolite(mnx_id, name)
    logger.info(f"Loaded {len(metabolites)} metabolites")

    for line in _iterate_tsv(_retrieve("chem_xref.tsv")):
        xref, mnx_id, _, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in metabolite_xrefs:
                metabolite_xrefs[mnx_id] = defaultdict(list)
            metabolite_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(metabolite_xrefs)} metabolite cross-references")


def _retrieve(filename):
    if os.environ.get("LOCAL_METANETX_DATA"):
        logger.debug(f"Reading data/{filename}")
        return open(f"data/{filename}")
    else:
        logger.debug(f"Downloading {filename} from external storage")
        r = requests.get(
            f"https://storage.googleapis.com/dd-decaf/metanetx/{filename}.gz"
        )
        return TextIOWrapper(GzipFile(fileobj=BytesIO(r.content)))


def _iterate_tsv(file_):
    with file_:
        for line in file_:
            if line.startswith("#"):
                continue
            yield line
