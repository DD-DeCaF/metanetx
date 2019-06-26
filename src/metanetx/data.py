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

import gzip
import json
import logging
import re
from collections import defaultdict


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
    metabolite_regex = re.compile(r"^([\d\.]+) (\w+)@(\w+)$")

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
        for substrate in substrates.split(" + "):
            match = re.match(Reaction.metabolite_regex, substrate)
            if not match:
                raise ValueError(f"Invalid metabolite format: {substrate}")
            coefficient, metabolite_id, compartment_id = match.groups()
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": float(coefficient) * -1,
                }
            )
        for product in products.split(" + "):
            match = re.match(Reaction.metabolite_regex, product)
            if not match:
                raise ValueError(f"Invalid metabolite format: {product}")
            coefficient, metabolite_id, compartment_id = match.groups()
            equation.append(
                {
                    "metabolite_id": metabolite_id,
                    "compartment_id": compartment_id,
                    "coefficient": float(coefficient),
                }
            )
        return equation


class Metabolite:
    def __init__(self, mnx_id, name, formula):
        self.mnx_id = mnx_id
        self.name = name
        self.formula = formula

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
    with gzip.open("data/reaction_names.json.gz", "rt") as file_:
        reaction_names = json.load(file_)
    logger.info(f"Loaded {len(reaction_names)} reaction name mappings")

    for line in _iterate_tsv(gzip.open("data/comp_prop.tsv.gz", "rt")):
        mnx_id, name, xref = line.rstrip("\n").split("\t")
        compartments[mnx_id] = Compartment(mnx_id, name, xref)
    logger.info(f"Loaded {len(compartments)} compartments")

    for line in _iterate_tsv(gzip.open("data/comp_xref.tsv.gz", "rt")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in compartment_xrefs:
                compartment_xrefs[mnx_id] = defaultdict(list)
            compartment_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(compartment_xrefs)} compartment cross-references")

    filtered_reaction_count = 0
    for line in _iterate_tsv(gzip.open("data/reac_prop.tsv.gz", "rt")):
        mnx_id, equation, _, _, ec, _ = line.rstrip("\n").split("\t")
        name = reaction_names.get(mnx_id)
        try:
            reactions[mnx_id] = Reaction(mnx_id, name, equation, ec)
        except ValueError:
            # At the time of this writing, this is expected to amount up to 315
            # reactions which have an inexact stoichiometry coefficient "(n)"
            # and hence are unusable in the context of a metabolic model.
            filtered_reaction_count += 1
    logger.info(
        f"Loaded {len(reactions)} reactions (filtered {filtered_reaction_count}"
        " unparseable equations)"
    )

    for line in _iterate_tsv(gzip.open("data/reac_xref.tsv.gz", "rt")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in reaction_xrefs:
                reaction_xrefs[mnx_id] = defaultdict(list)
            reaction_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(reaction_xrefs)} reaction cross-references")

    for line in _iterate_tsv(gzip.open("data/chem_prop.tsv.gz", "rt")):
        mnx_id, name, formula, _, _, _, _, _, _ = line.rstrip("\n").split("\t")
        metabolites[mnx_id] = Metabolite(mnx_id, name, formula)
    logger.info(f"Loaded {len(metabolites)} metabolites")

    for line in _iterate_tsv(gzip.open("data/chem_xref.tsv.gz", "rt")):
        xref, mnx_id, _, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            if mnx_id not in metabolite_xrefs:
                metabolite_xrefs[mnx_id] = defaultdict(list)
            metabolite_xrefs[mnx_id][namespace].append(reference)
    logger.info(f"Loaded {len(metabolite_xrefs)} metabolite cross-references")


def _iterate_tsv(file_):
    with file_:
        for line in file_:
            if line.startswith("#"):
                continue
            yield line
