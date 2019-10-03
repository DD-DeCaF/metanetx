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

from fuzzywuzzy import fuzz


logger = logging.getLogger(__name__)


class Compartment:
    def __init__(self, mnx_id, name, xref):
        self.mnx_id = mnx_id
        self.name = name
        self.xref = xref
        self.annotation = defaultdict(list)


class Reaction:
    metabolite_regex = re.compile(r"^([\d\.]+) (\w+)@(\w+)$")

    def __init__(self, mnx_id, name, equation_string, ec):
        self.mnx_id = mnx_id
        self.name = name
        self.equation_string = equation_string
        self.equation_parsed = Reaction.parse_equation(self.equation_string)
        self.ec = ec
        self.annotation = defaultdict(list)

    def match(self, query):
        """
        Match an arbitrary search string against this reaction.

        Use Levenshtein Distance to match the query on ID (both MNX and
        annotations), name or EC number. Returns a number between 0 and 100, 100
        being a perfect match.
        """
        return max(
            [
                fuzz.ratio(query, self.mnx_id),
                # Select the highest match among all annotations
                max(
                    max(fuzz.ratio(query, identifier) for identifier in db)
                    for db in self.annotation.values()
                )
                # Handle reactions with no annotations
                if self.annotation else 0,
                fuzz.partial_ratio(query, self.name),
                fuzz.ratio(query, self.ec),
            ]
        )

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
        self.annotation = defaultdict(list)

    def match(self, query):
        """
        Match an arbitrary search string against this metabolite.

        Levenshtein Distance is a bit too slow here, so simply check for
        substring matches case insensitively. Returns True if match.
        """
        return (
            query.lower() in self.mnx_id.lower()
            or query.lower() in self.name.lower()
        )


# MetaNetX data will be read into memory into the following dicts, keyed by ID.
compartments = {}
reactions = {}
metabolites = {}

# Reaction names are collected separately and merged
reaction_names = {}

# MetaNetX doesn't use miriam namespace identifiers. We do, so map them
# correctly.
compartment_namespace_map = {
    "name": "name",  # unconfirmed
    "cco": "cco",
    "go": "go",
    "bigg": "bigg.compartment",
    "seed": "seed.reaction",
}
reaction_namespace_map = {
    "reactome": "reactome",
    "sabiork": "sabiork.reaction",
    "deprecated": "deprecated",
    "kegg": "kegg.reaction",
    "rhea": "rhea",
    "metacyc": "metacyc.reaction",
    "bigg": "bigg.reaction",
    "seed": "seed.reaction",
}
metabolite_namespace_map = {
    "reactome": "reactome",
    "sabiork": "sabiork.compound",
    "chebi": "chebi",
    "metacyc": "metacyc.compound",
    "envipath": "envipath",  # unconfirmed
    "deprecated": "deprecated",
    "slm": "slm",  # unconfirmed
    "kegg": "kegg.compound",
    "hmdb": "hmdb",
    "lipidmaps": "lipidmaps",
    "bigg": "bigg.metabolite",
    "seed": "seed.compound",
}


def load_metanetx_data():
    with gzip.open("data/reaction_names.json.gz", "rt") as file_:
        reaction_names = json.load(file_)
    logger.info(f"Loaded {len(reaction_names)} reaction name mappings")

    for line in _iterate_tsv(gzip.open("data/comp_prop.tsv.gz", "rt")):
        mnx_id, name, xref = line.rstrip("\n").split("\t")
        compartments[mnx_id] = Compartment(mnx_id, name, xref)
    logger.info(f"Loaded {len(compartments)} compartments")

    compartment_xrefs = 0
    for line in _iterate_tsv(gzip.open("data/comp_xref.tsv.gz", "rt")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            namespace, reference = xref.split(":", 1)
            compartment = compartments[mnx_id]
            compartment.annotation[namespace].append(reference)
            compartment_xrefs += 1
    logger.info(f"Loaded {compartment_xrefs} compartment cross-references")

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

    reaction_xrefs = 0
    reaction_xrefs_missing = 0
    for line in _iterate_tsv(gzip.open("data/reac_xref.tsv.gz", "rt")):
        xref, mnx_id, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            try:
                reaction = reactions[mnx_id]
            except KeyError:
                # At the time of this writing, 815 references don't exist in the
                # reaction database. Some are deprecated, but some aren't there
                # for unknown reasons.
                reaction_xrefs_missing += 1
            else:
                namespace, reference = xref.split(":", 1)
                reaction.annotation[namespace].append(reference)
                reaction_xrefs += 1
    logger.info(
        f"Loaded {reaction_xrefs} reaction cross-references (ignored "
        f"{reaction_xrefs_missing} unknown references)"
    )

    for line in _iterate_tsv(gzip.open("data/chem_prop.tsv.gz", "rt")):
        mnx_id, name, formula, _, _, _, _, _, _ = line.rstrip("\n").split("\t")
        metabolites[mnx_id] = Metabolite(mnx_id, name, formula)
    logger.info(f"Loaded {len(metabolites)} metabolites")

    metabolite_xrefs = 0
    metabolite_xrefs_missing = 0
    for line in _iterate_tsv(gzip.open("data/chem_xref.tsv.gz", "rt")):
        xref, mnx_id, _, _ = line.rstrip("\n").split("\t")
        if ":" in xref:
            try:
                metabolite = metabolites[mnx_id]
            except KeyError:
                # At the time of this writing, there is only one deprecated
                # reference which doesn't exist in the metabolite database.
                metabolite_xrefs_missing += 1
            else:
                namespace, reference = xref.split(":", 1)
                metabolite.annotation[namespace].append(reference)
                metabolite_xrefs += 1
    logger.info(
        f"Loaded {metabolite_xrefs} metabolite cross-references (ignored "
        f"{metabolite_xrefs_missing} unknown references)"
    )


def _iterate_tsv(file_):
    with file_:
        for line in file_:
            if line.startswith("#"):
                continue
            yield line
