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

    def with_references(self):
        """
        Return an object with this reaction and referenced objects.

        Returns a dictionary with keys:
            "reaction": This reaction object
            "metabolites": List of referred metabolite objects
            "compartments": List of referred compartments
        """
        metabolite_ids = set(m["metabolite_id"] for m in self.equation_parsed)
        compartment_ids = set(m["compartment_id"] for m in self.equation_parsed)
        return {
            "reaction": self,
            "metabolites": [metabolites[m] for m in metabolite_ids],
            "compartments": [compartments[m] for m in compartment_ids],
        }

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

# These dictionaries also include names and identifiers from other namespaces,
# lowercased for quick case-insensitive lookup.
reaction_key_index = {}
metabolite_key_index = {}
