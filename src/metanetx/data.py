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

from collections import namedtuple
import logging

logger = logging.getLogger(__name__)

Compartment = namedtuple("Compartment", ["mnx_id", "name", "xref"])
Reaction = namedtuple("Reaction", ["mnx_id", "equation", "ec"])
Metabolite = namedtuple("Metabolite", ["mnx_id", "description"])

# MetaNetX data will be read into memory into the following dicts, keyed by ID.
compartments = {}
reactions = {}
metabolites = {}


def read_metanetx_files():
    with open("data/comp_prop.tsv") as file_:
        for line in file_:
            if line.startswith("#"):
                continue
            mnx_id, name, xref = line.rstrip("\n").split("\t")
            compartments[mnx_id] = Compartment(
                mnx_id=mnx_id, name=name, xref=xref
            )
    logger.info(f"Loaded {len(compartments)} compartments")

    with open("data/reac_prop.tsv") as file_:
        for line in file_:
            if line.startswith("#"):
                continue
            mnx_id, equation, _, _, ec, _ = line.rstrip("\n").split("\t")
            reactions[mnx_id] = Reaction(
                mnx_id=mnx_id, equation=equation, ec=ec
            )
    logger.info(f"Loaded {len(reactions)} reactions")

    with open("data/chem_prop.tsv", "rt") as file_:
        for line in file_:
            if line.startswith("#"):
                continue
            mnx_id, description, _, _, _, _, _, _, _ = line.rstrip("\n").split(
                "\t"
            )
            metabolites[mnx_id] = Metabolite(
                mnx_id=mnx_id, description=description
            )
    logger.info(f"Loaded {len(metabolites)} metabolites")
