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

"""Functions to read MetaNetX source files."""

import gzip
import json
import logging

from .data import (
    Compartment,
    Metabolite,
    Reaction,
    compartments,
    metabolite_key_index,
    metabolites,
    reaction_key_index,
    reactions,
)


logger = logging.getLogger(__name__)


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
            namespace, reference = _miriam_identifiers(
                "compartment", namespace, reference
            )
            compartment = compartments[mnx_id]
            compartment.annotation[namespace].append(reference)
            compartment_xrefs += 1
    logger.info(f"Loaded {compartment_xrefs} compartment cross-references")

    filtered_reaction_count = 0
    for line in _iterate_tsv(gzip.open("data/reac_prop.tsv.gz", "rt")):
        mnx_id, equation, _, _, ec, _ = line.rstrip("\n").split("\t")
        name = reaction_names.get(mnx_id)
        try:
            reaction = Reaction(mnx_id, name, equation, ec)
            reactions[mnx_id] = reaction
            reaction_key_index[mnx_id.lower()] = reaction
            # We haven't been able to map names for all reactions, so ignore the
            # missing ones.
            if name:
                reaction_key_index[name.lower()] = reaction
            reaction_key_index[ec.lower()] = reaction
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
                namespace, reference = _miriam_identifiers(
                    "reaction", namespace, reference
                )
                reaction.annotation[namespace].append(reference)
                reaction_key_index[reference.lower()] = reaction
                reaction_xrefs += 1
    logger.info(
        f"Loaded {reaction_xrefs} reaction cross-references (ignored "
        f"{reaction_xrefs_missing} unknown references)"
    )

    for line in _iterate_tsv(gzip.open("data/chem_prop.tsv.gz", "rt")):
        mnx_id, name, formula, _, _, _, _, _, _ = line.rstrip("\n").split("\t")
        metabolite = Metabolite(mnx_id, name, formula)
        metabolites[mnx_id] = metabolite
        metabolite_key_index[mnx_id.lower()] = metabolite
        metabolite_key_index[name.lower()] = metabolite
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
                namespace, reference = _miriam_identifiers(
                    "metabolite", namespace, reference
                )
                metabolite.annotation[namespace].append(reference)
                metabolite_key_index[reference.lower()] = metabolite
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


def _miriam_identifiers(type_, namespace, identifier):
    """
    Fix the MetaNetX identifiers into miriam equivalents.

    MetaNetX doesn't use correct miriam identifiers. This function maps the
    known namespace and entity identifiers used by MetaNetX to valid miriam
    identifiers.

    Parameters
    ----------
    type_ : string
        "compartment", "reaction" or "metabolite"
    namespace : string
        The MetaNetX namespace identifier
    identifier : string
        The object identifier

    Returns
    -------
    namespace : string
        The corrected namespace
    identifier : string
        The corrected identifier

    """
    if type_ == "compartment":
        ns_map = {
            "bigg": "bigg.compartment",
            "cco": "cco",
            "go": "go",
            "name": "name",  # unconfirmed
            "seed": "seed",
        }
        return (ns_map[namespace], identifier)
    elif type_ == "reaction":
        ns_map = {
            "bigg": "bigg.reaction",
            "deprecated": "metanetx.reaction",
            "kegg": "kegg.reaction",
            "metacyc": "metacyc.reaction",
            "reactome": "reactome",
            "rhea": "rhea",
            "sabiork": "sabiork.reaction",
            "seed": "seed.reaction",
        }
        return (ns_map[namespace], identifier)
    elif type_ == "metabolite":
        if namespace == "kegg":
            kegg_map = {
                "C": "kegg.compound",
                "D": "kegg.drug",
                "E": "kegg.environ",
                "G": "kegg.glycan",
            }
            return (kegg_map[identifier[0]], identifier)
        elif namespace == "slm":
            return ("swisslipid", f"SLM:{identifier}")
        elif namespace == "chebi":
            return (namespace, f"CHEBI:{identifier}")
        else:
            ns_map = {
                "bigg": "bigg.metabolite",
                "deprecated": "metanetx.chemical",
                "envipath": "envipath",  # unconfirmed
                "hmdb": "hmdb",
                "lipidmaps": "lipidmaps",
                "metacyc": "metacyc.compound",
                "reactome": "reactome",
                "sabiork": "sabiork.compound",
                "seed": "seed.compound",
            }
            return (ns_map[namespace], identifier)
