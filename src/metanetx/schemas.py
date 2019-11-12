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

"""Marshmallow schemas for marshalling the API endpoints."""

from marshmallow import Schema, fields
from webargs.fields import DelimitedList


class SearchSchema(Schema):
    query = fields.Str(required=True)


class BatchSearchSchema(Schema):
    query = DelimitedList(fields.Str(), required=True)


class CompartmentSchema(Schema):
    mnx_id = fields.Str()
    name = fields.Str()
    xref = fields.Str()
    annotation = fields.Raw()


class ReactionSchema(Schema):
    class EquationSchema(Schema):
        metabolite_id = fields.Str()
        compartment_id = fields.Str()
        coefficient = fields.Float()

    mnx_id = fields.Str()
    name = fields.Str()
    equation_string = fields.Str()
    equation_parsed = fields.Nested(EquationSchema, many=True)
    ec = fields.Str()
    annotation = fields.Raw()


class MetaboliteSchema(Schema):
    mnx_id = fields.Str()
    name = fields.Str()
    formula = fields.Str()
    annotation = fields.Raw()


class ReactionResponseSchema(Schema):
    reaction = fields.Nested(ReactionSchema)
    metabolites = fields.Nested(MetaboliteSchema, many=True)
    compartments = fields.Nested(CompartmentSchema, many=True)
