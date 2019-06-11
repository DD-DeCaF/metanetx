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


class ReactionSearchSchema(Schema):
    query = fields.Str(required=True)

    class Meta:
        strict = True


class CompartmentSchema(Schema):
    mnx_id = fields.Str()
    name = fields.Str()
    xref = fields.Str()


class ReactionSchema(Schema):
    mnx_id = fields.Str()
    equation = fields.Str()
    ec = fields.Str()


class MetaboliteSchema(Schema):
    mnx_id = fields.Str()
    description = fields.Str()


class ReactionResponseSchema(Schema):
    reaction = fields.Nested(ReactionSchema)
    metabolites = fields.Nested(MetaboliteSchema, many=True)
    compartments = fields.Nested(CompartmentSchema, many=True)
