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

"""Implement RESTful API endpoints using resources."""

from flask_apispec import MethodResource, marshal_with, use_kwargs
from flask_apispec.extension import FlaskApiSpec

from . import data
from .schemas import ReactionSearchSchema, ReactionResponseSchema


def init_app(app):
    """Register API resources on the provided Flask application."""

    def register(path, resource):
        app.add_url_rule(path, view_func=resource.as_view(resource.__name__))
        docs.register(resource, endpoint=resource.__name__)

    docs = FlaskApiSpec(app)
    app.add_url_rule("/healthz", view_func=healthz)
    register("/reactions", ReactionResource)


def healthz():
    """
    Return an empty, successful response for readiness checks.

    A successful response signals that the app is initialized and ready to
    receive traffic. The main use case is for apps with slow initialization, but
    external dependencies like database connections can also be tested here.
    """
    return ""


class ReactionResource(MethodResource):
    @use_kwargs(ReactionSearchSchema)
    @marshal_with(ReactionResponseSchema(many=True), code=200)
    def get(self, query):
        # Search through the data store for matching reactions.
        reactions = [r for r in data.reactions.values() if r.mnx_id == query]
        # Collect all unique references to metabolites and compartments, and
        # include the objects in the response.
        results = []
        for reaction in reactions:
            metabolite_ids = set(m["metabolite_id"] for m in reaction.equation_parsed)
            compartment_ids = set(m["compartment_id"] for m in reaction.equation_parsed)
            results.append(
                {
                    "reaction": reaction,
                    "metabolites": [
                        data.metabolites[m] for m in metabolite_ids
                    ],
                    "compartments": [
                        data.compartments[m] for m in compartment_ids
                    ],
                }
            )
        return results
