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

import warnings

from flask_apispec import MethodResource, marshal_with, use_kwargs
from flask_apispec.extension import FlaskApiSpec

from . import data
from .schemas import (
    BatchSearchSchema,
    MetaboliteSchema,
    ReactionResponseSchema,
    SearchSchema,
)


def init_app(app):
    """Register API resources on the provided Flask application."""

    def register(path, resource):
        app.add_url_rule(path, view_func=resource.as_view(resource.__name__))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            docs.register(resource, endpoint=resource.__name__)

    docs = FlaskApiSpec(app)
    app.add_url_rule("/healthz", view_func=healthz)
    register("/reactions", ReactionResource)
    register("/reactions/batch", ReactionBatchResource)
    register("/metabolites", MetaboliteResource)
    register("/metabolites/batch", MetaboliteBatchResource)


def healthz():
    """
    Return an empty, successful response for readiness checks.

    A successful response signals that the app is initialized and ready to
    receive traffic. The main use case is for apps with slow initialization, but
    external dependencies like database connections can also be tested here.
    """
    return ""


class ReactionResource(MethodResource):
    @use_kwargs(SearchSchema)
    @marshal_with(ReactionResponseSchema(many=True), code=200)
    def get(self, query):
        # Search through the data store for matching reactions.
        reactions = sorted(
            data.reactions.values(), key=lambda r: r.match(query), reverse=True
        )

        # Limit the results to the first 30.
        reactions = reactions[:30]

        # Collect all unique references to metabolites and compartments, and
        # include the objects in the response.
        return [reaction.with_references() for reaction in reactions]


class ReactionBatchResource(MethodResource):
    @use_kwargs(BatchSearchSchema)
    @marshal_with(ReactionResponseSchema(many=True), code=200)
    def get(self, query):
        # Search through the data store for multiple exact matching reactions.
        results = []
        for q in query:
            try:
                results.append(
                    data.reaction_key_index[q.lower()].with_references()
                )
            except KeyError:
                results.append(None)
        return results


class MetaboliteResource(MethodResource):
    @use_kwargs(SearchSchema)
    @marshal_with(MetaboliteSchema(many=True), code=200)
    def get(self, query):
        # Search through the data store for matching metabolites.
        metabolites = [m for m in data.metabolites.values() if m.match(query)]

        # Limit the results to the first 30.
        metabolites = metabolites[:30]
        return metabolites


class MetaboliteBatchResource(MethodResource):
    @use_kwargs(BatchSearchSchema)
    @marshal_with(MetaboliteSchema(many=True), code=200)
    def get(self, query):
        # Search through the data store for multiple exact matching reactions.
        return [data.metabolite_key_index.get(q.lower()) for q in query]
