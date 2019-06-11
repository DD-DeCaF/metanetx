# MetaNetX

[![Build Status](https://travis-ci.org/DD-DeCaF/metanetx.svg?branch=master)](https://travis-ci.org/DD-DeCaF/metanetx)
[![Codecov](https://codecov.io/gh/DD-DeCaF/metanetx/branch/master/graph/badge.svg)](https://codecov.io/gh/DD-DeCaF/metanetx/branch/master)

## Development

Run `make setup` first when initializing the project for the first time. Type
`make` to see all commands.

### Source files

The MetaNetX source files are stored in DD-DeCaF's google cloud storage
retrieved on startup. They're not stored in the repo as they're too large, and
we're not downloading them directly from metanetx.org because 1) our deployments
should not depend on the availability on the metanetx.org service, and 2) to
avoid compatibility issues if the source file format is updated.

To speed up source file reads during development, download them into the `data/`
folder and set `LOCAL_METANETX_DATA=1` in your `.env`.

Reaction names are not part of MetaNetX, but collected manually by running
./scripts/generate_reaction_names.py`. Note that the script takes several hours
to complete. Names are retrieved from cross referenced databases (currently
BiGG, kegg, ModelSEED and EC numbers are checked).

### Environment

Specify environment variables in a `.env` file. See `docker-compose.yml` for the
possible variables and their default values.

* Set `ENVIRONMENT` to either
  * `development`,
  * `testing`, or
  * `production`.
* `SECRET_KEY` Flask secret key. Will be randomly generated in development and testing environments.
* `SENTRY_DSN` DSN for reporting exceptions to
  [Sentry](https://docs.sentry.io/clients/python/integrations/flask/).
* `ALLOWED_ORIGINS`: Comma-seperated list of CORS allowed origins.

### Code style

In order of priority, code must adhere to the rules of the following tools:

1. [black](https://github.com/ambv/black)
2. [flake8](http://flake8.pycqa.org/en/latest/)
    * pycodestyle
    * pyflakes
    * mccabe
    * [pydocstyle](http://www.pydocstyle.org/en/2.1.1/index.html)
    * [bugbear](https://github.com/PyCQA/flake8-bugbear)
3. The [NumPy docstring standard](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard)
4. [isort](https://github.com/timothycrosley/isort)
