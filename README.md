# MetaNetX

![master Branch](https://img.shields.io/badge/branch-master-blue.svg)
[![master Build Status](https://travis-ci.org/DD-DeCaF/metanetx.svg?branch=master)](https://travis-ci.org/DD-DeCaF/metanetx)
[![master Codecov](https://codecov.io/gh/DD-DeCaF/metanetx/branch/master/graph/badge.svg)](https://codecov.io/gh/DD-DeCaF/metanetx/branch/master)

![devel Branch](https://img.shields.io/badge/branch-devel-blue.svg)
[![devel Build Status](https://travis-ci.org/DD-DeCaF/metanetx.svg?branch=devel)](https://travis-ci.org/DD-DeCaF/metanetx)
[![devel Codecov](https://codecov.io/gh/DD-DeCaF/metanetx/branch/devel/graph/badge.svg)](https://codecov.io/gh/DD-DeCaF/metanetx/branch/devel)

## Development

Run `make setup` first when initializing the project for the first time. Type
`make` to see all commands.

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
