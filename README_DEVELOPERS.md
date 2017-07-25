## Quick start

    make

That tests only the code which has no external dependencies. It relies
on `.` implicitly in `$PYTHONPATH`, so no installation is required.

## Bamboo
To do what Bamboo does, try

    make bamboo

That should take about a minute.

## Environment
To get a development environment, try

    source bamboo_setup.sh

That will use our Anaconda Python and some GNU modules,
plus some pbcore Python.

Then,

    make
    make itest
    make coverage
    make ctest

Getting coverage on the cram tests is possible, but not trivial. Contact @cdunn.
