# pytest can handle these nose flags. We will have to see about others.
# e.g. --junit-xml=nosetests.xml
MY_NOSE_FLAGS?=-v -s
MY_CRAM_FLAGS?=-v
THISDIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))

test-fast:
	# These use the local tree. No installation required. Extremely fast.
	py.test ${MY_NOSE_FLAGS} utests.py
doctest:
	py.test --doctest-modules pbsv_polish/
pylint:
	pylint --errors-only pbsv_polish/
wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
indep:
	py.test ${MY_NOSE_FLAGS} tests/itest/test_indep.py
test:
	# Relatively fast tests
	py.test ${MY_NOSE_FLAGS} tests/itest/test_*.py
autofmt:
	find pbsv_polish -type f -name '*.py' | xargs autoflake --in-place --remove-unused-variables --expand-star-imports
	find pbsv_polish -type f -name '*.py' | xargs autopep8  --in-place --max-line-length 120
vulture:
	vulture pbsv_polish/
clean:
	find . -type f -name '*.pyc' |xargs rm -f
