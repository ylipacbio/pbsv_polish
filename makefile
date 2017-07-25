# pytest can handle these nose flags. We will have to see about others.
# e.g. --junit-xml=nosetests.xml
MY_NOSE_FLAGS?=-v -s
MY_CRAM_FLAGS?=-v
THISDIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))

test-fast:
	# These use the local tree. No installation required. Extremely fast.
	py.test ${MY_NOSE_FLAGS} utests.py
doctest:
	py.test --doctest-modules pbsv/
pylint:
	pylint --errors-only pbsv/
wheel:
	# This is basically a syntax check.
	python setup.py bdist_wheel
indep:
	py.test ${MY_NOSE_FLAGS} tests/itest/test_indep.py
test:
	# Relatively fast tests
	py.test ${MY_NOSE_FLAGS} tests/itest/test_*.py
bamboo:
	bash bamboo_pbsv_build.sh
	bamboo_build_working_directory=$(shell pwd) bash bamboo_pbsv_test.sh
build:
	pip install --user --edit .
coverage:
	${MAKE} coverage-install
	# Test sitecustomize.
	python -c 'import sitecustomize as s; print s'
	${MAKE} coverage-clean
	${MAKE} coverage-test
	${MAKE} coverage-report
coverage-%:
	# Test sitecustomize.
	python -c 'import sitecustomize as s; print s'
	# make target, with coverage.
	COVERAGE_PROCESS_START=${THISDIR}/mycoverage.cfg ${MAKE} $*
coverage-report:
	ls -larth
	# We assume the .coverage* files are in ${PYTHONUSERBASE}/, from mycoverage.cfg
	coverage combine ${PYTHONUSERBASE}/
	ls -larth
	coverage xml -o coverage.xml
	sed -i -e 's@filename="@filename="./@g' coverage.xml
	coverage report -m
coverage-clean:
	rm -f ${PYTHONUSERBASE}/.coverage* .coverage* coverage.xml
coverage-install:
	# This is needed only if you run from a different directory, since ./sitecustomize.py
	# would not be in 'sys.path'.
	# Assume PYTHONUSERBASE is set.
	#pip install --user coverage
	mkdir -p ${PYTHONUSERBASE}/lib/python2.7/site-packages
	ln -sf ${THISDIR}/mysitecustomize.py ${PYTHONUSERBASE}/lib/python2.7/site-packages/sitecustomize.py
coverage-uninstall:
	rm -f ${PYTHONUSERBASE}/lib/python2.7/site-packages/sitecustomize.py*
