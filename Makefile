unittest unittests test-unittest test-unittests:
	nosetests --nocapture --nologcapture --all-modules --verbose --with-coverage  --cover-package=autotst 
	rm -r species
	rm -r ts