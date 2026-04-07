unittest unittests test-unittest test-unittests:
	pytest --verbose --tb=short autotst/
	rm -rf species ts test/species test/ts
