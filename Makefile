.PHONY: all test clean python-test

all:
	$(MAKE) -C cpp all

test:
	$(MAKE) -C cpp test

clean:
	$(MAKE) -C cpp clean

python-test:
	cd python && python3 -m pip install -e ".[dev]" && python3 -m pytest
