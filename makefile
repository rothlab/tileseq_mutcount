TOKEN:=~/.twine/token_tileseqmut

build:
	rm -rf dist build
	python3.8 setup.py bdist_wheel 

install: build
	pip install --user dist/TileSeqMut-*.whl

upload: build
	twine upload -u __token__ -p $$(cat $(TOKEN)) dist/TileSeqMut-*.whl

