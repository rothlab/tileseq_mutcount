TOKEN:=~/.twine/token_tileseqmut

wheel:
	rm -rf dist wheel
	python3.8 setup.py bdist_wheel 

install: wheel
	pip install --user dist/TileSeqMut-*.whl

upload: wheel
	twine upload -u __token__ -p $$(cat $(TOKEN)) dist/TileSeqMut-*.whl

