#!/bin/bash
if [ "$TRAVIS_OS_NAME" == 'osx' ]; then
  export PATH="$HOME/miniconda/bin:$PATH"
  source $HOME/miniconda/bin/activate
  conda activate travis
  PYTHON=$(which python)
  PYTEST=$(which pytest)
else
  PYTHON=${PYTHON:-python}
  PYTEST=${PYTEST:-"pytest"}
fi

echo "python: ${PYTHON}"

echo 'Running tests...'

# LOCAL TESTING:
# PYTHONPATH=$(pwd)/build/lib:PYTHONPATH pytest -rxXs test

${PYTEST} "test"
