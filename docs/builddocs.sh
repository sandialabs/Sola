#!/usr/bin/env bash

# builddocs.sh
# Build the sphinx documentation, enabling sphinxcontrib.matlab
# to document abstract MATLAB methods by temporarily rewriting
# them as non-abstract methods.

TEMPBRANCH="_temporarydocbuild"

# Stash any changes and switch to a temporary branch.
CURRENTBRANCH=`git rev-parse --abbrev-ref HEAD`
STASHED=`git stash push -m "temporary stash for documentation build"`
git switch --create ${TEMPBRANCH}

# Modify the MATLAB source code as needed.
python3 format_abstracts.py

# Build the documentation.
make html

# Switch back to the old git state.
git reset --hard
git switch ${CURRENTBRANCH}
git branch -D ${TEMPBRANCH}
if [[ ${STASHED} != "No local changes to save" ]]; then
    git stash pop
fi
