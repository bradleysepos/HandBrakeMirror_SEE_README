#! /bin/bash
#
# Retrieves git repository info for directory ${1}

REPO_DIR=${1}

if [ "x${REPO_DIR}" != "x" ]; then
    cd ${REPO_DIR}
fi

# check if there is a valid git repo here
HASH=$(git rev-parse --short HEAD)
ERR=$?
if [ ${ERR} -ne 0 ]; then
    # Proably not a valid repo, bail
    exit ${ERR}
fi

# Only write tag and rev if they exist. A fresh clone currently has no tags.
URL=$(git config remote.origin.url)
TAG=$(git describe --abbrev=0)
ERR=$?
if [ ${ERR} -eq 0 ]; then
    echo "TAG=${TAG}"
    REV=$(git rev-list ${TAG}.. --count)
    ERR=$?
    if [ ${ERR} -eq 0 ]; then
        echo "REV=${REV}"
    fi
fi

DATE=$(git log -1 --format="format:%ai")
BRANCH=$(git symbolic-ref -q --short HEAD)
UPSTREAM=$(git config branch.${BRANCH}.remote)
if [ "x${UPSTREAM}" != "x" ]; then
    REMOTE=$(git config remote.${UPSTREAM}.url)
else
    REMOTE=${URL}
fi

echo "URL=${URL}"
echo "HASH=${HASH}"
echo "DATE=${DATE}"
echo "BRANCH=${BRANCH}"
echo "REMOTE=${REMOTE}"
