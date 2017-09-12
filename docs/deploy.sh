#!/bin/bash
# Update documentation on gh-pages branch after a manual build.

# Ensure we've committed all recent changes
# https://stackoverflow.com/a/3879077
require_clean_work_tree () {

    # Update the index
    git update-index -q --ignore-submodules --refresh
    err=0

    # Disallow unstaged changes in the working tree
    if ! git diff-files --quiet --ignore-submodules --
    then
        echo >&2 "Cannot $1: you have unstaged changes."
        git diff-files --name-status -r --ignore-submodules -- >&2
        err=1
    fi

    # Disallow uncommitted changes in the index
    if ! git diff-index --cached --quiet HEAD --ignore-submodules --
    then
        echo >&2 "Cannot $1: your index contains uncommitted changes."
        git diff-index --cached --name-status -r --ignore-submodules HEAD -- >&2
        err=1
    fi

    if [ $err = 1 ]
    then
        echo >&2 "Please commit or stash them."
        exit 1
    fi
}

# Exit on errors
set -o errexit -o nounset

# Ensure we're good to publish docs
require_clean_work_tree "publish docs"

# Begin
branch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
echo "Building docs from ${branch} branch..."

# Get git hash
rev=$(git rev-parse --short HEAD)

# Copy the html folder to a temporary location, initialize
# a new git repo there, add the necessary files, and force-push 
# to planetplanet/gh-pages
cd .build
cp -r html tmp_html
cp ../title.png tmp_html/_images/title.png
cp ../title.gif tmp_html/_images/title.gif
cp ../PPOs.pdf tmp_html/PPOs.pdf
cd tmp_html
git init
touch .nojekyll
git add -f .nojekyll
git add -f *.html
git add -f *.js
git add -f _sources
git add -f _static
git add -f _images
git add -f scripts/*.html
git add -f api/*.html
git add -f PPOs.pdf
git -c user.name='sphinx' -c user.email='sphinx' commit -m "rebuild gh-pages at ${rev}"
git push -f https://github.com/rodluger/planetplanet.git HEAD:gh-pages

# Remove the temporary directory
cd ..
rm -rf tmp_html