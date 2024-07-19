#!/bin/bash

# This script is configured to run automatically by bumpver
# before it creates the release commit.
# We check for common mistakes, such as making a release commit
# in a wrong branch, or trying to push to a wrong remote.
#
# For now, only two checks are implemented:
#
# 1. Check that the current branch matches either release/* or support/*
#
# 2. Check that the remote 'origin' is pointing to the origin repository,
#    and not a fork. Note however that this assumes that origin is the default remote
#    where new branches are pushed. If the user configured a different default remote,
#    this check will not save them in the current implementation.
#
# Future work:
#  - make sure the main branch is up-to-date with origin/main
#  - make sure the HEAD commit was branched off of main branch,
#    although this rule should only apply to release/* branches, not support/* branches
#
# Ideally, some of these check would be handled by bumpver itself:
# Restricting releases from branch: https://github.com/mbarkhau/bumpver/issues/198
# Restricting releases to specified remote: https://github.com/mbarkhau/bumpver/issues/234

set -euo pipefail

ORIGIN="github\.com[:/]aiidalab/aiidalab-qe"

error=0

branch=$(git branch --show-current)

# Explicitly disallow master/main branch
if [[ $branch = "master" || $branch = "main"  ]];then
  echo "ERROR: You should not run bumpver from main/master branch!"
  echo "Make sure your main branch is up-to-date with origin ('git pull origin main')"
  echo "and create a release branch first, e.g. 'git switch -c release/v2.0.0'"
  error=1
fi

# Only allow release/* and support/* branches
if [[ ! $branch =~ 'release/' && ! $branch =~ 'support/' ]];then
  echo "ERROR: The branch name must be either release/<version> or support/<version>"
  error=1
fi

# TODO: We need to check which remote is actually configured for push!
origin_url=$(git remote get-url --push --all origin)
if [[ ! $origin_url =~ $ORIGIN ]];then
  echo "ERROR: Wrong default repo remote set!"
  echo "got: $origin_url"
  echo "expected: $ORIGIN"
  error=1
fi

if [[ $error != 0 ]];then
  exit 1
fi
