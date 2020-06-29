#!/bin/bash -eux

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/../cpp/main_evo $@

export PIPENV_PIPFILE=$script_dir/Pipfile
pipenv run python $script_dir/plot_abundance.py
pipenv run python $script_dir/plot_cooperation_level.py
