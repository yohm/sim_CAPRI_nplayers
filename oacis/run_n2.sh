#!/bin/bash -eux

script_dir=$(cd $(dirname $BASH_SOURCE); pwd)
$script_dir/../cpp/main_evo_n2.out $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}

export PIPENV_PIPFILE=$script_dir/Pipfile
pipenv run python $script_dir/plot_abundance.py
pipenv run python $script_dir/plot_cooperation_level_n2.py
