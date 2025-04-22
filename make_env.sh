set -e # force failure if anyting fails in script - ensuring completion
set -o errexit
ENVNAME="dandi_venv"
if [ -d $ENVNAME ]
then
    echo "Removing previous environment: $ENVNAME"
    set +e
    rsync -aR --remove-source-files $ENVNAME ~/.Trash/ || exit 1
    set -e
    rm -R $ENVNAME
else
    echo "No previous environment - ok to proceed"
fi
PYTHONVER="python3.13"
$PYTHONVER  -m venv $ENVNAME
source $ENVNAME/bin/activate

pip3 install --upgrade pip  # be sure pip is up to date in the new env.
python -m ensurepip --upgrade
pip3 install wheel  # seems to be missing (note singular)
pip3 install setuptools
pip3 install -r requirements.txt
source $ENVNAME/bin/activate
python3 --version

