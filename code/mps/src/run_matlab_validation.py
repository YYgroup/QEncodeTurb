"""Optionally rerun the untouched upstream routines in MATLAB or GNU Octave."""
import argparse
import os
from pathlib import Path
import subprocess
import shutil
import sys
from paths import ROOT, PACKAGE

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--octave',help='Path or command for octave-cli.')
    group.add_argument('--matlab',help='Path or command for MATLAB.')
    args=parser.parse_args()
    for relative in ('validation/original_inputs.mat','recheck/joint_original_inputs.mat'):
        target=ROOT/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():
            shutil.copy2(PACKAGE/'reference'/relative,target)
    folder=(PACKAGE/'matlab').as_posix().replace("'","''")
    expression=f"addpath('{folder}'); validate_original; validate_joint_original;"
    env=os.environ.copy()
    env['MPS_WORKDIR']=str(ROOT)
    if args.octave:
        command=[args.octave,'--quiet','--no-gui','--eval',expression]
    else:
        command=[args.matlab,'-batch',expression]
    subprocess.run(command,check=True,env=env,cwd=PACKAGE)
    subprocess.run([sys.executable,str(PACKAGE/'src/validate_original.py')],check=True,env=env)
    subprocess.run([sys.executable,str(PACKAGE/'src/validate_joint_original.py')],check=True,env=env)

if __name__=='__main__':
    main()
