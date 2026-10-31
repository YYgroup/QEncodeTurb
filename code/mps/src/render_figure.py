"""Execute the portable plotting script in the generated figure directory."""
import json
import os
import runpy
import shutil
from paths import ROOT, PACKAGE

def main():
    folder = ROOT/'figure'
    folder.mkdir(exist_ok=True)
    for name in ('plot_mps_sixpanel.py','plot_mps_sixpanel.ipynb'):
        shutil.copy2(PACKAGE/'templates'/name,folder/name)
    os.chdir(folder)
    runpy.run_path(str(folder/'plot_mps_sixpanel.py'),run_name='__main__')
    print(json.dumps(dict(pdf=str(folder/'mps_sixpanel_17p9cm.pdf'),
                         notebook=str(folder/'plot_mps_sixpanel.ipynb')),indent=2))

if __name__ == '__main__':
    main()
