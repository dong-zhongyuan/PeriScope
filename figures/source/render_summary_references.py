"""Render the user's supplied wheel structure with labelled illustrative data only.
The actual Figure 7D reads evidence_wheel.csv and never reads this template data.
"""
from pathlib import Path
import sys,types,runpy
import numpy as np,pandas as pd,matplotlib.pyplot as plt,fitz
W=Path(__file__).resolve().parents[1];REF=W/'references'
targets=['CD62P','CD22','CD35','CD20','CD36','CD226','CD45RB']
layers=[f'Evidence {i}'for i in range(1,10)]
data=pd.DataFrame([[int((i*3+j*2)%7>1)for j in range(9)]for i in range(7)],index=targets,columns=layers)
data.to_csv(REF/'illustrative_template_data_ONLY.csv')
S=types.ModuleType('style');S.LAYER_NAMES=layers;S.LAYER_COLORS=['#599CAF','#D0AE67','#C66C7E','#879C66','#9D88AE','#6D91AF','#4A9D8B','#9C887B','#8F95B9'];S.FS_SMALL=8;S.C_MICRO='#599783';S.mm_w=lambda x:7.5;S.apply=lambda:plt.rcParams.update({'font.family':'Arial','savefig.bbox':None});S.dump_csv=lambda *a:None
def save(fig,name):
    fig.text(.5,.99,'ILLUSTRATIVE TEMPLATE DATA',va='top',ha='center',size=12,color='#8A5B39')
    fig.savefig(REF/'user_wheel_structure.pdf',bbox_inches='tight')
S.save=save;sys.modules['style']=S
C=types.ModuleType('fig5_common');C.TARGETS=targets;C.build=lambda:data;sys.modules['fig5_common']=C
runpy.run_path(str(REF/'user_fig5_e_wheel.py'),run_name='__main__')
with fitz.open(REF/'user_wheel_structure.pdf') as d:d[0].get_pixmap(matrix=fitz.Matrix(1.6,1.6)).save(REF/'user_wheel_structure.png')
