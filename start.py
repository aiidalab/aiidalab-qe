import ipywidgets as ipw
from IPython.core.display import HTML
import os
#HTML(open("/project/apps/calcexamples/custom1.css").read())
#HTML(os.getcwd())

def get_start_widget(appbase, jupbase, notebase):
    #http://fontawesome.io/icons/
    template = """
    <table style="width:100%">
      <tr>
        <th></th>
        <th><h4><b>Please choose the code:</h4></b></th> 
        <th></th>
      </tr>
      <tr>
        <th> <div align="center"> <a href="{appbase}/cp2k.ipynb" target="_blank"> <img src="https://www.cp2k.org/_media/cp2k_logo_300.png" height="120px" width="120px"> </a> </div> </th>
        <th> <div align="center"><font size="+2">OR</font></div></th> 
        <th><div align="center"> <a href="{appbase}/qe.ipynb" target="_blank"> <img src="https://gitlab.com/QEF/q-e/raw/develop/logo.jpg" height="120px" width="243px"> </a> </div> </th>
      </tr>
    </table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase, notebase=notebase)
    return ipw.HTML(html)
    
#EOF
