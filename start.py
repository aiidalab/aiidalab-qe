import ipywidgets as ipw
from IPython.core.display import HTML
import urllib2
import os
#HTML(open("/project/apps/calcexamples/custom1.css").read())
#HTML(os.getcwd())

def get_start_widget(appbase, jupbase, notebase):
    #http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="width:200px" rowspan=2></th>
        <th style="width:300px" rowspan=2></th>
        <th style="width:150px" rowspan=2></th>
        <th style="width:150px" rowspan=2></th>
        <th style="width:150px" rowspan=2></th>
        <th style="width:200px" rowspan=2></th>
    </tr>

    <tr>
    <tr>
    <td valign="center" rowspan="2">
    <ul style="list-style-type:none">
        <li><a href="{appbase}/structures.ipynb" target="_blank">
        <i class="fa fa-upload" style="color:#337ab7;font-size:8em;"></i></a></li>
        <li> Upload structure </li>
    </ul>
    </td>
    </td>


    <td valign="center">
    <ul style="list-style-type:none">
        <li> <img src="https://gitlab.com/QEF/q-e/raw/develop/logo.jpg" width="100%" height="100%"> </li>
    </ul>
    </td>

    <td valign="center">
    <ul style="list-style-type:none">
        <li> <i class="fa fa-arrow-right fa-4x" style="color:#337ab7;"></i> </li>
    </ul>
    </td>

    <td valign="center">
    <ul style="list-style-type:none">
        <li><big><center>Geometry optimization</center></big></li>
        <li><a href="{appbase}/qe/cell_opt.ipynb" target="_blank">
        <i class="fa fa-desktop" style="color:#337ab7;font-size:8em;" ></i> </a></li>
    </ul>
    </td>
    
    <td valign="center">
    <ul style="list-style-type:none">
        <i class="fa fa-arrow-right fa-4x" style="color:#337ab7;"></i>
    </ul>
    </td>

    <td valign="center">
    <ul style="list-style-type:none">
        <li><center><big>Equation of States</big></center></li>
        <a href="{appbase}/qe/equation_of_states.ipynb" target="_blank">
        <li><i class="fa fa-refresh fa-spin-hover" style="color:#337ab7;font-size:8em;" ></i></li></a>
    </ul>
    </td>
    

    </tr>
    
    <tr>
    <td valign="center">
    <ul style="list-style-type:none">
        <li> <img src="https://www.cp2k.org/_media/cp2k_logo_300.png" width="100%" height="100%"> </li>
    </ul>
    </td>
    
    <td valign="center">
    <ul style="list-style-type:none">
        <li> <i class="fa fa-arrow-right fa-4x" style="color:#337ab7;"></i> </li>
    </ul>
    </td>

    <td valign="center">
    <ul style="list-style-type:none">
        <li><a href="{appbase}/cp2k/cell_opt.ipynb" target="_blank">
        <i class="fa fa-desktop" style="color:#337ab7;font-size:8em;" ></i></a></li>
    </ul>
    </td>
    
    <td valign="center">
    <ul style="list-style-type:none">
        <i class="fa fa-arrow-right fa-4x" style="color:#337ab7;"></i>
    </ul>
    </td>

    <td valign="center">
    <ul style="list-style-type:none">
        <a href="{appbase}/cp2k/equation_of_states.ipynb" target="_blank">
        <li><i class="fa fa-refresh fa-spin-hover" style="color:#337ab7;font-size:8em;" ></i></li></a>
    </ul>
    </td>
    

    </tr>
    
    </table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase, notebase=notebase)
    return ipw.HTML(html)
    
#EOF
