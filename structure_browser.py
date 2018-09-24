from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.structure import StructureData
from aiida.orm.data.cif import CifData
from aiida.orm.calculation.job import JobCalculation
from aiida.orm.calculation.work import WorkCalculation
from aiida.orm import Node
    
from collections import OrderedDict
import ipywidgets as ipw
import datetime

class StructureBrowser(ipw.VBox):
    
    def __init__(self, mymode='uploaded'):
        # Find all process labels
        qb = QueryBuilder()
        qb.append(WorkCalculation,
                  project="attributes._process_label",
                  filters={
                      'attributes': {'!has_key': 'source_code'}
                  }
        )
        qb.order_by({WorkCalculation:{'ctime':'desc'}})

        process_labels = []
        for i in qb.iterall():
            if i[0] not in process_labels:
                process_labels.append(i[0])

        layout = ipw.Layout(width="900px")

        self.mode = mymode
                
        hr = ipw.HTML('<hr>')
        
        self.results = ipw.Dropdown(layout=layout)
        self.search()
        super(StructureBrowser, self).__init__([hr, self.results])

    
    def preprocess(self):
        qb = QueryBuilder()
        filters = {'extras': {'!has_key': 'formula'}}
        filters['or'] = [{'type':CifData._plugin_type_string},{'type':StructureData._plugin_type_string}]
        qb.append(Node, filters=filters)
        for n in qb.all(): # iterall() would interfere with set_extra()
            try:
                formula = n[0].get_formula()
            except:
                formula = n[0].get_formulae()[0]
            n[0].set_extra("formula", formula)

    
    def search(self, c=None):
        self.preprocess()
        
        qb = QueryBuilder()


        filters = {}
        filters['or'] = [{'type':CifData._plugin_type_string},{'type':StructureData._plugin_type_string}]
        if self.mode== "uploaded":
            qb2 = QueryBuilder()
            qb2.append(StructureData, project=["id"])
            qb2.append(Node, input_of=StructureData)
            processed_nodes = [n[0] for n in qb2.all()]
            if processed_nodes:
                filters['id'] = {"!in":processed_nodes}
            qb.append(StructureData, filters=filters)

        elif self.mode == "calculated":
            qb.append(JobCalculation)
            qb.append(Node, output_of=JobCalculation, filters=filters)
        else:
            pass

        matches = set([n[0] for n in qb.iterall()])
        matches = sorted(matches, reverse=True, key=lambda n: n.ctime)
        
        c = len(matches)
        options = OrderedDict()
        for n in matches:
            label  = "PK: %d" % n.pk
            label += " | " + n.ctime.strftime("%Y-%m-%d %H:%M")
            label += " | " + n.get_extra("formula")
            label += " | " + n.description
            options[label] = n

        self.results.options = options
