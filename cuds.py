import json
from urllib.parse import urlsplit, parse_qs

import ase
import requests
import traitlets
import ipywidgets as ipw

from osp.core.utils import deserialize
from marketplace import MarketPlace


OPTIMADE_GATEWAY_APP='a99f9f38-90b0-497b-8dcc-dc2d7926d534'


class CUDSImporter(ipw.HBox):
    
    title = 'CUDS'
    structure = traitlets.Instance(ase.Atoms, allow_none=True)
    
    def __init__(self, **kwargs):
        self._marketplace = MarketPlace()

        self.query_edit = ipw.Text(placeholder='OPTIMADE query')
        self.query_edit.on_submit(self._on_submit)
        self.callback_indicator = ipw.Valid(description='osp-callback')

        super().__init__(children=[self.query_edit, self.callback_indicator], **kwargs)
        
    @staticmethod
    def _to_ase(cuds_object):
        from osp.wrappers.simase import SimaseSession
        from osp.core.namespaces import simase_ontology

        session = SimaseSession()
        ase_wrapper = simase_ontology.simase_wrapper(session=session)
        ase_wrapper.add(cuds_object)
        session.run()
        return session._atoms

    def _get_from_optimade_gateway(self, query):
        response = self._marketplace.get(f'/api/proxy/proxy/{OPTIMADE_GATEWAY_APP}/{query}')
        response.raise_for_status()
        cuds_object = deserialize(response.json())
        self.structure = self._to_ase(cuds_object)

    def _on_submit(self, query_edit):
        self._get_from_optimade_gateway(query_edit.value)

    def from_callback(self, url):
        with self.hold_trait_notifications():
            try:
                query = parse_qs(urlsplit(url).query)
                optimade_query = query['optimade_query'][0]
            except (KeyError, IndexError):
                self.callback_indicator.value = False
            else:
                self._get_from_optimade_gateway(optimade_query)
                self.callback_indicator.value = True
