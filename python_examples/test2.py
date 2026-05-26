# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:31:00 2026

@author: Modesto
"""

import pandas as pd

from netcoin import dichotomize as di
from netcoin import simEdge
from netcoin import netcoin



df = pd.DataFrame({
    "texto": ["me llamo;isabel", "hola,mundo", "hola.isabel", None],
    "texto2": ["hola", "hola qué tal", "hola me llamo jose", "quieres"]
})

incidencia=di(df, variables="texto2")
aa = netcoin.coin(incidencia)
di(df, "texto2")
bb = simEdge.sim(aa)
simEdge.edge_list(aa)
cc =simEdge.distant(bb, True)
print(cc)
