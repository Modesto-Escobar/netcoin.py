# -*- coding: utf-8 -*-
"""
Created on Mon May  4 12:31:00 2026

@author: Modesto
"""

import pandas as pd

from netcoin import coin, dichotomize, edge_list, sim
from netcoin.similarity import distant


df = pd.DataFrame({
    "texto": ["me llamo;isabel", "hola,mundo", "hola.isabel", None],
    "texto2": ["hola", "hola qué tal", "hola me llamo jose", "quieres"]
})

incidencia = dichotomize(df, variables="texto2", sep=" ")
aa = coin(incidencia)
print(aa)
bb = sim(aa)
print(edge_list(aa, criteria="Jaccard", procedures="Jaccard", min=0))
cc = distant(bb, True)
print(cc)
