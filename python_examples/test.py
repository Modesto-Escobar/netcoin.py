import netcoin as nc
#from netcoin import dichotomize
#from netcoin import sim
#from netcoin import coin
#from netcoin import edge_list

Nodes = {
  'Gender':["Man","Women","Undet."],
  'frequency':[3,3,1],
  'fx':[1,3,2],
  'fy':[1,1,2]
}

Edges = {
  'Source':["Man","Women"],
  'Target':["Undet.","Undet."],
  'value':[1,2],
  'Haberman':[0.66,0.66],
  'p(Z)':[0.27,0.27]
}

net = nc.netCoin(Nodes, Edges, name = "Gender", lwidth="value", defaultColor = "firebrick", main = "This is a title", directory = "gender", language='es')

net.plot()

#dichotomize test

import pandas as pd

df = pd.DataFrame({
    "texto": ["me llamo;isabel", "hola,mundo", "hola.isabel", None],
    "texto2": ["hola", "hola qué tal", "hola me llamo jose", "quieres"]
})

incidencia=nc.dichotomize(df, variables="texto2")

aa = nc.coin(incidencia)

bb = nc.sim(aa)

nc.edge_list(aa, procedures=["j", "d"])

nc.sim(aa, procedures="h")
