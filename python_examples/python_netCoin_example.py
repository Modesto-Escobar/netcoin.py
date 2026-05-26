from netcoin import netCoin, multigraphCreate
import csv, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))


l = ['Source','Target','Haberman','Z']
x = [[6.0, 8.0, 14.38, 5.08e-37],
[8.0, 9.0, 14.11, 5.19e-36],
[5.0, 8.0, 13.48, 1.35e-33],
[5.0, 9.0, 13.40, 2.70e-33],
[6.0, 9.0, 12.80, 4.75e-31],
[7.0, 8.0, 12.53, 4.87e-30],
[5.0, 6.0, 12.43, 1.19e-29],
[7.0, 9.0, 12.00, 4.47e-28],
[4.0, 8.0, 11.92, 8.56e-28],
[5.0, 7.0, 11.90, 9.93e-28],
[6.0, 7.0, 11.36, 8.87e-26],
[2.0, 8.0, 11.22, 2.93e-25],
[4.0, 7.0, 10.85, 5.65e-24],
[4.0, 9.0, 10.65, 2.96e-23],
[2.0, 6.0, 10.47, 1.21e-22],
[2.0, 7.0, 10.21, 9.00e-22],
[4.0, 6.0, 9.80, 2.22e-20],
[2.0, 5.0, 9.68, 5.65e-20],
[4.0, 5.0, 8.91, 1.77e-17],
[2.0, 9.0, 8.49, 3.41e-16],
[2.0, 4.0, 7.49, 3.13e-13],
[3.0, 9.0, 7.01, 6.41e-12],
[9.0, 10.0, 6.23, 6.98e-10],
[3.0, 7.0, 6.20, 8.31e-10],
[2.0, 3.0, 5.88, 4.79e-09],
[3.0, 5.0, 5.87, 5.08e-09],
[3.0, 8.0, 5.34, 8.64e-08],
[5.0, 10.0, 5.08, 3.03e-07],
[6.0, 10.0, 4.79, 1.21e-06],
[3.0, 4.0, 4.55, 3.72e-06],
[8.0, 10.0, 4.54, 3.94e-06],
[7.0, 10.0, 4.15, 2.05e-05],
[4.0, 10.0, 3.72, 1.15e-04],
[3.0, 6.0, 3.55, 2.16e-04],
[2.0, 10.0, 2.47, 6.88e-03],
[1.0, 6.0, 2.11, 1.77e-02],
[1.0, 5.0, 1.91, 2.81e-02],
[1.0, 8.0, 1.72, 4.24e-02],
[1.0, 9.0, 1.23, 0.10],
[1.0, 7.0, 1.06, 0.14],
[1.0, 3.0, 0.96, 0.16],
[1.0, 2.0, 0.86, 0.19],
[3.0, 10.0, 0.76, 0.22]]

listaVar = ['variable', 'Sexo', 'Edad', 'Ambito']
nodes = {
  'variable': ['Unamuno', 'Lizarraga', 'Fernando', 'Pablo', 'Salome', 'Felisa', 'Jose', 'Maria', 'Rafael', 'Ramon', 'Jugo', 'MariaH', 'MiguelN', 'SalomeN', 'CarminaN', 'MercedesN', 'MiguelinN', 'ConchaN', 'MTeresaN', 'Molina', 'Mugica', 'Soriano', 'TJ', 'Onis', 'Landa', 'Legendre', 'Castinera', 'Ortega', 'Perez', 'Velarde', 'Bastianini'],
  'Sexo': ['Hombre', 'Mujer', 'Hombre', 'Hombre', 'Mujer', 'Mujer', 'Hombre', 'Mujer', 'Hombre', 'Hombre', 'Mujer', 'Mujer', 'Hombre', 'Mujer', 'Mujer', 'Mujer', 'Hombre', 'Mujer', 'Mujer', 'Mujer', 'Hombre', 'Hombre', 'Mujer', 'Hombre', 'Hombre', 'Hombre', 'Hombre', 'Hombre', 'Hombre', 'Hombre', 'Mujer'],
  'Edad': ['Contemporáneos', 'Contemporáneos', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Descendientes-1', 'Ascendientes', 'Contemporáneos', 'Descendientes-2', 'Descendientes-2', 'Descendientes-2', 'Descendientes-2', 'Descendientes-2', 'Descendientes-2', 'Descendientes-2', 'Contemporáneos', 'Contemporáneos', 'Contemporáneos', 'Histórico', 'Contemporáneos', 'Ascendientes', 'Contemporáneos', 'Contemporáneos', 'Contemporáneos', 'Contemporáneos', 'Contemporáneos', 'Descendientes-1'],
  'Ambito': ['Personaje', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Familia', 'Amistad', 'Amistad', 'Amistad', 'Otros', 'Amistad', 'Familia', 'Amistad', 'Amistad', 'Amistad', 'Amistad', 'Amistad', 'Amistad']
}

net = netCoin.fromMatrix(x,nodes,l)
net.netCoinArgs(nodeAttrNames=listaVar,linkAttrNames=l,name='variable',color="Edad",shape="Sexo",main='Fotos de Unamuno',layout="kamada_kawai",background="./icons/usal.png",language="es")
net.plot(directory="tmp/unamuno")
print(net)


# Gender
N = {
  'Gender':["Man","Women","Undet."],
  'type':['male','female','none'],
  'frequency':[3,3,1],
  'icons':["./icons/male.png","./icons/female.png","./icons/undet.png"],
  'icons2':["./icons/bioinfo.png","./icons/bioinfo.png","./icons/nucleus.png"],
  'fx':[1,3,2],
  'fy':[1,1,2]
}

E = {
  'value':[1,2],
  'Haberman':[0.66,0.66],
  'p(Z)':[0.27,0.27]
}

net2 = netCoin([(0,2), (1,2)],vertex_attrs=N, edge_attrs=E, name = "Gender", lwidth="value", defaultColor = "firebrick", background="#e8e3d0", main = "This is a title", limits = [0,0,4,3], image = ["icons","icons2"], directory = "tmp/gender")


# Dice
ncolnames = ['name','frequency']
nodes = {
  'name':['1','2','3','4','5','6','odd','even','small','large'],
  'frequency':[15, 13, 26, 18, 13, 15, 54, 46, 54, 46]
}

lcolnames = ['Source','Target','Haberman','Z']
x = [
  [1,7,3.88,9.46e-05],
  [1,9,3.88,9.46e-05],
  [2,8,4.19,3.03e-05],
  [2,9,3.57,2.77e-04],
  [3,7,5.47,1.66e-07],
  [3,9,5.47,1.66e-07],
  [4,8,5.08,8.91e-07],
  [4,10,5.08,8.92e-07],
  [5,7,3.57,2.77e-04],
  [5,10,4.19,3.03e-05],
  [6,8,4.55,7.51e-06],
  [6,10,4.55,7.51e-06],
  [7,9,4.77,3.19e-06],
  [8,10,4.77,3.19e-06]
]
tree = [['odd','1'],['even','2'],['odd','3'],['even','4'],['odd','5'],['even','6']]

net3 = netCoin.fromMatrix(x,nodes,lcolnames)
net3.netCoinArgs(nodeAttrNames=ncolnames,linkAttrNames=lcolnames,tree=tree,lwidth="Haberman")
net3.summary()
net3.plot(directory="tmp/dice")

nets = []
edges = [list(map(lambda x: x-1,row[:2])) for row in x]
nets.append(netCoin(edges, vertex_attrs=nodes, tree=tree))
for i in range(7,15):
  nodes['name'].append(str(i))
  nodes['frequency'].append(0)
  edges.append([nodes['name'].index(str(i-1)),nodes['name'].index(str(i))])
  nets.append(netCoin(edges, vertex_attrs=nodes, tree=tree))

multigraphCreate(nets, mode="f", directory="tmp/dice_incremental")


#italian
filepath = "./italian/"

families = csv.DictReader(open(filepath+"families.tsv"),delimiter="\t")
marriage = csv.DictReader(open(filepath+"marriage.tsv"),delimiter="\t")
netM = netCoin.DictList(vertices=families, edges=marriage, edge_foreign_keys=("Source", "Target"))
netM.netCoinArgs(layout="mds",size="f.Marriages",color="seat", main="Marriage Links beetween Italian families", note="Data source: Padgett & Ansell (1983)")

families = csv.DictReader(open(filepath+"families.tsv"),delimiter="\t")
business = csv.DictReader(open(filepath+"business.tsv"),delimiter="\t")
netB = netCoin.DictList(vertices=families, edges=business, edge_foreign_keys=("Source", "Target"))
netB.netCoinArgs(layout="mds",size="f.Business",color="seat", main="Business Links beetween Italian families", note="Data source: Padgett & Ansell (1983)")

multigraphCreate([netM,netB], names=['Marriages','Business'], directory="tmp/multgraphItalian")
multigraphCreate([netM,netB], mode="p", directory="tmp/splitScreenItalian")
multigraphCreate(["tmp/multgraphItalian","tmp/splitScreenItalian"], directory="tmp/multi_iframes")


#sanderson
filepath = "./sanderson/"

finches = csv.reader(open(filepath+"finches.tsv"),delimiter="\t")
ncolnames = next(finches)
finches = [row for row in finches]
finches = zip(*finches)
finches = dict(zip(ncolnames,finches))
finches['species'] = list(map(lambda x: "./sanderson/images/"+x,finches['species']))

galapagos = csv.reader(open(filepath+"galapagos.tsv"),delimiter="\t")
lcolnames = next(galapagos)
galapagos = [row for row in galapagos]

sanderson = netCoin.TupleList(edges=galapagos, edge_attrs=lcolnames[2:])
for col in finches:
  sanderson.vs[col] = finches[col]

sanderson.netCoinArgs(nodeAttrNames=ncolnames,linkAttrNames=lcolnames,lwidth="Haberman",lweight="Haberman", size="frequency", image="species", layout="mds", main="Species coincidences in Galapagos Islands", note="Data source: Sanderson (2000)")
sanderson.plot()


#incremental
nets = []

N = {'name': ['Man','Woman','Undet0']}
E = [(0,1),(1,2)]

nets.append(netCoin(E, vertex_attrs=N, label="", repulsion=100))

for i in range(100):
  N['name'].append('Undet'+str(i+1))
  E.append((i+2,i+3))
  nets.append(netCoin(E, vertex_attrs=N, label="", repulsion=99-i))

multigraphCreate(nets, mode="f", directory="tmp/incremental")


#arte
filepath = "./Arte/"

nodes = list(csv.DictReader(open(filepath+"Arte_nodes.tsv"),delimiter="\t"))
links = list(csv.DictReader(open(filepath+"Arte_links.tsv"),delimiter="\t"))
art = netCoin.DictList(vertices=nodes, edges=links, vertex_name_attr="artista", edge_foreign_keys=("Source", "Target"))
art.netCoinArgs(name="artista",lwidth="Haberman",layout="ka",size="importance")
art.plot()

fechas = {}
for i in range(len(nodes)):
  artista = nodes[i]['artista']
  birth = int(nodes[i]['birth'])
  fechas[artista] = birth

nets = []
for i in range(127,191):
  n = i*10
  sublinks = []
  subnodes = []
  for j in range(len(nodes)):
    if fechas[nodes[j]['artista']]<=n:
      subnodes.append(nodes[j])
  for j in range(len(links)):
    if fechas[links[j]['Source']]<=n and fechas[links[j]['Target']]<=n:
      sublinks.append(links[j])
  if len(sublinks) and len(subnodes):
    net = netCoin.DictList(vertices=subnodes, edges=sublinks, vertex_name_attr="artista", edge_foreign_keys=("Source", "Target"))
    net.netCoinArgs(name="artista",lwidth="Haberman")
    nets.append(net)

multigraphCreate(nets, mode="f", directory="tmp/Arte_time")


#Adjacency matrix
adj = [[0,1,0,0,1,0],[1,0,1,0,1,0],[0,1,0,1,0,0],[0,0,1,0,1,1],[1,1,0,1,0,0],[0,0,0,1,0,0]]
net = netCoin.Adjacency(adj)
net.plot()


#graphml
net = netCoin.Read_GraphML("simple.graphml")
net.netCoinArgs(name="id")
net.plot()


#votarem
filepath = "./Proces/"
nets = []
for i in range(80):
  nodes = csv.DictReader(open(filepath+"S"+str(i+1)+"_nodes.tsv"),delimiter="\t")
  links = csv.DictReader(open(filepath+"S"+str(i+1)+"_links.tsv"),delimiter="\t")
  net = netCoin.DictList(vertices=nodes, edges=links, edge_foreign_keys=("Source", "Target"))
  net.netCoinArgs(color="Walktrap",size="tuits")
  nets.append(net)

multigraphCreate(nets, mode="f", directory="tmp/votarem")

