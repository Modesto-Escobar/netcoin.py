import os, shutil, re, json

www = os.path.join(os.path.dirname(__file__), 'www')

PROJECT_MARK = '<!--netCoin Project-->'
_FONT_STYLES = ("styles.css", "styles2.css")
_FONTS = ("Roboto-Regular-latin.woff2", "Roboto-Regular-latin-ext.woff2")


def _is_netcoin_page(html):
  """True for a page written by this package, including those of earlier versions.

  Earlier versions did not write the project mark, but their pages always held
  the JSON data and loaded the engine scripts from ``scripts/``.
  """
  if PROJECT_MARK in html:
    return True
  return ('id="data"' in html and 'src="scripts/' in html
          and any(('scripts/' + js) in html for js in ('network.js', 'barplot.js', 'timeline.js',
                                                       'multigraph.js', 'polygraph.js', 'pie.js')))


def createHTML(directory, styles, scripts, data=None, title=None):
  """Write ``directory/index.html`` with its styles, scripts and the JSON ``data``.

  Like rD3plot, an existing directory is only replaced when it holds a page made
  by this package (its index.html carries the netCoin mark); an empty one is
  used as it is, and anything else raises ``FileExistsError`` instead of being
  deleted.
  """
  directory = str(directory)
  index_path = os.path.join(directory, "index.html")
  if os.path.exists(directory):
    if os.path.isfile(index_path):
      with open(index_path, encoding="utf-8", errors="replace") as f:
        if not _is_netcoin_page(f.read()):
          raise FileExistsError("directory: '%s' already exists and is not a netCoin graph "
                                "(delete it or choose another one)" % directory)
      shutil.rmtree(directory)
    elif os.path.isdir(directory) and not os.listdir(directory):
      pass
    else:
      raise FileExistsError("directory: '%s' already exists" % directory)
  os.makedirs(directory, exist_ok=True)

  with open(os.path.join(www, 'template.html'), encoding="utf-8") as f:
    html = f.read()
  name = title if title is not None else os.path.basename(os.path.normpath(directory))
  html = html.replace('titulo', name)

  head = ''
  if styles:
    os.makedirs(os.path.join(directory, 'styles'), exist_ok=True)
  for style in styles or []:
    head += '\n<link rel="stylesheet" type="text/css" href="styles/' + style + '"></link>'
    shutil.copy(os.path.join(www, style), os.path.join(directory, 'styles'))
    if style in _FONT_STYLES:
      for font in _FONTS:
        shutil.copy(os.path.join(www, font), os.path.join(directory, 'styles'))
  if scripts:
    os.makedirs(os.path.join(directory, 'scripts'), exist_ok=True)
  for script in scripts or []:
    head += '\n<script src="scripts/' + script + '"></script>'
    shutil.copy(os.path.join(www, script), os.path.join(directory, 'scripts'))
  html = html.replace('<!--scripts-->', '<!--scripts-->' + head)

  if data is not None:
    if callable(data):
      data = data()
    # '</' is escaped so that the data can never close the script element
    payload = json.dumps(data, ensure_ascii=False).replace('</', '<\\/')
    html = html.replace('<!--json-->', '<script type="application/json" id="data">' + payload + '</script>')

  with open(index_path, "w", encoding="utf-8") as f:
    f.write(html)
  print("The graph has been generated in the '%s' folder." % directory)


def isnumeric(x):
  if type(x) is list:
    return all(map(lambda z : isnumeric(z),x))
  else:
    return isinstance(x, (int, float, complex)) and not isinstance(x, bool)

def getLanguageScript(obj):
  """The language file ("en.js", "es.js" or "ca.js") for a netCoin object or a net dict."""
  options = obj.get('options') if isinstance(obj, dict) else getattr(obj, 'options', None)
  if isinstance(options, dict):
    language = options.get("language")
    return (language if language in ("es", "ca") else "en") + ".js"
  return None

def transpose(m):
  return [[m[j][i] for j in range(len(m))] for i in range(len(m[0]))]



# ------------------------------------------------------------- languages ----
languages = ['en', 'es', 'ca']

nameList = dict(zip(languages, ['name', 'nombre', 'nom']))
frequencyList = ["frequency", "frecuencia", "freqüència", "%"]
frequencyNames = dict(zip(languages, frequencyList[:3]))
labelList = dict(zip(languages, ['label', 'etiqueta', 'etiqueta']))


def getByLanguage(varlist, language):
  """The entry of ``varlist`` for ``language`` ('en' if it is unknown)."""
  if language not in varlist:
    language = "en"
  return varlist[language]


groupList = dict(zip(languages, ['Group', 'Grupo', 'Grup']))
groupsList = dict(zip(languages, ['Groups', 'Grupos', 'Grups']))
casesList = dict(zip(languages, ['Cases', 'Casos', 'Casos']))
restList = dict(zip(languages, ['rest', 'resto', 'resta']))
