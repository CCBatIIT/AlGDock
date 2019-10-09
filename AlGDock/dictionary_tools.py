# merge_dictionaries, convert_dictionary_relpath, and dict_view

import os
import numpy as np

from collections import OrderedDict


def merge_dictionaries(dicts, required_consistency=[]):
  """
  Merges a list of dictionaries, giving priority to items in descending order.
  Items in the required_consistency list must be consistent with one another.
  """
  merged = OrderedDict()
  for a in range(len(dicts)):  # Loop over all dictionaries,
    # giving priority to the first
    if not isinstance(dicts[a], dict):
      continue
    for key in dicts[a].keys():
      if isinstance(dicts[a][key], dict):  # item is a dictionary
        merged[key] = merge_dictionaries(
            [
                dicts[n][key]
                for n in range(len(dicts)) if key in dicts[n].keys()
            ],
            required_consistency=required_consistency)
      elif (key not in merged.keys()):
        # Merged dictionary will contain value from
        # first dictionary where key appears
        merged[key] = dicts[a][key]
        # Check for consistency with other dictionaries
        for b in (range(a) + range(a + 1, len(dicts))):
          if isinstance(dicts[b],dict) and (key in dicts[b].keys()) \
              and (dicts[a][key] is not None) and (dicts[b][key] is not None):
            if (isinstance(dicts[b][key], np.ndarray)):
              inconsistent_items = (dicts[b][key] != dicts[a][key]).any()
            else:
              inconsistent_items = (dicts[b][key] != dicts[a][key])
            if inconsistent_items:
              if key in required_consistency:
                print 'Dictionary items for %s are inconsistent:' % key
                print dicts[a][key]
                print dicts[b][key]
                raise Exception('Items must be consistent!')
      elif (merged[key] is None):  # Replace None
        merged[key] = dicts[a][key]
  return merged


def convert_dictionary_relpath(d, relpath_o=None, relpath_n=None):
  """
  Converts all file names in a dictionary from one relative path to another.
  If relpath_o is None, nothing is joined to the original path.
  If relpath_n is None, an absolute path is used.
  """
  converted = OrderedDict()
  for key in d.keys():
    if d[key] is None:
      pass
    elif isinstance(d[key], dict):
      converted[key] = convert_dictionary_relpath(d[key],
                                                  relpath_o=relpath_o,
                                                  relpath_n=relpath_n)
    elif isinstance(d[key], str):
      if d[key] == 'default':
        converted[key] = 'default'
      else:
        if relpath_o is not None:
          p = os.path.abspath(os.path.join(relpath_o, d[key]))
        else:
          p = os.path.abspath(d[key])
        # if os.path.exists(p): # Only save file names for existent paths
        if relpath_n is not None:
          converted[key] = os.path.relpath(p, relpath_n)
        else:
          converted[key] = p
  return converted


def dict_view(dict_c, indent=2, relpath=None, show_None=False):
  view_string = ''
  for key in dict_c.keys():
    if dict_c[key] is None:
      if show_None:
        view_string += ' ' * indent + key + ': None\n'
    elif isinstance(dict_c[key], dict):
      subdict_string = dict_view(dict_c[key], indent + 2, relpath=relpath)
      if subdict_string != '':
        view_string += ' ' * indent + key + ':\n' + subdict_string
    elif isinstance(dict_c[key], str):
      view_string += ' ' * indent + key + ': '
      if relpath is not None:
        view_string += os.path.relpath(dict_c[key], relpath) + '\n'
        if not os.path.exists(dict_c[key]):
          view_string += ' ' * (indent + len(key)) + 'DOES NOT EXIST\n'
      else:
        view_string += dict_c[key] + '\n'
    else:
      view_string += ' ' * indent + key + ': ' + repr(dict_c[key]) + '\n'
  return view_string
