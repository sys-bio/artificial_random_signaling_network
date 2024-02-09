# This script was written by Jin Xu and available on Github
# https://github.com/SunnyXu/artificial_random_signaling_network

#import libsbml
#print(libsbml.LIBSBML_VERSION)

from libsbml import *
import os
import math

def getSymbols(kinetic_law):
  """
  Finds the parameters and species names for the
  kinetics law. Exposing this information requires
  a recursive search of the parse tree for the
  kinetics expression.
  :return list-str:
  """
  global cur_depth
  MAX_DEPTH = 20
  cur_depth = 0
  def augment(ast_node, result):
    global cur_depth
    cur_depth += 1
    if cur_depth > MAX_DEPTH:
      pass
    for idx in range(ast_node.getNumChildren()):
      child_node = ast_node.getChild(idx)
      if child_node.getName() is None:
        additions = augment(child_node, result)
        result.extend(additions)
      else:
        if child_node.isFunction():
          additions = augment(child_node, result)
          result.extend(additions)
        else:
          result.append(child_node.getName())
    return result

  ast_node = kinetic_law.getMath()
  if ast_node.getName() is None:
    result = []
  else:
    result = [ast_node.getName()]
  return augment(ast_node, result)

files = os.listdir('.')
files.remove('degree_dist.py')
#files.remove('degree_dist.txt')
#files = ['feedback.xml']
xml_num = len(files)
#print(files)

reaction_num_list = []
species_num_list = []
degree_dist_list = []
xml_real_num = 0
for i in range(xml_num): 

  reader = SBMLReader()
  document = reader.readSBMLFromFile(files[i])
  model = document.getModel()

  reaction_num = model.getNumReactions()
  #print("reaction_num:", reaction_num)

  species_num = model.getNumSpecies()
  #print("sepecies_num:", species_num)

  error_flag = 0
  for j in range(reaction_num):
    reaction = model.getReaction(j)
    kinetic_law = reaction.getKineticLaw()
    try: 
      species_parameter_list = list(dict.fromkeys(getSymbols(kinetic_law)))
    except:
      error_flag = 1
      print("Memory issue")
      break

  if reaction_num == 0:
    error_flag = 1
    print("There are no reactions")

  if species_num == 0:
    error_flag = 1
    print("There are no species.")


  if error_flag == 1:
    print(files[i]) 
  else:
    xml_real_num += 1

    #reaction number involved for each species
    count_react_num = [0]*species_num

    species_list = []
    for j in range(species_num):
      species = model.getSpecies(j)
      species_id = species.getId()
      species_list.append(species_id)

    #print(species_list)

    for j in range(reaction_num):
      species_in_kinetic_law = []
      reactant_list = []
      product_list = []
      species_parameter_list = []
      reaction = model.getReaction(j)
      reactant_num = reaction.getNumReactants()
      product_num = reaction.getNumProducts()
      for k in range(reactant_num):
        reactant = reaction.getReactant(k).getSpecies()
        reactant_list.append(reactant)
      for k in range(product_num):
        product = reaction.getProduct(k).getSpecies()
        product_list.append(product)
      kinetic_law = reaction.getKineticLaw()
      #print(kinetic_law.getFormula())
      species_parameter_list = list(dict.fromkeys(getSymbols(kinetic_law)))

      for k in range(len(species_parameter_list)):
        if species_parameter_list[k] in species_list:
          species_in_kinetic_law.append(species_parameter_list[k])
      species_in_reaction = reactant_list + product_list + species_in_kinetic_law 
      species_in_reaction = list(dict.fromkeys(species_in_reaction))
      for k in range(species_num):
        if species_list[k] in species_in_reaction:
          count_react_num[k] += 1

    #order follows species_list 
    #print(count_react_num)    

    degree_dist = [0]*(reaction_num+1)

    for j in range(len(degree_dist)):
      for k in range(species_num):
        if count_react_num[k] == j:
          degree_dist[j] += 1

    #print(degree_dist)
    for j in range(len(degree_dist)):
      degree_dist[j] = degree_dist[j]/species_num

    #0, 1, 2, 3, ..., reaction_num (maximum degree)
    degree_dist_list.append(degree_dist)
    reaction_num_list.append(reaction_num)
    species_num_list.append(species_num)

#print(reaction_num_list)
#print(species_num_list)
if xml_real_num == 0:
  print("There are no valid sbml files.")

else:
  reaction_num_avg = 0
  species_num_avg  = 0
  for i in range(xml_real_num):
    reaction_num_avg += reaction_num_list[i]
    species_num_avg  += species_num_list[i]

  reaction_num_avg = reaction_num_avg/xml_real_num
  species_num_avg  = species_num_avg/xml_real_num

  reaction_num_sdv = 0
  species_num_sdv = 0
  for i in range(xml_real_num):
    reaction_num_sdv += (reaction_num_list[i]-reaction_num_avg)**2
    species_num_sdv  += (species_num_list[i]-species_num_avg)**2

  reaction_num_sdv = math.sqrt(reaction_num_sdv/xml_real_num)/math.sqrt(xml_real_num)
  species_num_sdv = math.sqrt(species_num_sdv/xml_real_num)/math.sqrt(xml_real_num)

  print("valid sbml files:", xml_real_num)
  print("reaction_num:", int(reaction_num_avg), int(reaction_num_sdv))
  print("species_num:", int(species_num_avg), int(species_num_sdv))

  longest_list_len = 0
  for i in range(len(degree_dist_list)):
    if len(degree_dist_list[i]) > longest_list_len:
      longest_list_len = len(degree_dist_list[i])


  degree_dist_list_avg = [0]*longest_list_len
  degree_dist_list_sdv = [0]*longest_list_len


  for i in range(longest_list_len):
    for j in range(len(degree_dist_list)):
      if len(degree_dist_list[j]) > i:
        degree_dist_list_avg[i] += degree_dist_list[j][i]

  degree_dist_list_avg = [x / len(degree_dist_list) for x in degree_dist_list_avg]
  #print(degree_dist_list_avg)

  for i in range(longest_list_len):
    for j in range(len(degree_dist_list)):
      if len(degree_dist_list[j]) > i:
        degree_dist_list_sdv[i] += (degree_dist_list[j][i]-degree_dist_list_avg[i])**2

  degree_dist_list_sdv = [x / len(degree_dist_list) for x in degree_dist_list_sdv] 
  degree_dist_list_sdv = [math.sqrt(x) for x in degree_dist_list_sdv]
  degree_dist_list_sdv = [x/math.sqrt(xml_real_num) for x in degree_dist_list_sdv]
  #print(degree_dist_list_sdv)


  with open('degree_dist.txt', 'w+') as file:
      for listitem in degree_dist_list_avg:
          file.write('%f\t' % listitem)
      file.write('\n')
      for listitem in degree_dist_list_sdv:
          file.write('%f\t' % listitem)
      file.write('\n')

  file.close()
