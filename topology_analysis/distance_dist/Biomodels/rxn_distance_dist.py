# This script was written by Jin Xu and available on Github
# https://github.com/SunnyXu/artificial_random_signaling_network


#import libsbml
#print(libsbml.LIBSBML_VERSION)

from libsbml import *
from collections import OrderedDict
import math
import matplotlib.pyplot as plt
import zipfile

def getSymbols(kinetic_law):
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

def common_specs(a, b):
    a_set = set(a)
    b_set = set(b)

    if (a_set & b_set):
        return len(a_set & b_set)
    else:
        return 0

#files = os.listdir('.')
#files.remove('distance_dist.py')
#files.remove('distance_dist.txt')
#files = ['feedback.xml']

    
with zipfile.ZipFile('Biomodels.zip') as zip_file:
    zip_file.extractall()

with zipfile.ZipFile('Biomodels.zip') as zip_file:
  file_list = zip_file.namelist()
  file_list.remove('Biomodels/')
  xml_num = len(file_list)
  files = file_list

distance_dist_list = []
xml_real_num = 0
for i in range(xml_num):
  incomplete_network_flag = 0
  reader = SBMLReader()
  #print("files:", files[i])
  document = reader.readSBMLFromFile(files[i])
  model = document.getModel()

  species_num = model.getNumSpecies()

  species_list = []
  for j in range(species_num):
    species = model.getSpecies(j)
    species_id = species.getId()
    species_list.append(species_id)


  reaction_num = model.getNumReactions()
  species_in_reaction = [[] for j in range(reaction_num)]

  error_flag = 0
  for j in range(reaction_num):
    reaction = model.getReaction(j)
    kinetic_law = reaction.getKineticLaw()
    #print(kinetic_law.getFormula())
    try: 
      species_parameter_list = list(dict.fromkeys(getSymbols(kinetic_law)))
    except:
      error_flag = 1
      print("Memory Issue")
      break

  if reaction_num == 0:
    error_flag = 1
    print("There are no reactions")

  if species_num == 0:
    error_flag = 1
    print("There are no species")

  if error_flag == 1:
    print(files[i])
  else:
    xml_real_num += 1

    for j in range(reaction_num):
      species_in_kinetic_law = []
      reactant_list = []
      product_list = []
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
      species_in_reaction[j] = reactant_list + product_list + species_in_kinetic_law 
      species_in_reaction[j] = list(dict.fromkeys(species_in_reaction[j]))

    #print("reaction")
    #for j in range(reaction_num):
    #  print(species_in_reaction[j])

    species_list = []
    for j in range(reaction_num):
      species_list.extend(species_in_reaction[j])

    species_list = list(dict.fromkeys(species_list))
    #print("species_list:", species_list)
    species_num = len(species_list)

    pair_size = int(species_num*(species_num-1)/2)

    species_pair = [[] for j in range(pair_size)]
    #print("species_pair:", species_pair)
    count = 0
    for j in range(species_num):
      for k in range(species_num):
        if j<k:
          species_pair[count] = [species_list[j], species_list[k]]
          count += 1

    #print("pair")
    #for j in range(pair_size):
    #  print(species_pair[j])


    distance = [-1]*pair_size
    #for j in range(pair_size):
    # print(j)
    # print(species_pair[j])

    for j in range(pair_size):
      for k in range(reaction_num):
        if species_pair[j][0] in species_in_reaction[k]:
          if species_pair[j][1] in species_in_reaction[k]:
            distance[j] = 0
            
      if distance[j] != 0:
        species_twin1 = []
        species_twin2 = []
        for k in range(reaction_num):
          if species_pair[j][0] in species_in_reaction[k]:
            for l in range(len(species_in_reaction[k])):
              species_twin1.append(species_in_reaction[k][l]) 

        for k in range(reaction_num):
          if species_pair[j][1] in species_in_reaction[k]:
            for l in range(len(species_in_reaction[k])):
              species_twin2.append(species_in_reaction[k][l])

        if common_specs(species_twin1, species_twin2) != 0:
          distance[j] = 1

      if distance[j] != 0 and distance[j] != 1:
        reaction_involved = []
        species_twin1 = []
        species_twin2 = []
        distance[j] = 1
        for k in range(reaction_num):
          if species_pair[j][0] in species_in_reaction[k]:
            reaction_involved.append(k)
            species_twin1.extend(species_in_reaction[k])
        species_twin1 = list(OrderedDict.fromkeys(species_twin1))

        species = []
        for k in range(reaction_num):
          if species_pair[j][1] in species_in_reaction[k]:
            reaction_involved.append(k)
            species.extend(species_in_reaction[k])
        species_twin2 = list(OrderedDict.fromkeys(species))
        reaction_involved_len = len(list(OrderedDict.fromkeys(reaction_involved)))
        common_check = common_specs(species_twin1, species_twin2)

        while common_check == 0 and reaction_involved_len < reaction_num and distance[j] <= reaction_num:
          #species = []
          for k in range(reaction_num):
            for l in range(len(species_twin2)):
              if species_twin2[l] in species_in_reaction[k]:
                reaction_involved.append(k)
                species.extend(species_in_reaction[k])
          species_twin2 = list(OrderedDict.fromkeys(species))
          reaction_involved_len = len(list(OrderedDict.fromkeys(reaction_involved)))
          common_check = common_specs(species_twin1, species_twin2)
          distance[j] += 1
        if common_check == 0 and reaction_involved_len == reaction_num:
          incomplete_network_flag = 1
        if distance [j] >= reaction_num:
          incomplete_network_flag = 1
      
      #print(distantce[j])
    if incomplete_network_flag != 1:
      distant_dist = [0]*reaction_num
      for j in range(reaction_num):
        for k in range(pair_size):
          if distance[k] == j:
            distant_dist[j] += 1 

      #print(distant_dist)
      for j in range(reaction_num):
      #from 1 to reaction_num:
        distant_dist[j] = distant_dist[j]/pair_size
        #print(distant_dist[i])

      #print(distant_dist)
      distance_dist_list.append(distant_dist)
    else:
      xml_real_num -= 1
      print("Incomplete network")
      print(files[i])

if xml_real_num == 0:
  print("There are no valid sbml files.")

else:
  print("valid sbml files:", xml_real_num)
  longest_list_len = 0
  for i in range(len(distance_dist_list)):
    if len(distance_dist_list[i]) > longest_list_len:
      longest_list_len = len(distance_dist_list[i])


  distance_dist_list_avg = [0]*longest_list_len
  distance_dist_list_sdv = [0]*longest_list_len


  for i in range(longest_list_len):
    for j in range(len(distance_dist_list)):
      if len(distance_dist_list[j]) > i:
        distance_dist_list_avg[i] += distance_dist_list[j][i]

  distance_dist_list_avg = [x / len(distance_dist_list) for x in distance_dist_list_avg]
  #print(distance_dist_list_avg)

  for i in range(longest_list_len):
    for j in range(len(distance_dist_list)):
      if len(distance_dist_list[j]) > i:
        distance_dist_list_sdv[i] += (distance_dist_list[j][i]-distance_dist_list_avg[i])**2

  distance_dist_list_sdv = [x / len(distance_dist_list) for x in distance_dist_list_sdv]
  distance_dist_list_sdv = [math.sqrt(x) for x in distance_dist_list_sdv]
  distance_dist_list_sdv = [x/math.sqrt(xml_real_num) for x in distance_dist_list_sdv]
  #print(distance_dist_list_sdv)


  with open('distance_dist.txt', 'w+') as file:
      for listitem in distance_dist_list_avg:
          file.write('%f\t' % listitem)
      file.write('\n')
      for listitem in distance_dist_list_sdv:
          file.write('%f\t' % listitem)
      file.write('\n')

  file.close()

  #dist_list_len = len(degree_dist_list_avg)
  dist_list_len = 11

  x = []
  for i in range(dist_list_len):
    x.append(i)

  y = distance_dist_list_avg[:dist_list_len]
  e = distance_dist_list_sdv[:dist_list_len]

  fig, ax = plt.subplots(figsize = (7, 5))
  ax.errorbar(x, y, yerr=e, fmt="o", c='b', marker="o", label='Biomodels')
  plt.legend(loc='upper right')

  ax.set_xlabel('Reaction distance')
  ax.set_ylabel('Probability')
  ax.set_title('Reaction distance distribution')

  plt.savefig("Distance_Distribution.pdf", format="pdf")
  plt.show()



