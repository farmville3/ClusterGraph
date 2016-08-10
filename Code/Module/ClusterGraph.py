import copy
import pickle
import os
import time
import datetime
from collections import defaultdict
import argparse
import sys
#----------------------Imports-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Début de la procédure. On initialise le temps.
start_time = time.time()


#----------------------OptionParser----------------------------------------------------------------------------------------------------------------------------------------
class OptionParser():
    def __init__(self, args):
        self.parser = argparse.ArgumentParser(description="Code.")


        ##########         parser       #########################
        self.parser.add_argument('-i', type=str, help='Input file.Le fichier en input doit correspondre au fichier .clstr de cd-hit.',default=None,required=True)

        self.parser.add_argument('-find', type=str, help="L'input doit être le nom d'un gène qui est présent dans le graph. La fonction va alors retourner le cluster respectif du gène.", default=None)

        self.parser.add_argument('-lop', type=str, help='La fonction retourne tous les chemins de longueur "x" partant du nom du cluster que vous devez donner en parametètre.'
                                                        ' Utilisez -x pour donner la longueur des chemins désirés.', default=None)

        self.parser.add_argument('-x', type=int, help="Longueur des chemin désiré lors de l'utilisation de -lop ou -mn ou les deux.", default=-1)

        self.parser.add_argument('-sg', type=bool, help="Construit un sous-graph à partir des données reçues par -lop et -x"
                                 '|  True ou False.', default=False)

        self.parser.add_argument('-cyto', type=bool, help="Permet de visualiser le graph dans cytoscape à l'aide du fichier cytoscape.txt"
                                 'Celui-ci doit se trouvé dans le répertoire Cytoscape du projet.'
                                                          '|  True ou False. ', default=False)

        self.parser.add_argument('-j', type=bool, help="Permet de visualiser le graph dans une page html de façon interactive."
                                 "Il ne faut qu'ouvrir le fichier index.html dans javascript/hcls-dataset-description-master/type-graphs-html."
                                                       "|  True ou False.", default=False)

        self.parser.add_argument('-g', type=str, help="Convertit les fichiers obtenus à l'aide de greps en ligne de commande de façon"
                                                      "à ce que les données soient convertible dans un fichier .xml", default=None)

        self.parser.add_argument('-xml', type=bool, help="Créer un fichier xml qui peut être ouvert dans Microsoft Excel."
                                                         "|  True ou False.", default=False)

        self.parser.add_argument('-y', type=int, help="Longueur des chemin désiré lors de l'utilisation de -xml.", default=-1)

        self.parser.add_argument('-s', type=str, help="Affiche différentes statistiques sur le graph.", default=False )

        self.parser.add_argument('-save', type=str, help="Sauvegarde le graph dans le fichier donné en paramètre.", default=None)
        self.parser.add_argument('-r', type=str, help="Load le graph à partir du fichier de sauvegarde deonné en paramètre", default=None)
        self.parser.add_argument('-mn', type=str, help="Créer des fichiers cytoscape pour plusieurs gènes à la fois."
                                 'Utilisez -x quand vous utilisez -mn pour spécifier la longueur des chemins.', default=None)

        #
        #parse the args...
        self.Arguments = vars(self.parser.parse_args(args))

    #Needed to get the args
    def getArguments(self):
        return self.Arguments

    def getParser(self):
        return self.parser


#---------------------------------NODE-------------------------------------------------------------------------------------------------------------------------------------
#Le constructeur Node construit des objets de type Node qui ont comme attribut un identifiant de cluster, une liste de genes liés à ce cluster et une liste
class Node:
    def __init__(self,cluster_id, gene_list=None, links=None):
        self.cluster_id = cluster_id
        if gene_list is None:
            self.gene_list=[]
        else:
            self.gene_list = gene_list
        if links is None:
            self.links=[]
        else:
            self.links = links

#-------------------GRAPH---------------------------------------------------------------------------------------------------------------------------------------------------
#Un graph possède comme attribut plusieurs noeux.
class Graph:
    def __init__(self):
        self.nodes = {}
        self.echantillons = defaultdict(list)
        self.gene_dict = {}

#--------------------- Order by asc-----------------------------------------------------------------------------------------------------------------------------------------
#Classer en ordre les gènes dans les différents échantillons pour pouvoir retrouvé facilement les gènes avant et apres un gène en particulier
    def order_by_asc(self):
        for sample_id, gene in self.echantillons.items():
            i = 0
            while i < len(self.echantillons[sample_id]):
                ((self.echantillons[sample_id])[i]) = int((((self.echantillons[sample_id])[i]).split('_')[-1]))
                i=i+1
            self.echantillons[sample_id].sort()
            i=0
            while i < len(self.echantillons[sample_id]):
                ((self.echantillons[sample_id])[i]) = str(sample_id) + '_' + str((self.echantillons[sample_id])[i])
                i = i + 1
#-----------------------TIME------------------------------------------------------------------------------------------------------------------------------------------------
    #Calcule le temps passé entre le début du programme et sa fin
    @staticmethod
    def function_time(_start_time):
        print('--------------------------------------------------------------------------------------------------------')
        print("\t" + "\t" + "\t" + "\t" + '\t'+'\t' + 'Time')
        print('--------------------------------------------------------------------------------------------------------')
        end_time = time.time() - start_time
        return(datetime.timedelta(seconds=end_time))

#-----------------------------Append in graph-------------------------------------------------------------------------------------------------------------
#Méthode simple qui nous perment d'ajouter rapidement un noeux a un graph.
    def append(self, node):
        if isinstance(node,Node):
            #Verifier si node.cluster_id n'est pas vide
            if node.cluster_id:
                self.nodes[node.cluster_id]=node
        else:
            raise TypeError('''L'objet que vous essayé d'ajouter dans le graph n'est pas un objet de type noeux!''')



#---------------FIND CLUSTER -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Méthode simple qui permet de récupérer le cluster respectif d'un gène.
    def find_cluster(self, gene_id):
        try:
            return self.gene_dict[gene_id]
        except:
            raise KeyError('''Le gene entre n'existe pas.''')

#---------------------LOAD GRAPH----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Créer un graph selon le fichier .clstr de cd-hit
    def load_graph(self, file):
        #Créer un graph à partir du fichier output cd-hit
        f1 = open(file, 'r')
        line = f1.readline()
        line_list=[]
        # Remplir le dictionnaire gene_dict qui a la forme Gene#:Cluster#
        while line != '':
            if line.startswith('>Cluster'):
                cluster_name = line.lstrip('>').rstrip('\n')
                line_list.append(line)
                line = f1.readline()
                while line.startswith('>Cluster') == False:
                    if line !='':
                        gene = (line.split('>')[1].split('...')[0])
                        line_list.append(line)
                        line = f1.readline()
                        self.gene_dict[gene] = cluster_name
                    else:
                        break
        f1.close()

        # Ajouter les gènes à la liste pour laquel le bon identifiant de échantillon est attaché
        for gene_ids in self.gene_dict.keys():
            sample_id = (gene_ids.rstrip(gene_ids.split('_')[-1])).rstrip('_')
            self.echantillons[sample_id].append(gene_ids)

        # Classer en ordre les gènes dans les différents échantillons

        self.order_by_asc()

        #Pour tous les objets nodes, on leur attribut un Cluster# et une liste de gène qui sont attachés à ce cluster#.

        i=0
        while i<len(line_list):
            if line_list[i].startswith('>Cluster'):
                cluster_name = (line_list[i]).lstrip('>').rstrip('\n')
                i+=1
                cluster = Node(cluster_name)
                while i<len(line_list) and ((line_list[i].startswith('>Cluster')) == False):
                        gene = ((line_list[i]).split('>')[1].split('...')[0])
                        if gene not in cluster.gene_list:
                            cluster.gene_list.append(gene)
                        i += 1
                if ((len(cluster.gene_list)!=0) or i<len(line_list) and line_list[i].startswith('>Cluster') ):
                #Remplir la liste des links pour chaque node
                    for genes_in_list in cluster.gene_list:
                        sample_id = (genes_in_list.rstrip(genes_in_list.split('_')[-1])).rstrip('_')
                        if genes_in_list in self.echantillons[sample_id]:
                            #Liste de longueur 1
                            if len(self.echantillons[sample_id])==1:
                                continue
                            #Le gene est le dernier gene dans la liste des genes dans le dict echantillon
                            elif ((self.echantillons[sample_id]).index(genes_in_list))==len((self.echantillons[sample_id]))-1:
                                previous_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list))-1]
                                previous_gene_cluster = self.gene_dict[previous_gene]
                                if previous_gene_cluster not in cluster.links:
                                    (cluster.links).append(previous_gene_cluster)
                            #Le gene est le premier gene dans la liste des genes dans le dict echantillon
                            elif ((self.echantillons[sample_id]).index(genes_in_list))==0:
                                next_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list))+1]
                                next_gene_cluster = self.gene_dict[next_gene]
                                if next_gene_cluster not in cluster.links:
                                    (cluster.links).append(next_gene_cluster)
                            else:
                                previous_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list)) - 1]
                                next_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list)) + 1]
                                previous_gene_cluster = self.gene_dict[previous_gene]
                                next_gene_cluster = self.gene_dict[next_gene]
                                if previous_gene_cluster not in cluster.links:
                                    (cluster.links).append(previous_gene_cluster)
                                if next_gene_cluster not in cluster.links:
                                    (cluster.links).append(next_gene_cluster)
                    self.append(cluster)
        del i
        f1.close()
#---------------------------GRAPH CYTOSCAPE---------------------------------------------------------------------------------------------------------------------------
    #Écrit un fichier qui permet la visualisation du graph avec cytoscape
    def cytoscape(self, file=(os.path.join(os.path.dirname(__file__), '../../Cytoscape/cytoscape.txt'))):
        print('Cytoscape ...')
        f2 = open(file, 'w')
        dict={}
        # Cytoscape
        #Créer un fichier tabulé pour cytoscape
        for nodes in self.nodes.values():
            for links in nodes.links:
                links_object = self.nodes[links]
                for genes in nodes.gene_list:
                    plus=genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1])+1)
                    moins=genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1])-1)
                    if plus in links_object.gene_list or moins in links_object.gene_list:
                        debut=str(nodes.cluster_id)
                        fin=str(links_object.cluster_id)
                        lien=(genes.rstrip(genes.split('_')[-1])).rstrip('_')
                        if (fin,lien,debut) not in dict.keys():
                            dict[(debut,lien,fin)]=''

        for path in dict.keys():
                        f2.writelines((path[0]) + '\t' + path[1] + '\t' + path[2] +'\n')
        f2.close()
        print('Cytoscape done.')

#------------------------GRAPH JAVACRIPT------------------------------------------------------------------------------------------------------------------------------
    def graph_javascript(self):
        # Écrit un fichier qui permet la visualisation du graph avec un script javascript
        print('Javascript ...')
        #index.html
        #Créer une liste cluster_start_end pour laquelle chaque élément aura la forme: Cluster_de_départ_du_edge,Cluster_de_fin_du_edge.
        dir = os.path.dirname(__file__)
        file = os.path.join(dir, '../../javascript/hcls-dataset-description-master/type-graphs-html/graph.js')
        f3 = open(file, 'w')
        weights = defaultdict(int)
        for nodes in self.nodes.values():
            #node_object = self.nodes[nodes]
            for links in nodes.links:
                links_object = self.nodes[links]
                for genes in nodes.gene_list:
                    plus = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) + 1)
                    moins = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) - 1)
                    if plus in links_object.gene_list or moins in links_object.gene_list:
        # La liste weights est la liste de tous les edges du graph incluant les duplicates.
                        weights[(str(nodes.cluster_id),str(links_object.cluster_id))] += 1


        #write nodes in f3
        f3.writelines('var nodes_array = [' + '\n')
        dict_edges = {}
        counter = 0
        buffer = ''
        for noeux in self.nodes.values():
            maximum=0
            for links in noeux.links:
                maximum = max(weights[(noeux.cluster_id,links)],maximum)
                # On ajoute tous les edges du graph dans dict_edges. Ce dictionnaire ne contient pas les duplicates edges.
                if (noeux.cluster_id, links) not in dict_edges and (links, noeux.cluster_id) not in dict_edges:
                    dict_edges[(noeux.cluster_id, links)]=''

            buffer+=('{id: ' + ((noeux.cluster_id).split(' ')[1])+', weight: '+ str(maximum)+', label: '+ "'" + '{}'.format(noeux.cluster_id) + "'" + '},' + '\n')
            counter+=1
            if counter%50000==0:
                f3.writelines(buffer)
                buffer = ''
        f3.writelines(buffer)
        f3.writelines(']' + '\n')



        #write edges in f3
        f3.writelines('var edges_array = [' + '\n')
        i=0
        buffer=''
        for tuples in dict_edges:
            buffer+=('{id: ' + str(i) + ', from: ' + str((tuples[0]).split(' ')[1]) + ', to: ' + str((tuples[1]).split(' ')[1])+', weight: '+ str(weights[tuples]) + ', label: ' + "''"+'},' + '\n')
            i += 1
            if i%50000==0:
                f3.writelines(buffer)
                buffer=''
        f3.writelines(buffer)
        f3.writelines(']' + '\n')

        f3.close()
        del weights
        print('Javascript done.')

#-------------------SOUS_GRAPH---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def sous_graph(self,paths_list):
    #Permet de créer un sous graph selon un cluster_id et un path_lenght
        print('sous-graph ...')
        graph = Graph()
        for paths in paths_list:
            i=0
            while i<len(paths):
                if i == len(paths)-1:
                    cluster_depart = Node(paths[-1])
                    cluster_precedent = paths[i - 1]
                    if cluster_depart.cluster_id in graph.nodes:
                        ((graph.nodes[cluster_depart.cluster_id]).links).append(cluster_precedent)
                    else:
                        (cluster_depart.links).append(cluster_precedent)
                        graph.append(cluster_depart)
                    break
                else:
                    if i==0:
                        cluster_depart = Node(paths[i])
                        cluster_arrive = paths[i + 1]
                        if cluster_depart.cluster_id in graph.nodes:
                            ((graph.nodes[cluster_depart.cluster_id]).links).append(cluster_arrive)
                        else:
                            (cluster_depart.links).append(cluster_arrive)
                            graph.append(cluster_depart)
                        i+=1
                    else:
                        cluster_precedent = paths[i - 1]
                        cluster_depart = Node(paths[i])
                        cluster_arrive = paths[i + 1]
                        if cluster_depart.cluster_id in graph.nodes:
                            ((graph.nodes[cluster_depart.cluster_id]).links).append(cluster_arrive)
                            ((graph.nodes[cluster_depart.cluster_id]).links).append(cluster_precedent)
                        else:
                            (cluster_depart.links).append(cluster_arrive)
                            (cluster_depart.links).append(cluster_precedent)
                            graph.append(cluster_depart)
                        i += 1


        for noeux in graph.nodes.values():
            try:
                noeux.gene_list = (self.nodes[noeux.cluster_id]).gene_list

            except:
                raise ValueError('''Le noeux n'existe pas''')

        print('sous-graph done.')
        return graph

#-----------------------FIND PATH--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #Trouver tous les chemins de longueur 'path_lenght' en partant du noeux ayant comme identifiant de cluster, le cluster_id donné en parametre.
    def find_path(self,cluster_id,path_lenght):

        list_of_current_paths=[]
        list_of_future_paths=[]
        try:
            cluster_id_node = (self.nodes)[cluster_id]
            for links in cluster_id_node.links:
                list_of_current_paths.append([cluster_id_node.cluster_id,links])
            i=1
            while i < path_lenght:
                for lists in list_of_current_paths:
                    #Si les noeux forment un cercle
                    for links in self.nodes[lists[-1]].links:
                        if links in lists:
                            verification='yes'
                        else:
                            verification='no'
                            break
                    if verification=='no':
                        for links in self.nodes[lists[-1]].links:
                            list_of_links=self.nodes[lists[-1]].links
                            #Si le lien n'est pas dans la lists
                            if links not in lists:
                                list_copy = copy.copy(lists)
                                list_copy.append(links)
                                if list_copy not in list_of_future_paths:
                                    list_of_future_paths.append(list_copy)
                            #Si nos clusters et nos edges forment un cercle
                            elif links in lists:
                                pass

                    elif verification=='yes':
                        list_copy = copy.copy(lists)
                        if list_copy not in list_of_future_paths:
                            list_of_future_paths.append(list_copy)
                list_of_current_paths=list_of_future_paths
                list_of_future_paths=[]
                i+=1
            del i
            return list_of_current_paths

        except:
            raise TypeError('''Le cluster entré en paramètre n'existe pas.''')
#---------------------------------------------------------------------------------------------------------------------------------------
    #Retourne les séquences passants sur les paths retournés par find_path()
    def sequences_in_find_path(self,list_of_paths):
        sequence_paths = {}
        if isinstance(self,Graph):
            if isinstance(list_of_paths, list):
                for paths in list_of_paths:
                    return_list = []
                    test_list=[]
                    node_object = self.nodes[paths[0]]
                    for genes in node_object.gene_list:
                        plus = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) + 1)
                        moins = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) - 1)
                        if plus in (self.nodes[paths[1]]).gene_list:
                            i = 1
                            while i<len(paths):
                                cluster_working_on = self.nodes[paths[i]]
                                #dernier cluster
                                if i==(len(paths)-1) and (plus in cluster_working_on.gene_list) and ((plus.rstrip(plus.split('_')[-1]).rstrip('_')) in test_list):
                                    return_list.append((plus.rstrip(plus.split('_')[-1]).rstrip('_')))
                                    break
                                #premier cluster
                                elif i == 1 and (plus in cluster_working_on.gene_list):
                                    test_list.append((plus.rstrip(plus.split('_')[-1]).rstrip('_')))
                                    i += 1
                                    plus = plus.rstrip(plus.split('_')[-1]) + str(int(plus.split('_')[-1]) + 1)
                                #autres clusters
                                elif (plus in cluster_working_on.gene_list) and (i!=(len(paths)-1)) and ((plus.rstrip(plus.split('_')[-1]).rstrip('_')) in test_list):
                                    i += 1
                                    plus = plus.rstrip(plus.split('_')[-1]) + str(int(plus.split('_')[-1]) + 1)
                                else:
                                    break
                        elif moins in (self.nodes[paths[1]]).gene_list:
                            i = 1
                            while i < len(paths):
                                cluster_working_on = self.nodes[paths[i]]
                                # dernier cluster
                                if i == (len(paths) - 1) and (moins in cluster_working_on.gene_list) and ((moins.rstrip(moins.split('_')[-1]).rstrip('_')) in test_list):
                                    return_list.append((moins.rstrip(moins.split('_')[-1]).rstrip('_')))
                                    break
                                # premier cluster
                                elif i == 1 and (moins in cluster_working_on.gene_list):
                                    test_list.append((moins.rstrip(moins.split('_')[-1]).rstrip('_')))
                                    i += 1
                                    moins = moins.rstrip(moins.split('_')[-1]) + str(int(moins.split('_')[-1]) - 1)
                                # autres clusters
                                elif (moins in cluster_working_on.gene_list) and (i != (len(paths) - 1)) and ((moins.rstrip(moins.split('_')[-1]).rstrip('_')) in test_list):
                                    i += 1
                                    moins = moins.rstrip(moins.split('_')[-1]) + str(int(moins.split('_')[-1]) - 1)
                                else:
                                    break

                    #Retourner le tout dans un dictionnaire
                    for samples in return_list:
                        if not return_list:
                            pass
                        else:
                            if samples not in sequence_paths.keys():
                                sequence_paths[samples] = [paths]
                            else:
                                sequence_paths[samples].append(paths)

            else:
                raise TypeError('''La liste que vous avez tentez d'entrer n'est pas un objet de type list.''')
        else:
            raise TypeError('''Le graph que vous avez tentez d'entrer n'est pas un objet de type graph.''')


        return sequence_paths

# -------------------------Print paths by samples------------------------------------------------------------------------------------------------------------------------------------
    def show_path_by_samples(self, sequences_paths):
    #Imprime les paths pour chaque séquence
        if isinstance(sequences_paths,dict):
            for sequences in sequences_paths.keys():
                for cluster_list in sequences_paths.values():
                    print('\033[1m'+'Sequences:'+'\033[0m', sequences)

                    print('\033[1m'+'Paths:'+'\033[0m')
                    for paths in cluster_list:
                        print(paths)
                    print('')
        else:
            raise TypeError('''L'objet entré en paramètre n'est pas un objet de type dictionnaire''')

#------------------------COMPARE GRAPHS------------------------------------------------------------------------------------------------------------------------------------
    def compare_sequences_excel(self,grep_file, path_lenght):
        #Créer un fichier .xml qui s'ouvre facilement avec excel qui affiche les chemins de certains gène en paticulier
        f1=open(grep_file, 'r')
        f2 = open('/home/saiant01/PycharmProjects/Git/Data/XML/gene_compare.xml', 'w')
        line=f1.readline()
        f2.writelines('<?xml version="1.0" encoding="UTF-8"?>'+'\n')
        f2.writelines('<ComparingSequences ComparingSequences="{}">'.format('')+'\n')
        dict_of_tuples={}
        while line!='':
            gene_id = (line.split()[0])
            cluster_number =  self.find_cluster(gene_id)
            list_of_paths = self.find_path(cluster_number, path_lenght)
            sequence_path = self.sequences_in_find_path(list_of_paths)
            for sequences,  list_of_paths in sequence_path.items():
                if sequences not in dict_of_tuples.keys():
                    dict_of_tuples[(gene_id,sequences)]=list_of_paths
            line = f1.readline()

        for tuples, list_of_paths in dict_of_tuples.items():
            if len(list_of_paths) != 0:
                f2.writelines('\t' + '<Sequence name="{}">'.format(str(tuples[1])) + '\n')
                f2.writelines('\t' + '\t'  '<Gene_id GeneNumber="{}">'.format(str(tuples[0].split('_')[-1])) + '\n')
                for paths in list_of_paths:
                    f2.writelines('\t' + '\t' + '\t' '<Paths>' + '\n')
                    i = 0
                    while i < len(paths):
                        f2.writelines(
                            ((4) * '\t') + '<Node{}>'.format(i) + str(paths[i]) + '</Node{}>'.format(i) + '\n')
                        i += 1
                    f2.writelines('\t' + '\t' + '\t' '</Paths>' + '\n')
                f2.writelines('\t' + '\t' '</Gene_id>' + '\n')
                f2.writelines('\t' + '</Sequence>' + '\n')

        f2.writelines('</ComparingSequences>'+'\n')
        f2.close()
        f1.close()

#------------------------Grep fred-----------------------------------------------------------------------------------------------------------------------------------------
    def beta_lactam_file(self,file):
        #Convertis les données reçu grace aux 'greps' de fred
        f1=open(file,'r')
        dir = os.path.dirname(__file__)
        write_file = os.path.join(dir, '/home/saiant01/PycharmProjects/Git/Data/XML/compare_samples.txt')
        f2=open(write_file,'w')
        line = f1.readline()
        try:
            while line!='':
                f2.writelines(line.split('\t')[0]+'_'+line.split('\t')[1]+'\t'+line.split('\t')[6]+'\n')
                line=f1.readline()
        except:
            f1.close()
            f1 = open(file, 'r')
            dir = os.path.dirname(__file__)
            write_file = os.path.join(dir, '/home/saiant01/PycharmProjects/Git/Data/XML/compare_samples.txt')
            f2 = open(write_file, 'w')
            line = f1.readline()
            while line!='':
                f2.writelines(line.split('\t')[0]+'\n')
                line=f1.readline()
        f1.close()
        f2.close()

#--------------------------Save graph--------------------------------------------------------------------------------------------------------------------------------------
    def save_graph(self,file):
        #Enregistre un graph
        f1=open(file,'wb')
        pickle.dump(self.nodes,f1)
        pickle.dump(self.echantillons,f1)
        pickle.dump(self.gene_dict, f1)
        f1.close()

#---------------------------Load graph-------------------------------------------------------------------------------------------------------------------------------------
    def reload_graph(self,file):
        #Load un graph enregistré
        f1 = open(file,'rb')
        nodes = pickle.load(f1)
        échantillons = pickle.load(f1)
        gene_dict = pickle.load(f1)
        self.nodes = nodes
        self.echantillons = échantillons
        self.gene_dict = gene_dict
        f1.close()


#-------------------------Compare graphs------------------------------------------------------------------------------------------------------------------------------------
    def stats_graphs(self):
        #Imprime certains statistiques de graph
        print('--------------------------------------------------------------------------------------------------------')
        print("\t"+"\t"+"\t"+"\t"+'\t'+'Some stats on the graph')
        print('--------------------------------------------------------------------------------------------------------')
        print(" "+str(len(self.nodes))+ ' noeux')
        max=0
        min=100

        same_sample=[]
        for nodes in self.nodes.values():
            if len(nodes.gene_list)>max:
                max=len(nodes.gene_list)

            if len(nodes.gene_list)<min:
                min=len(nodes.gene_list)
            else:
                pass
            for genes in nodes.gene_list:
                debut = genes.rstrip(genes.split('_')[-1])
                counter = 0
                for genes in nodes.gene_list:
                    if genes.startswith(debut):
                        counter +=1
                    if counter==2:
                        if nodes.cluster_id not in same_sample:
                            same_sample.append(nodes.cluster_id)
        print(" "+'Max: '+ str(max))
        print(" "+'Min: '+ str(min))
        print(" "+str(len(same_sample)))


#----------------Cytoscape for more than one gene/cluster---------------------------------------------------------------------------------------------------------

    def many_genes_cytoscape(self, file, lenght):
        f1 = open(file,'r')
        line = f1.readline()

        while line!='':
            gene = line.rstrip('\n')
            cluster = self.find_cluster(gene)
            list_of_paths = self.find_path(cluster,lenght)
            sous_graph = self.sous_graph(list_of_paths)
            dir = os.path.dirname(__file__)
            write_file = os.path.join(dir, '/home/saiant01/PycharmProjects/ClusterGraph/Cytoscape/many_genes_cytoscape/{}.txt'.format(gene))
            sous_graph.cytoscape(write_file)
            line=f1.readline()



#-----------------------MAIN------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print('========================================================================================================')
    print(
        "\t" + "\t"+ "\t"+ "\t"  + "\t" + "\t" + '\033[4m' + 'start' + '\033[0m')
    print('========================================================================================================')

    parser = OptionParser(sys.argv[1:])
    arg = parser.getArguments()
    if (arg['i'])!=None:
        load_file=(arg['i'])
    else:
        load_file = None
    find = (arg['find'])
    node_for_list_of_paths = (arg['lop'])
    lenght_lop = int(arg['x'])
    sous_graph = bool(arg['sg'])
    cytoscape = bool(arg['cyto'])
    javacript = bool(arg['j'])
    convert_greps_fred = (arg['g'])
    xml = bool(arg['xml'])
    lenght_xml = int(arg['y'])
    many = arg['mn']
    stats=(arg['s'])
    save = (arg['save'])
    reload = (arg['r'])

    # Load graph
    graph = Graph()

    if load_file !=None:
        # graph.load_graph(('/home/saiant01/cat_Sample_P4Jx-Assembly_100.fa.clstr'))

        # P4
        graph.load_graph((load_file))

        if save != None:
            graph.save_graph(save)

        if reload != None:
            graph = Graph()
            graph.reload_graph(save)

        # P4J0 vs P4J7
        # graph.load_graph(('/home/saiant01/Desktop/cat_prodigal_cd-hit_p0p7.fasta.clstr'))

        # Local
        # graph.load_graph(('/home/saiant01/Desktop/cat_Sample_P4Jx-Assembly_cd-hit.fa.clstr'))



        # Server
        # graph.load_graph(('/home/saiant01/cat_Sample_PxJy-Assembly_cd-hit.fasta.clstr'))
        # print('Nombre de noeux:'+str(len(graph.nodes)))

        # DataTest
        # graph.load_graph(('/home/saiant01/PycharmProjects/Git/Data/cat_prodigal-cd-hit.fasta.clstr'))


        if find != None:
            # Trouver les clusters pour lequel le gene appartient
            print(graph.find_cluster(find))
        else:
            pass



            # Liste des chemins partant du cluster en faisant au maximum n pas.
        if lenght_lop != -1 and node_for_list_of_paths != None:
            list_of_paths = graph.find_path(node_for_list_of_paths, lenght_lop)

            # Coloring
            sequence_path = graph.sequences_in_find_path(list_of_paths)
            graph.show_path_by_samples(sequence_path)
        else:
            pass

            # Sous-graph
        if sous_graph != False:
            sous_graph = graph.sous_graph(list_of_paths)
        else:
            pass

            # Visualisation
        if cytoscape == True and sous_graph == True:
            sous_graph.cytoscape()
        elif cytoscape == True:
            graph.cytoscape()
        if javacript == True and sous_graph == True:
            sous_graph.graph_javascript()
        elif cytoscape == True:
            graph.graph_javascript()


            # Compare
        if convert_greps_fred != None:
            file = graph.beta_lactam_file(convert_greps_fred)
        else:
            pass

        if xml != False and lenght_xml != -1:
            (graph.compare_sequences_excel(file, lenght_xml))
        else:
            pass

        print('')
        print(" "+'Done!')
        print('')

        if many!=None:
            graph.many_genes_cytoscape(many,lenght_lop)

            # Stats
        if stats == 'True':
            graph.stats_graphs()
        else:
            pass



            # Time
        print(" "+'Time:', graph.function_time(time), '/  H:M:S')


    else:
        raise FileNotFoundError("Aucun fichier n'a été donné en entré!")



