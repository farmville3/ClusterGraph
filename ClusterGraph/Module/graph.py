import copy
import os
import time
import datetime
from collections import defaultdict
#----------------------Imports-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Début de la procédure. On initialise le temps.
start_time = time.time()

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

#--------------------- Classer en ordre les gènes dans les différents échantillons pour pouvoir retrouvé facilement les gènes avant et apres un gène en particulier------------------------------------------
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
    @staticmethod
    def function_time(_start_time):
        end_time = time.time() - start_time
        return(datetime.timedelta(seconds=end_time))

#-----------------------------Méthode qui nous perment d'ajouter rapidement un noeux a un graph.-------------------------------------------------------------------------------------------------------------
    def append(self, node):
        if isinstance(node,Node):
            #Verifier si node.cluster_id n'est pas vide
            if node.cluster_id:
                self.nodes[node.cluster_id]=node
        else:
            raise TypeError('''L'objet que vous essayé d'ajouter dans le graph n'est pas un objet de type noeux!''')


#---------------------LOAD GRAPH----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def load_graph(self, file):
        #Créer un graph à partir du fichier output cd-hit
        f1 = open(file, 'r')
        line = f1.readline()
        gene_dict = {}
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
                        gene_dict[gene] = cluster_name
                    else:
                        break
        f1.close()

        # Création de la liste contenant tous les identifiants des échantillons
        sample_id_list = []
        for gene_ids in gene_dict.keys():
            sample_id = (gene_ids.rstrip(gene_ids.split('_')[-1])).rstrip('_')
            if sample_id not in sample_id_list:
                sample_id_list.append(sample_id)

        # Ajouter les gènes à la liste pour laquel le bon identifiant de échantillon est attaché
        for sample_ids in sample_id_list:
            self.echantillons[sample_ids] = []
            for gene_ids in gene_dict.keys():
                if (gene_ids).startswith(sample_ids):
                    (self.echantillons.get(sample_ids)).append(gene_ids)
                else:
                    pass

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
                            if ((self.echantillons[sample_id]).index(genes_in_list))==len((self.echantillons[sample_id]))-1:
                                previous_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list))-1]
                                previous_gene_cluster = gene_dict[previous_gene]
                                if previous_gene_cluster not in cluster.links:
                                    (cluster.links).append(previous_gene_cluster)
                            elif ((self.echantillons[sample_id]).index(genes_in_list))==0:
                                next_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list))+1]
                                next_gene_cluster = gene_dict[next_gene]
                                if next_gene_cluster not in cluster.links:
                                    (cluster.links).append(next_gene_cluster)
                            else:
                                previous_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list)) - 1]
                                next_gene = (self.echantillons[sample_id])[((self.echantillons[sample_id]).index(genes_in_list)) + 1]
                                previous_gene_cluster = gene_dict[previous_gene]
                                next_gene_cluster = gene_dict[next_gene]
                                if previous_gene_cluster not in cluster.links:
                                    (cluster.links).append(previous_gene_cluster)
                                if next_gene_cluster not in cluster.links:
                                    (cluster.links).append(next_gene_cluster)
                    self.append(cluster)
        del i
        f1.close()
#---------------------------GRAPH CYTOSCAPE---------------------------------------------------------------------------------------------------------------------------
    def cytoscape(self,file='/home/saiant01/PycharmProjects/Git/Cytoscape/cytoscape.txt'):
        f2 = open(file, 'w')
        # Cytoscape
        #Créer un fichier tabulé pour cytoscape
        for nodes in self.nodes.values():
            for links in nodes.links:
                links_object = self.nodes[links]
                for genes in nodes.gene_list:
                    plus=genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1])+1)
                    moins=genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1])-1)
                    if plus in links_object.gene_list or moins in links_object.gene_list:
                        f2.writelines((str(nodes.cluster_id) + '\t' + (genes.rstrip(genes.split('_')[-1])).rstrip('_')+ '\t' + str(links_object.cluster_id)+'\n'))
        f2.close()
#------------------------GRAPH JAVACRIPT------------------------------------------------------------------------------------------------------------------------------
    def graph_javascript(self):
        #index.html
        #Créer une liste cluster_start_end pour laquelle chaque élément aura la forme: Cluster_de_départ_du_edge,Cluster_de_fin_du_edge.
        cluster_start_end = []
        dir = os.path.dirname(__file__)
        file = os.path.join(dir, '../../javascript/hcls-dataset-description-master/type-graphs-html/graph.js')
        f3 = open(file, 'w')
        for nodes in self.nodes.values():
            for links in nodes.links:
                links_object = self.nodes[links]
                for genes in nodes.gene_list:
                    plus = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) + 1)
                    moins = genes.rstrip(genes.split('_')[-1]) + str(int(genes.split('_')[-1]) - 1)
                    if plus in links_object.gene_list or moins in links_object.gene_list:
                        cluster_start_end.append(str(nodes.cluster_id) + '\t' + (genes.rstrip(genes.split('_')[-1])).rstrip('_') + '\t' + str(links_object.cluster_id) + '\n')

            #Edges Weight
            weights = []
            for triplets in cluster_start_end:
                #La liste weights est la liste de tous les edges du graph incluant les duplicates.
                weights.append([triplets.split('\t')[0], triplets.split('\t')[2].rstrip('\n')])



            #write nodes in f3
            f3.writelines('var nodes_array = [' + '\n')
            for nodes in self.nodes.values():
                max=0
                for links in nodes.links:
                    if weights.count([nodes.cluster_id,links]) > max:
                        max = weights.count([nodes.cluster_id,links])
                f3.writelines('{id: ' + ((nodes.cluster_id).split(' ')[1])+', weight: '+ str(max)+', label: '+ "'" + '{}'.format(nodes.cluster_id) + "'" + '},' + '\n')
            f3.writelines(']' + '\n')

            #On ajoute tous les edges du graph dans list_edges. Cette liste ne contient pas les duplicates edges.
            list_edges = []
            for nodes in self.nodes.values():
                for links in nodes.links:
                    if [nodes.cluster_id, links] not in list_edges and [links, nodes.cluster_id] not in list_edges:
                        list_edges.append([nodes.cluster_id, links])

            #write edges in f3
            f3.writelines('var edges_array = [' + '\n')
            i=1
            for tuples in list_edges:
                f3.writelines('{id: ' + str(i) + ', from: ' + str((tuples[0]).split(' ')[1]) + ', to: ' + str((tuples[1]).split(' ')[1])+', weight: '+ str(weights.count(tuples)) + ', label: ' + "''"+'},' + '\n')
                i+=1
            f3.writelines(']' + '\n')
        f3.close()



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
                    cluster_id_node = self.nodes[lists[-1]]
                    for links in cluster_id_node.links:
                        if len(cluster_id_node.links)==1 and (cluster_id_node.links)[0] in lists:
                            list_copy = copy.copy(lists)
                            list_of_future_paths.append(list_copy)
                        elif links not in lists:
                            list_copy = copy.copy(lists)
                            list_copy.append(links)
                            list_of_future_paths.append(list_copy)
                        else:
                            pass

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
                        i=1
                        if plus in (self.nodes[paths[1]]).gene_list:
                            while i<len(paths):
                                cluster_working_on = self.nodes[paths[i]]
                                name = cluster_working_on.cluster_id
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
                                elif (moins in cluster_working_on.gene_list) and (
                                    i != (len(paths) - 1)) and (
                                    (moins.rstrip(plus.split('_')[-1]).rstrip('_')) in test_list):
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

#-----------------------MAIN------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    graph = Graph()

    graph.load_graph(('/home/saiant01/PycharmProjects/ClusterGraph/Data/cat_prodigal-cd-hit.fasta.clstr'))

    list_of_paths=graph.find_path('Cluster 10', 4)

    graph.cytoscape()

    graph.graph_javascript()

    print('')
    print(graph.sequences_in_find_path(list_of_paths))
    print(' ')

    print('Time:',graph.function_time(time), '/  H:M:S')






