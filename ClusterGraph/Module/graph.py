import copy
import pickle
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
        self.gene_dict = {}

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



#---------------FIND CLUSTER -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    def find_cluster(self, gene_id):
        try:
            return self.gene_dict[gene_id]
        except:
            raise KeyError('''Le gene entre n'existe pas.''')

#---------------------LOAD GRAPH----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
    def cytoscape(self,file='/home/saiant01/PycharmProjects/Git/Cytoscape/cytoscape.txt'):
        print('Cytoscape ...')
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
        print('Cytoscape done.')

#------------------------GRAPH JAVACRIPT------------------------------------------------------------------------------------------------------------------------------
    def graph_javascript(self):
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

        print('sous-graph done ...')
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
                                    hello=((moins.rstrip(moins.split('_')[-1]).rstrip('_')))
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
        f1=open(grep_file, 'r')
        f2 = open('/home/saiant01/Desktop/gene_compa graph = Graph()re.xml', 'w')
        line=f1.readline()
        f2.writelines('<?xml version="1.0" encoding="UTF-8"?>'+'\n')
        f2.writelines('<ComparingSequences ComparingSequences="{}">'.format('')+'\n')
        dict_of_tuples={}
        while line!='':
            #if int(line.split()[0])>=hit_number:
                #gene_id = (line.split()[1])
                gene_id = (line.split()[0])
                cluster_number =  self.find_cluster(gene_id)
                list_of_paths = self.find_path(cluster_number, path_lenght)
                sequence_path = self.sequences_in_find_path(list_of_paths)
                for sequences,  list_of_paths in sequence_path.items():
                    if sequences not in dict_of_tuples.keys():
                        dict_of_tuples[(gene_id,sequences)]=list_of_paths
                line = f1.readline()
            #elif int(line.split()[0])<hit_number:
                #line=f1.readline()
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

#------------------------Grep fred-----------------------------------------------------------------------------------------------------------------------------------------
    def fred(self,file):
        f1=open(file,'r')
        f2=open('/home/saiant01/Desktop/compare_samples.txt','w')

        line = f1.readline()
        while line!='':
            f2.writelines(line.split('\t')[0]+'_'+line.split('\t')[1]+'\t'+line.split('\t')[6]+'\n')
            line=f1.readline()

        f1.close()
        f2.close()

#--------------------------Save graph--------------------------------------------------------------------------------------------------------------------------------------
    def save_graph(self,file):
        f1=open(file,'wb')
        pickle.dump(self.nodes,f1)
        pickle.dump(self.echantillons,f1)
        pickle.dump(self.gene_dict, f1)

#---------------------------Load graph-------------------------------------------------------------------------------------------------------------------------------------
    def reload_graph(self,file):
        f1 = open(file,'rb')
        nodes = pickle.load(f1)
        échantillons = pickle.load(f1)
        gene_dict = pickle.load(f1)
        self.nodes = nodes
        self.echantillons = échantillons
        self.gene_dict = gene_dict
        f1.close()



#-----------------------MAIN------------------------------------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print('\033[4m'+'start'+ '\033[0m'+' ...')
    print('')
    #Load graph
    graph = Graph()


    #graph.load_graph(('/home/saiant01/Desktop/cat_prodigal_cd_hit_P4.fasta.clstr'))

    #P4J0 vs P4J7
    #graph.load_graph(('/home/saiant01/Desktop/cat_prodigal_cd-hit_p0p7.fasta.clstr'))

    #Local
    #graph.load_graph(('/home/saiant01/Desktop/cat_Sample_P4Jx-Assembly_cd-hit.fa.clstr'))

    #Server
    #graph.load_graph(('/home/saiant01/cat_Sample_PxJy-Assembly_cd-hit.fasta.clstr'))
    #print('Nombre de noeux:'+str(len(graph.nodes)))

    #DataTest
    graph.load_graph(('/home/saiant01/PycharmProjects/Git/Data/cat_prodigal-cd-hit.fasta.clstr'))

    graph.save_graph('/home/saiant01/PycharmProjects/Git/Data/hi.txt')
    graph = Graph()
    graph.reload_graph('/home/saiant01/PycharmProjects/Git/Data/hi.txt')

    #Trouver les clusters pour lequel le gene appartient
    #print(graph.find_cluster('Sample_P4J7-FOX-ANA-Assembly.fa_contig-3000007_12'))


    #Liste des chemins partant du cluster en faisant au maximum n pas.
    #list_of_paths=graph.find_path('Cluster 0', 10)

    #Coloring
    #sequence_path = graph.sequences_in_find_path(list_of_paths)
    #graph.show_path_by_samples(sequence_path)


    #Sous-graph
    #sous_graph = graph.sous_graph(list_of_paths)

    #Visualisation
    graph.graph_javascript()
    #graph.cytoscape()
    #sous_graph.graph_javascript()
    #sous_graph.cytoscape()


    #Compare
    #graph.fred('/home/saiant01/Desktop/Beta_Lactam_Fred_P4J0-7-90')
    #(graph.compare_sequences_excel('/home/saiant01/Desktop/compare_samples.txt',9))



    print('Time:',graph.function_time(time), '/  H:M:S')



